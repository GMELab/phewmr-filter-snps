use std::{
    collections::{HashMap, HashSet},
    io::Write,
};

use clap::Parser;

#[derive(Parser)]
#[command(version)]
pub struct Cli {
    #[arg(long)]
    outcome: String,
    #[arg(long)]
    eth_outcome: String,
    #[arg(long)]
    outcome_source: String,
    #[arg(long)]
    sumstat_file: String,
    #[arg(long)]
    output_file: String,
    #[arg(long)]
    list_snps_exposures: String,
    #[arg(long)]
    outcome_panel: String,
    #[arg(long)]
    outcome_assay: String,
    #[arg(long)]
    outcome_gene: String,
    #[arg(long)]
    sample_size_outcome: String,
    #[arg(long)]
    n_case_outcome: String,
    #[arg(long)]
    n_control_outcome: String,
}

fn main() {
    let _ = env_logger::Builder::from_env(env_logger::Env::default().filter_or("RUST_LOG", "warn"))
        .try_init();

    let Cli {
        mut outcome,
        eth_outcome,
        outcome_source,
        sumstat_file,
        output_file,
        list_snps_exposures,
        outcome_panel,
        outcome_assay,
        outcome_gene,
        sample_size_outcome,
        n_case_outcome,
        n_control_outcome,
    } = Cli::parse();

    let mut col_num_out = HashMap::new();
    let file = std::fs::File::open(sumstat_file).unwrap();
    let reader = flate2::bufread::GzDecoder::new(std::io::BufReader::new(file));
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(reader);
    let headers_rec = rdr.headers().unwrap();
    let headers = headers_rec
        .into_iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>();
    // drop(headers_rec);
    for (i, header) in headers.into_iter().enumerate() {
        col_num_out.insert(header, i);
    }

    let snps_exposures = std::fs::read_to_string(list_snps_exposures).unwrap();
    let lines = snps_exposures
        .lines()
        .map(|x| x.split('\t'))
        .collect::<Vec<_>>();

    let mut snps_exposures: HashSet<(&str, &str, &str, &str)> = HashSet::new();
    for mut line in lines {
        let chr = line.next().unwrap();
        let pos = line.next().unwrap();
        let ref_ = line.next().unwrap();
        let alt = line.next().unwrap();
        snps_exposures.insert((chr, pos, ref_, alt));
        snps_exposures.insert((chr, pos, alt, ref_));
    }

    let mut formatted_rows = HashSet::new();

    let ref_i = col_num_out.get("ref").unwrap();
    let alt_i = col_num_out.get("alt").unwrap();
    let effect_size = if outcome_source == "PURE_BM" {
        col_num_out.get("beta_raw").unwrap()
    } else {
        col_num_out.get("effect_size").unwrap()
    };
    if outcome_source == "PURE_BM" {
        outcome = format!("{eth_outcome}_{outcome_panel}_{outcome_assay}_{outcome_gene}");
    }
    let chr_i = if outcome_source == "PURE_BM" {
        col_num_out.get("chr").unwrap()
    } else {
        col_num_out.get("chr_hg19").unwrap()
    };
    let pos_i = if outcome_source == "PURE_BM" {
        col_num_out.get("pos").unwrap()
    } else {
        col_num_out.get("pos_hg19").unwrap()
    };
    let se_i = if outcome_source == "PURE_BM" {
        col_num_out.get("se").unwrap()
    } else {
        col_num_out.get("standard_error").unwrap()
    };
    let af_i = if eth_outcome == "EUR" || eth_outcome == "EUROPEAN" {
        if outcome_source == "PURE_BM" {
            col_num_out.get("compiled_EUR_AF").unwrap()
        } else {
            col_num_out.get("gnomAD_AF_EUR").unwrap()
        }
    } else if outcome_source == "PURE_BM" {
        col_num_out.get("study_AF").unwrap()
    } else {
        col_num_out.get("EAF").unwrap()
    };
    let pvalue_i = col_num_out.get("pvalue").unwrap();
    #[derive(PartialEq, Eq)]
    enum N {
        Pure,
        Na,
        SampleSizeOutcome,
    }
    let n = if outcome_source == "PURE_BM" {
        N::Pure
    } else if sample_size_outcome == "NA" && n_case_outcome == "NA" && n_control_outcome == "NA" {
        N::Na
    } else {
        N::SampleSizeOutcome
    };
    let (n_total, n_case, n_ctrl) = if n == N::Na {
        (
            col_num_out.get("N_total").unwrap(),
            col_num_out.get("N_case").unwrap(),
            col_num_out.get("N_ctrl").unwrap(),
        )
    } else {
        (&0, &0, &0)
    };

    for i in rdr.records() {
        let record = i.unwrap();
        let chr = record.get(*chr_i).unwrap();
        let pos = record.get(*pos_i).unwrap();
        let ref_ = record.get(*ref_i).unwrap();
        let alt = record.get(*alt_i).unwrap();
        if snps_exposures.contains(&(chr, pos, ref_, alt)) {
            let effect_size = record.get(*effect_size).unwrap();
            if ref_.len() + alt.len() == 2 && effect_size != "0" {
                let se = record.get(*se_i).unwrap();
                let af = record.get(*af_i).unwrap();
                let pvalue = record.get(*pvalue_i).unwrap();
                let (n_total, n_case, n_ctrl) = if n == N::Na {
                    (
                        record.get(*n_total).unwrap(),
                        record.get(*n_case).unwrap(),
                        record.get(*n_ctrl).unwrap(),
                    )
                } else if n == N::SampleSizeOutcome {
                    (
                        sample_size_outcome.as_str(),
                        n_case_outcome.as_str(),
                        n_control_outcome.as_str(),
                    )
                } else {
                    (sample_size_outcome.as_str(), "NA", "NA")
                };
                formatted_rows.insert(format!("{outcome}\t{chr}_{pos}_{ref_}_{alt}\t{chr}\t{pos}\t{effect_size}\t{se}\t{alt}\t{ref_}\t{af}\t{pvalue}\tNA\tNA\t{n_total}\t{n_case}\t{n_ctrl}"));
            }
        }
    }

    let mut file = std::fs::File::create(output_file).unwrap();
    for row in formatted_rows {
        writeln!(file, "{}", row).unwrap();
    }
}
