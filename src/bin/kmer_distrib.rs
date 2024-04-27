use clap::Parser;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Parser, Clone, Debug)]
#[clap(
    version,
    about = "Evaluates genome read distribution and estimates reads per species for specific taxonomy IDs.",
    long_about = "Analyzes each genome's read distribution and calculates the expected number of reads for each species that correspond to a given taxonomy ID."
)]
pub struct Args {
    /// Kraken counts file for each genome mapped to the overall database.
    #[clap(short, long, required = true)]
    input: PathBuf,

    /// Output file containing each classified taxonomy ID and the
    /// kmer distributions of all genomes with this classification.
    #[clap(short, long, required = true)]
    output: PathBuf,
}

fn parse_single_genome(curr_str: &str) -> (String, u32, HashMap<String, u32>) {
    let split_str: Vec<&str> = curr_str.trim().split('\t').collect();
    if split_str.len() < 4 {
        return ("0".to_string(), 0, HashMap::new());
    }

    let genome_taxid = split_str[1].to_string();
    let mut mapped_id_kmers = HashMap::new();
    let mut total_kmers = 0;

    for kmers in split_str[3].split_whitespace() {
        let pair: Vec<&str> = kmers.split(':').collect();
        if pair.len() != 2 {
            continue;
        }
        let (curr_m_id, curr_kmers) = (pair[0], pair[1].parse::<u32>().unwrap_or(0));
        total_kmers += curr_kmers;
        *mapped_id_kmers.entry(curr_m_id.to_string()).or_insert(0) += curr_kmers;
    }

    if mapped_id_kmers.is_empty() {
        return ("0".to_string(), 0, HashMap::new());
    }

    (genome_taxid, total_kmers, mapped_id_kmers)
}

pub fn run(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    let input_path = args.input;
    let file = File::open(input_path)?;
    let reader = BufReader::new(file);

    let mut genome_dict: HashMap<String, HashMap<String, u32>> = HashMap::new();
    let mut genome_dict_totalkmers: HashMap<String, u32> = HashMap::new();
    let mut num_genomes = 0u32;

    for line in reader.lines() {
        let line = line?;
        let (genome_taxid, total_kmers, mapped_taxids_kmers) = parse_single_genome(&line);

        if genome_taxid == "0" {
            continue;
        }

        let counter = genome_dict_totalkmers
            .entry(genome_taxid.clone())
            .or_insert(0);
        *counter += total_kmers;
        if !genome_dict.contains_key(&genome_taxid) {
            num_genomes += 1;
        }

        let sub_map = genome_dict
            .entry(genome_taxid.clone())
            .or_insert(HashMap::new());
        for (m_taxid, count) in mapped_taxids_kmers {
            *sub_map.entry(m_taxid).or_insert(0) += count;
        }
    }
    println!(
        "...{} total genomes read from kraken output file",
        num_genomes
    );

    let mut mapped_taxids_dict: HashMap<String, HashMap<String, u32>> = HashMap::new();
    for (genome, sub_map) in &genome_dict {
        for (m_taxid, count) in sub_map {
            mapped_taxids_dict
                .entry(m_taxid.clone())
                .or_insert_with(HashMap::new)
                .insert(genome.clone(), *count);
        }
    }

    let mut output_file = File::create(args.output)?;
    writeln!(
        output_file,
        "mapped_taxid\tgenome_taxids:kmers_mapped:total_genome_kmers"
    )?;

    for (m_taxid, sub_map) in &mapped_taxids_dict {
        let mut line = format!("{}\t", m_taxid); // 以基因组 ID 开头

        for (genome_taxid, count) in sub_map {
            let total_kmers = genome_dict_totalkmers.get(genome_taxid).unwrap_or(&0);
            line.push_str(&format!("{}:{}:{} ", genome_taxid, count, total_kmers));
            // 构建同一行的多个条目
        }

        writeln!(output_file, "{}", line.trim_end())?; // 写入整行，去除末尾的空格
    }

    Ok(())
}

#[allow(dead_code)]
fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Application error: {}", e);
    }
}
