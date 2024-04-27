use bracken::kraken;
use chrono::{DateTime, Local};
use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser, Clone, Debug)]
#[clap(
    version,
    about = "Estimates species or genus level abundance from Kraken outputs using Bayesian methods.",
    long_about = "Estimates species or genus level abundance based on assigned reads and expected kmer distributions from Kraken report outputs."
)]
pub struct Args {
    /// Kraken counts file for each genome mapped to the overall database.
    #[clap(short, long, required = true)]
    input: PathBuf,

    /// Kmer distribution file.
    #[clap(short, long, required = true)]
    kmer_distr: PathBuf,

    /// Output modified kraken report file with abundance estimates
    #[clap(short, long, required = true)]
    output: PathBuf,

    /// Level to push all reads to [default: S].
    #[clap(short, long, default_value = "S")]
    level: String,

    /// Threshold for the minimum number of reads kraken must assign
    /// to a classification for that classification to be considered in the
    /// final abundance estimation.
    #[clap(short, long, default_value_t = 10)]
    threshold: usize,
}

fn check_and_parse(input: &str) -> Result<usize, &'static str> {
    let mut chars = input.chars();

    // 检查第一个字符是否是字母
    if let Some(first_char) = chars.next() {
        if !first_char.is_alphabetic() {
            return Err("The first character is not a letter.");
        }
    } else {
        return Err("Input is empty.");
    }

    let remaining: String = chars.collect();

    if remaining.is_empty() {
        return Ok(0);
    }

    // 尝试将剩余的字符串解析为数字
    remaining
        .parse::<usize>()
        .map_err(|_| "Failed to parse the remaining characters as a number.")
}

fn check_report_file(input_file: &PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    eprintln!(">> Checking report file: {:?}", input_file);
    let r_file = File::open(input_file)?;
    let mut reader = BufReader::new(r_file);
    let mut first_line = String::new();

    // 读取第一行
    if reader.read_line(&mut first_line)? == 0 {
        eprintln!("File is empty");
        return Ok(());
    }

    if let Some(first_char) = first_line.chars().next() {
        if first_char == 'C' || first_char == 'U' {
            eprintln!("\tERROR: Bracken does not use the Kraken default output.");
            eprintln!(
                "\t       Bracken requires the Kraken report file (--report option with Kraken)"
            );
            return Err(Box::new(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid file format for Bracken",
            )));
        }
    }

    // 检查是否是 mpa 风格报告
    if first_line.split('\t').count() == 2 {
        eprintln!("\tERROR: Bracken is not compatible with mpa-style reports.");
        eprintln!("\t       Bracken requires the default Kraken report format");
        return Err(Box::new(io::Error::new(
            io::ErrorKind::InvalidData,
            "Invalid file format for Bracken",
        )));
    }

    Ok(())
}

pub fn run(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    let mut lvl_dict: HashMap<String, &str> = HashMap::new();

    let now: DateTime<Local> = Local::now();
    let time = now.format("%m-%d-%Y %H:%M:%S").to_string();
    println!("PROGRAM START TIME: {}", time);

    lvl_dict.insert("D".into(), "domains");
    lvl_dict.insert("P".into(), "phylums");
    lvl_dict.insert("O".into(), "orders");
    lvl_dict.insert("C".into(), "classes");
    lvl_dict.insert("F".into(), "families");
    lvl_dict.insert("G".into(), "genuses");
    lvl_dict.insert("S".into(), "species");

    let abundance_lvl = lvl_dict
        .get(&args.level)
        .unwrap_or(&&args.level.as_str())
        .to_string();

    let branch = check_and_parse(&args.level)?;

    // 定义主级别的数组
    let main_lvls = ['R', 'K', 'D', 'P', 'C', 'O', 'F', 'G', 'S'];

    // 查找给定级别的索引
    let branch_lvl = main_lvls
        .iter()
        .position(|&x| x == args.level.chars().next().unwrap())
        .unwrap();

    let input_file = args.input.clone();
    check_report_file(&input_file)?;

    let mut stats = kraken::Stats::default();
    kraken::parse_kraken_report(
        &input_file,
        &args.level,
        args.threshold,
        branch,
        branch_lvl,
        &mut stats,
    )?;

    let kmer_distr = kraken::read_kmer_distribution(&args.kmer_distr, &stats);

    kraken::dfs_iterative(1, &mut stats, &args.level, kmer_distr);

    // For all genomes, map reads up to level
    for (_, value) in stats.map2lvl_taxids.iter() {
        stats.lvl_taxids.get_mut(&value.0).unwrap().3 += value.2;
    }

    // Sum all of the reads for the desired level -- use for fraction of reads
    let mut sum_all_reads = 0;
    for (_, value) in stats.lvl_taxids.iter() {
        let new_all_reads = value.1 + value.3;
        sum_all_reads += new_all_reads;
    }
    if sum_all_reads == 0 {
        panic!("Error: no reads found. Please check your Kraken report");
    }

    let mut file = BufWriter::new(File::create(&args.output)?);

    writeln!(file,
        "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads"
    )?;

    for (taxid, value) in stats.lvl_taxids.iter() {
        let new_all_reads = value.1 + value.3;
        let fraction_total_reads = new_all_reads as f64 / sum_all_reads as f64;
        writeln!(file,
            "{name}\t{taxid}\t{level}\t{kraken_assigned_reads}\t{added_reads}\t{tnew_est_reads}\t{tfraction_total_reads:.5}",
            name=value.0,
            taxid=taxid,
            level=args.level,
            kraken_assigned_reads=value.1,
            added_reads=value.3,
            tnew_est_reads=new_all_reads,
            tfraction_total_reads=fraction_total_reads
        )?;
    }

    println!("BRACKEN SUMMARY (Kraken report: {:?})", args.input);
    println!("    >>> Threshold: {} ", args.threshold);
    println!(
        "    >>> Number of {:?} in sample: {:?} ",
        abundance_lvl, stats.n_lvl_total
    );
    println!(
        "\t  >> Number of {:} with reads > threshold: {:} ",
        abundance_lvl, stats.n_lvl_est
    );
    println!(
        "\t  >> Number of {} with reads < threshold: {} ",
        abundance_lvl, stats.n_lvl_del
    );
    println!("    >>> Total reads in sample: {}", stats.total_reads);
    println!(
        "\t  >> Total reads kept at {} level (reads > threshold): {}",
        abundance_lvl, stats.kept_reads
    );
    println!(
        "\t  >> Total reads discarded ({} reads < threshold): {}",
        abundance_lvl, stats.ignored_reads
    );
    println!("\t  >> Reads distributed: {}", stats.distributed_reads);
    println!(
        "\t  >> Reads not distributed (eg. no {} above threshold): {}",
        abundance_lvl, stats.nondistributed_reads
    );
    println!("\t  >> Unclassified reads: {:}", stats.u_reads);
    println!("BRACKEN OUTPUT PRODUCED: {:?}", &args.output.display());

    let now: DateTime<Local> = Local::now();
    let time = now.format("%m-%d-%Y %H:%M:%S").to_string();
    println!("PROGRAM END TIME: {}", time);

    Ok(())
}

#[allow(dead_code)]
fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("Application error: {}", e);
    }
}
