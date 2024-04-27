use bracken::ctime::timeval_subtract;
use bracken::{kmer2read_distr, taxonomy};
use clap::Parser;
use std::path::PathBuf;
use std::time::SystemTime;

#[derive(Parser, Clone, Debug)]
#[clap(version, about = "bracken kmer2read_distr", long_about)]
pub struct Args {
    /// seqid2taxid file generated during
    /// Kraken database building process
    #[clap(long)]
    seqid2taxid: PathBuf,

    /// taxonomy folder containing the nodes.dmp file
    /// (typically downloaded with the Kraken taxonomy)
    #[clap(long = "taxonomy")]
    taxonomy_dir: PathBuf,

    /// kraken file of all classifications of all library
    /// sequences (typically database.kraken)
    #[clap(long)]
    kraken: PathBuf,

    /// name of an output file to print read distributions to
    /// (suggested name: databaseXmers.kraken_cnts)
    #[clap(long)]
    output: PathBuf,

    /// kmer length used to build Kraken database
    /// (default = 31)
    #[clap(short = 'k', default_value_t = 31)]
    kmer_len: usize,

    /// read length (evaluate every l-length read)
    /// (default = 100)
    #[clap(short = 'l', default_value_t = 100)]
    read_len: usize,

    /// number of threads
    /// (default = 1)
    #[clap(short = 't', default_value_t = 1)]
    threads: usize,
}

pub fn run(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    let ta = SystemTime::now();
    println!("\t>>STEP 0: PARSING COMMAND LINE ARGUMENTS");
    let taxonomy_dir = args.taxonomy_dir;
    let json_file = taxonomy_dir.join("nodes.json");
    let dmp_file = taxonomy_dir.join("nodes.dmp");
    if json_file.exists() {
        println!("\t\tTaxonomy nodes file: {:}", json_file.display());
    } else {
        println!("\t\tTaxonomy nodes file: {:}", dmp_file.display());
    }
    println!("\t\tSeqid file:          {:}", args.seqid2taxid.display());
    println!("\t\tNum Threads:         {:?}", args.threads);
    println!("\t\tKmer Length:         {:?}", args.kmer_len);
    println!("\t\tRead Length:         {:?}", args.read_len);

    let seq_tax_map = kmer2read_distr::get_seqid2taxid(args.seqid2taxid)?;
    let taxo = taxonomy::load_taxonomy(taxonomy_dir)?;

    kmer2read_distr::evaluate_kfile(
        args.kraken,
        args.output,
        seq_tax_map,
        args.read_len,
        args.kmer_len,
        &taxo,
    )?;

    let tb = SystemTime::now();

    // 计算时间差
    match timeval_subtract(ta, tb) {
        Ok(duration) => {
            let total_seconds = duration.as_secs();
            let minutes = total_seconds / 60; // 得到分钟数
            let seconds = total_seconds % 60; // 得到剩余秒数
            let microseconds = duration.subsec_micros(); // 得到微秒数

            println!(
                "\tTime Elapsed: {} minutes, {} seconds, {:.5} microseconds",
                minutes, seconds, microseconds as f64
            );
            println!("\t=============================");
        }
        Err(e) => println!("{}", e),
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
