use clap::{Parser, Subcommand};

mod est_abundance;
mod kmer2read_distr;
mod kmer_distrib;

#[derive(Subcommand, Debug)]
enum Commands {
    Kmer2readDistr(kmer2read_distr::Args),
    KmerDistrib(kmer_distrib::Args),
    EstAbundance(est_abundance::Args),
}

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(subcommand)]
    cmd: Commands,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    match args.cmd {
        Commands::EstAbundance(cmd_args) => {
            est_abundance::run(cmd_args)?;
        }
        Commands::Kmer2readDistr(cmd_args) => {
            kmer2read_distr::run(cmd_args)?;
        }
        Commands::KmerDistrib(cmd_args) => {
            kmer_distrib::run(cmd_args)?;
        }
    }
    Ok(())
}
