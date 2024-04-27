# Bracken-rust

## 0.Installation Instructions
To install bracken-rust, follow these steps:

1. Download the appropriate version for your system:

Navigate to the Releases page of the bracken-rust GitHub repository.
Select the release suitable for your operating system. For example, if you are using CentOS 7, download bracken-rust-${VERSION}-centos7.tar.gz, where ${VERSION} is the version number of the release you wish to install.

2. Extract the downloaded archive:

Open a terminal.
Use the tar command to extract the files from the archive

```bash
tar -xvf bracken-rust-${VERSION}-centos7.tar.gz
```


## 1. Get Started

```bash
$ ./bracken -h
Usage: bracken <COMMAND>

Commands:
  kmer2read-distr  bracken kmer2read_distr
  kmer-distrib     Evaluates genome read distribution and estimates reads per species for specific taxonomy IDs.
  est-abundance    Estimates species or genus level abundance from Kraken outputs using Bayesian methods.
  help             Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```
