use crate::taxonomy::NCBITaxonomy;
use dashmap::DashMap;
use memmap2::MmapOptions;
use rayon::prelude::*;
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Result, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

/// 读取 seqid2taxid.map 文件。为了裁剪 ncbi 的 taxonomy 树
pub fn get_seqid2taxid<P: AsRef<Path>>(filename: P) -> Result<HashMap<String, u32>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut id_map = HashMap::new();
    let mut s_count = 0;
    print!("\t>>STEP 1: READING SEQID2TAXID MAP\n");
    print!("\t\t0 sequences read");
    for line in reader.lines() {
        let line = line?;
        s_count += 1;
        if s_count % 1000 == 0 {
            print!("\r\t\t{:?} sequences read", s_count);
        }

        let parts: Vec<&str> = line.trim().split_whitespace().collect();
        if parts.len() < 2 {
            continue;
        }
        let seq_id = parts[0].to_string();
        if let Ok(taxid) = parts[1].parse::<u32>() {
            id_map.insert(seq_id, taxid);
        }
    }
    print!("\r\t\t{:?} total sequences read\n", s_count);

    Ok(id_map)
}

fn convert_line(
    line: &str,
    seqid2taxid: &HashMap<String, u32>,
    n_kmers: usize,
    taxo: &NCBITaxonomy,
) -> Option<(String, String)> {
    // 处理行的逻辑，替换为适合你需求的处理过程
    // let mut taxids_mapped: HashMap<u32, usize> = HashMap::new();
    let taxid_map: DashMap<u32, usize> = DashMap::new();

    let fields: Vec<_> = line.split("\t").collect();
    if fields.len() < 5 {
        return None;
    }
    let seqid = fields[1].trim();
    let taxid = seqid2taxid.get(seqid).unwrap_or(&0);
    let mut output = String::new();
    output += &format!("{}\t{}\t\t", seqid, taxid);

    let curr_ks: Vec<u32> = fields[4]
        .trim()
        .split(" ")
        .flat_map(|item| {
            let pair: Vec<_> = item.trim().split(":").collect();
            let taxid = pair[0].parse::<u32>().unwrap_or(0);
            let count = pair[1].parse::<usize>().unwrap_or(0);
            if count >= n_kmers {
                *taxid_map.entry(taxid).or_insert(0) += count - n_kmers + 1;
                std::iter::repeat(taxid)
                    .take(n_kmers - 1)
                    .collect::<Vec<_>>()
                    .into_iter()
            } else {
                std::iter::repeat(taxid)
                    .take(count)
                    .collect::<Vec<_>>()
                    .into_iter()
            }
        })
        .collect();

    let mut curr_kmers = VecDeque::new();
    let mut taxid2kmers = HashMap::new();
    let mut pre_mer: Option<u32> = None;
    let mut pre_taxid: u32 = 0;
    for kmer in curr_ks.iter() {
        curr_kmers.push_back(kmer);
        *taxid2kmers.entry(*kmer).or_insert_with(|| 0) += 1;
        if curr_kmers.len() == n_kmers {
            if pre_mer == Some(*kmer) {
                *taxid_map.entry(pre_taxid).or_insert(0) += 1;
            } else {
                let mapped_taxid = get_classification(&taxid2kmers, taxo);
                pre_taxid = mapped_taxid;
                *taxid_map.entry(mapped_taxid).or_insert(0) += 1;
            }
            if let Some(cur) = curr_kmers.pop_front() {
                pre_mer = Some(*cur);
                let count = taxid2kmers.entry(*cur).or_default();
                *count -= 1;
                if *count == 0 {
                    taxid2kmers.remove(&cur);
                }
            }
        }
    }
    // let k_mers = Kmers::new(curr_ks, n_kmers);
    // k_mers.par_bridge().for_each(|taxid2kmers| {
    //     let mapped_taxid = get_classification(&taxid2kmers, taxo);
    //     *taxid_map.entry(mapped_taxid).or_insert(0) += 1;
    // });

    output += &taxid_map
        .iter()
        .map(|item| format!("{}:{}", item.key(), item.value()))
        .collect::<Vec<_>>()
        .join(" ");

    output += "\n";
    Some((seqid.to_string(), output))
}

fn get_classification(taxid2kmers: &HashMap<u32, usize>, taxo: &NCBITaxonomy) -> u32 {
    if taxid2kmers.len() == 1 {
        if let Some((&taxid, _)) = taxid2kmers.iter().next() {
            return taxid;
        }
    }

    let mut max_score = 0;
    let mut max_taxid = 0;

    for (&taxid, &count) in taxid2kmers.iter() {
        if taxid == 0 {
            continue;
        }

        let score = if let Some(node) = taxo.get_node(&taxid) {
            node.path_to_root
                .iter()
                .filter_map(|&ancestor| taxid2kmers.get(&ancestor))
                .sum::<usize>()
                + count
        } else {
            count
        };

        if score > max_score {
            max_score = score;
            max_taxid = taxid;
        } else if score == max_score && max_taxid != 0 {
            max_taxid = taxo.lca(max_taxid, taxid);
        }
    }

    max_taxid
}

const BATCH_SIZE: usize = 100;

pub fn evaluate_kfile<P: AsRef<Path>>(
    k_file: P,
    o_file: P,
    seqid2taxid: HashMap<String, u32>,
    read_len: usize,
    kmer_len: usize,
    taxo: &NCBITaxonomy,
) -> Result<()> {
    print!("\t>>STEP 3: CONVERTING KMER MAPPINGS INTO READ CLASSIFICATIONS:\n");
    print!(
        "\t\t{}mers, with a database built using {}mers\n",
        read_len, kmer_len,
    );

    let file = File::open(k_file)?;
    let mmap = unsafe { MmapOptions::new().map(&file)? };
    let data = unsafe { std::str::from_utf8_unchecked(&mmap) };

    let outfile = File::create(o_file)?;
    let writer = Arc::new(Mutex::new(BufWriter::new(outfile)));

    /*Initialize variables for getting read mappings instead of kmer mappings */
    let n_kmers = read_len - kmer_len + 1;
    let counter = AtomicUsize::new(1);

    print!("\t\t0 sequences converted...");

    let buffer = Arc::new(Mutex::new(Vec::new()));
    data.par_lines().for_each(|line| {
        if let Some((seqid, output)) = convert_line(line, &seqid2taxid, n_kmers, taxo) {
            let count = counter.fetch_add(1, Ordering::SeqCst);
            print!("\r\t\t{} sequences converted (finished: {})", count, seqid);
            let mut buffer = buffer.lock().unwrap();
            buffer.extend_from_slice(output.as_bytes());

            if buffer.len() >= BATCH_SIZE {
                // Acquire the lock and write the buffer contents
                let mut write = writer.lock().unwrap();
                write.write_all(&buffer).expect("write data error");
                write.flush().expect("flush writer error");

                // Clear the buffer
                buffer.clear();
            }
        }
    });

    let buffer = buffer.lock().unwrap();
    if !buffer.is_empty() {
        let mut write = writer.lock().unwrap();
        write.write_all(&buffer).expect("write data error");
        write.flush().expect("flush writer error");
    }
    Ok(())
}
