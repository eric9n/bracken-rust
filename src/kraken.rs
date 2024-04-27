use core::str;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

pub const MAIN_LVLS: &[char; 9] = &['R', 'K', 'D', 'P', 'C', 'O', 'F', 'G', 'S'];

#[derive(Clone, Debug)]
pub struct Node {
    pub name: String,
    pub taxid: u32,
    pub level_id: String,
    pub level_num: usize,
    pub all_reads: usize,
    pub lvl_reads: usize,
    pub children: Vec<usize>,
    pub parent: Option<usize>,
}

impl Node {
    pub fn from_str(curr_str: &str) -> Result<Self, &'static str> {
        let split_str: Vec<&str> = curr_str.trim().split('\t').collect();
        if split_str.len() < 5 {
            return Err("Input string has too few elements.");
        }

        let mut level_id = split_str[split_str.len() - 3].to_string();
        if level_id.trim().is_empty() {
            level_id = "-".to_string();
        }

        let lvl_reads = split_str[2]
            .parse::<usize>()
            .map_err(|_| "Invalid number for lvl_reads")?;

        let mut name = split_str[split_str.len() - 1].to_string();
        let spaces = name.chars().take_while(|&c| c == ' ').count();
        name = name.trim_start().to_string();

        // 尝试将读数转换为整数
        let all_reads = split_str[1]
            .parse::<usize>()
            .map_err(|_| "Invalid number for all_reads")?;

        let taxid = split_str[split_str.len() - 2]
            .parse::<u32>()
            .map_err(|_| "Invalid number for taxid")?;

        // 计算层级编号
        let level_num = spaces / 2;

        Ok(Self {
            name,
            taxid,
            level_id,
            level_num,
            all_reads,
            lvl_reads,
            children: Vec::new(),
            parent: None,
        })
    }

    // 添加子节点
    pub fn add_child(&mut self, indx: usize) {
        self.children.push(indx);
    }
}

pub fn parse_kraken_report(
    input_file: &PathBuf,
    level: &str,
    threshold: usize,
    branch: usize,
    branch_lvl: usize,
    stats: &mut Stats,
) -> Result<(), Box<dyn std::error::Error>> {
    let r_file = File::open(&input_file)?;
    let reader = BufReader::new(r_file);

    let mut prev_node_index: usize = 0; // 使用索引而非引用
                                        // let mut leaf_nodes = Vec::new();
                                        // let mut nodes: HashMap<usize, Node> = HashMap::new();
    for (indx, line) in reader.lines().enumerate().into_iter() {
        let line = line?;
        if line.is_empty() || line.starts_with("#") || line.starts_with("%") {
            continue;
        }
        if let Ok(mut node) = Node::from_str(&line) {
            stats.total_reads += node.lvl_reads;
            if node.level_id == "U" || node.name == "unclassified" {
                stats.u_reads = node.lvl_reads;
                continue;
            }

            stats.nodes.insert(indx, node.clone());
            if node.taxid == 1 {
                prev_node_index = indx; // 保存索引
                continue;
            }

            stats.nodes.get_mut(&indx).unwrap().parent = Some(prev_node_index);

            let mut prev_node = stats.nodes.get(&prev_node_index).unwrap();
            if node.level_num != prev_node.level_num + 1 {
                stats.leaf_nodes.push(prev_node.clone());
            }
            while node.level_num != prev_node.level_num + 1 {
                if let Some(pre) = stats.nodes.get(&prev_node_index) {
                    prev_node = pre;
                    prev_node_index = pre.parent.expect("node has no parent");
                } else {
                    break;
                }
            }

            let (level_id, test_branch) = correct_level_id(&node, &prev_node);
            abundance_est(
                &node,
                level,
                threshold,
                test_branch,
                branch,
                branch_lvl,
                stats,
            );
            node.level_id = level_id;
            stats
                .nodes
                .get_mut(&prev_node_index)
                .unwrap()
                .add_child(indx);
            stats.nodes.get_mut(&indx).unwrap().parent = Some(prev_node_index);
            prev_node_index = indx;
        }
    }

    stats
        .leaf_nodes
        .push(stats.nodes.get(&prev_node_index).unwrap().clone());
    Ok(())
}

fn correct_level_id(node: &Node, prev_node: &Node) -> (String, usize) {
    let mut level_id = node.level_id.clone();
    let mut test_branch = 0;
    if level_id == "-" || level_id.len() > 1 {
        if MAIN_LVLS.contains(&prev_node.level_id.chars().next().unwrap()) {
            level_id = format!("{}1", prev_node.level_id);
            test_branch = 1;
        } else {
            if let Some(last_char) = prev_node.level_id.chars().last() {
                if let Some(num) = last_char.to_digit(10) {
                    let new_num = num + 1;
                    test_branch = new_num as usize;
                    level_id = format!(
                        "{}{}",
                        &prev_node.level_id[..prev_node.level_id.len() - 1],
                        new_num
                    );
                }
            }
        }
    }

    (level_id, test_branch)
}

#[derive(Debug)]
/// name,all_reads,level_reads,add_reads
pub struct LvlValue(pub String, pub usize, pub usize, pub usize);

impl LvlValue {
    pub fn from_node(node: &Node) -> Self {
        Self(node.name.clone(), node.all_reads, node.lvl_reads, 0)
    }
}

#[derive(Debug)]
/// taxid,level_reads,add_reads
pub struct Map2LvlValue(pub u32, pub usize, pub usize);

impl Map2LvlValue {
    pub fn from_node(node: &Node) -> Self {
        Self(node.taxid, node.lvl_reads, 0)
    }
}

fn main_lvl_index(level_id: &str) -> usize {
    MAIN_LVLS
        .iter()
        .position(|&x| x == level_id.chars().next().unwrap())
        .unwrap()
}

#[derive(Debug)]
pub struct Stats {
    pub u_reads: usize,
    pub n_lvl_total: usize,
    pub n_lvl_del: usize,
    pub ignored_reads: usize,
    pub total_reads: usize,
    pub n_lvl_est: usize,
    pub kept_reads: usize,
    pub last_taxid: isize,
    pub lvl_taxids: HashMap<u32, LvlValue>,
    pub map2lvl_taxids: HashMap<u32, Map2LvlValue>,
    pub leaf_nodes: Vec<Node>,
    pub nodes: HashMap<usize, Node>,
    pub nondistributed_reads: usize,
    pub distributed_reads: usize,
}

impl Default for Stats {
    fn default() -> Self {
        Self {
            u_reads: 0,
            n_lvl_total: 0,
            n_lvl_del: 0,
            ignored_reads: 0,
            n_lvl_est: 0,
            total_reads: 0,
            kept_reads: 0,
            last_taxid: -1,
            lvl_taxids: HashMap::new(),
            map2lvl_taxids: HashMap::new(),
            leaf_nodes: Vec::new(),
            nodes: HashMap::new(),
            nondistributed_reads: 0,
            distributed_reads: 0,
        }
    }
}

pub fn abundance_est(
    node: &Node,
    level: &str,
    threshold: usize,
    test_branch: usize,
    branch: usize,
    branch_lvl: usize,
    stats: &mut Stats,
) {
    let level_id = node.level_id.clone();

    let mut should_insert_map2lvl = false;

    if level_id == level {
        stats.n_lvl_total += 1;
        if node.all_reads < threshold {
            stats.n_lvl_del += 1;
            stats.ignored_reads += node.all_reads;
            stats.last_taxid = -1;
        } else {
            stats.n_lvl_est += 1;
            stats.kept_reads += node.all_reads;
            stats
                .lvl_taxids
                .insert(node.taxid.clone(), LvlValue::from_node(node));
            stats.last_taxid = node.taxid as isize; // 假设 taxid 是 usize，这里需要确保类型一致
            should_insert_map2lvl = true;
        }
    } else if branch > 0 && test_branch > branch {
        if stats.last_taxid != -1 {
            should_insert_map2lvl = true;
        }
    } else if main_lvl_index(&node.level_id) >= branch_lvl {
        if stats.last_taxid != -1 {
            should_insert_map2lvl = true;
        }
    }

    // 根据标志进行单次插入
    if should_insert_map2lvl {
        stats
            .map2lvl_taxids
            .insert(node.taxid, Map2LvlValue::from_node(node));
    }
}

fn process_kmer_distribution(
    curr_str: &str,
    stats: &Stats,
) -> Option<(u32, HashMap<u32, Vec<f32>>)> {
    let split_str: Vec<&str> = curr_str.trim().split('\t').collect();
    let mut temp_dict: HashMap<u32, Vec<f32>> = HashMap::new();

    if let Ok(mapped_taxid) = split_str[0].parse::<u32>() {
        if split_str.len() > 1 {
            for genome_str in split_str[1].split_whitespace() {
                let parts: Vec<&str> = genome_str.split(':').collect();
                if parts.len() == 3 {
                    if let (Ok(g_taxid), Ok(mkmers), Ok(tkmers)) = (
                        parts[0].parse::<u32>(),
                        parts[1].parse::<f32>(),
                        parts[2].parse::<f32>(),
                    ) {
                        let fraction = mkmers / tkmers;
                        if stats.lvl_taxids.contains_key(&g_taxid)
                            || stats.map2lvl_taxids.contains_key(&g_taxid)
                        {
                            temp_dict
                                .entry(g_taxid)
                                .or_insert_with(Vec::new)
                                .push(fraction);
                        }
                    }
                }
            }
        }
        return Some((mapped_taxid, temp_dict));
    }
    None
}

pub fn read_kmer_distribution(
    filename: &PathBuf,
    stats: &Stats,
) -> HashMap<u32, HashMap<u32, Vec<f32>>> {
    let file = File::open(filename).expect("Unable to open file");
    let reader = BufReader::new(file);
    let mut kmer_distr: HashMap<u32, HashMap<u32, Vec<f32>>> = HashMap::new();

    for line in reader.lines().skip(1) {
        if let Ok(line) = line {
            if let Some((mapped_taxid, mapped_taxid_dict)) = process_kmer_distribution(&line, stats)
            {
                if !mapped_taxid_dict.is_empty() {
                    kmer_distr.insert(mapped_taxid, mapped_taxid_dict);
                }
            }
        }
    }
    kmer_distr
}

pub fn dfs_iterative(
    root_index: usize,
    stats: &mut Stats,
    level: &str,
    kmer_distr: HashMap<u32, HashMap<u32, Vec<f32>>>,
) {
    let mut stack = vec![root_index];

    while let Some(node_index) = stack.pop() {
        if let Some(node) = stats.nodes.get(&node_index) {
            if node.level_id == level {
                continue;
            }
            // 为保持深度优先顺序，需要逆序推入子节点
            for child_index in node.children.iter().rev() {
                stack.push(*child_index);
            }
            // No reads to distribute
            if node.lvl_reads == 0 {
                continue;
            }
            // No genomes produce this classification
            if !kmer_distr.contains_key(&node.taxid) {
                stats.nondistributed_reads += node.lvl_reads;
                continue;
            }

            stats.distributed_reads += node.lvl_reads;
            let curr_dict = kmer_distr.get(&node.taxid).unwrap();
            let mut all_genome_reads = 0;
            let mut probability_dict_prelim = HashMap::<u32, (f32, usize)>::new();
            for (genome, value) in curr_dict {
                // Get the fraction of kmers of the genome expected to map to this node
                let fraction: f32 = value[0];
                // Determine the number of reads classified by Kraken uniquely for the genome
                // and the fraction of the genome that is unique
                let num_classified_reads = stats.map2lvl_taxids.get(genome).unwrap().1;

                let lvl_fraction = if kmer_distr.contains_key(genome)
                    && kmer_distr.get(genome).unwrap().contains_key(genome)
                {
                    kmer_distr.get(genome).unwrap().get(genome).unwrap()[0]
                } else {
                    1.0
                };

                let est_genome_reads = (num_classified_reads as f64 / lvl_fraction as f64) as usize;
                all_genome_reads += est_genome_reads;
                probability_dict_prelim.insert(*genome, (fraction, est_genome_reads));
            }
            if all_genome_reads == 0 {
                continue;
            }
            // # Get final probabilities
            // # P_R_A = probability that a read is classified at the node given that it belongs to genome A
            // # P_A = probability that a randomly selected read belongs to genome A
            // # P_A_R = probability that a read belongs to genome A given that its classified at the node
            let mut total_probability = 0.0;
            let mut probability_dict_final = HashMap::new();
            for (genome, value) in probability_dict_prelim.iter() {
                let p_a = value.1 as f64 / all_genome_reads as f64;
                let p_a_r = value.0 as f64 * p_a;
                probability_dict_final.insert(genome, p_a_r);
                total_probability += p_a_r;
            }

            // Find the normalize probabilty and Distribute reads accordingly
            for (genome, value) in probability_dict_final.iter() {
                let add_fraction = value / total_probability;
                let add_reads = (add_fraction / node.lvl_reads as f64) as usize;
                stats.map2lvl_taxids.get_mut(&genome).unwrap().2 += add_reads;
            }
        } else {
            println!("Node with index {} not found.", node_index);
        }
    }
}
