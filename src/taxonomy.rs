use memmap2::Mmap;
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::HashMap;
use std::error::Error;
use std::fmt::{self, Debug};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct TaxonomyError {
    pub message: String,
}

impl TaxonomyError {
    pub fn new(message: &str) -> Self {
        TaxonomyError {
            message: message.to_string(),
        }
    }
}

impl fmt::Display for TaxonomyError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "TaxonomyError: {}", self.message)
    }
}

impl Error for TaxonomyError {}

impl From<std::io::Error> for TaxonomyError {
    fn from(_: std::io::Error) -> Self {
        TaxonomyError::new("io Error: Failed to handle nodes file")
    }
}

impl From<std::num::ParseIntError> for TaxonomyError {
    fn from(_: std::num::ParseIntError) -> Self {
        TaxonomyError::new("Failed to parse nodes file line")
    }
}

impl From<serde_json::Error> for TaxonomyError {
    fn from(_: serde_json::Error) -> Self {
        TaxonomyError::new("Failed to handle serde_json error")
    }
}

#[derive(Serialize, Deserialize)]
pub struct BiMap<T>
where
    T: Debug + Copy + Eq + PartialEq + Hash,
{
    forward: HashMap<T, T>,
    backward: HashMap<T, T>,
}

impl<T> BiMap<T>
where
    T: Debug + Copy + Eq + PartialEq + Hash,
{
    pub fn new() -> Self {
        BiMap {
            forward: HashMap::new(),
            backward: HashMap::new(),
        }
    }

    pub fn insert(&mut self, key: T, value: T) {
        self.forward.insert(key, value);
        self.backward.insert(value, key);
    }

    pub fn get_by_key(&self, key: &T) -> Option<&T> {
        self.forward.get(key)
    }

    pub fn get_by_value(&self, value: &T) -> Option<&T> {
        self.backward.get(value)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonomyNode {
    pub taxid: u32,
    pub parent: u32,
    pub rank: String,
    pub depth: u32,
    pub path_to_root: Vec<u32>,
}

impl TaxonomyNode {
    fn new(taxid: u32, parent: u32, rank: String, depth: u32) -> Result<Self, TaxonomyError> {
        if taxid != 1 && taxid == parent {
            return Err(TaxonomyError::new(
                "taxid should not be equal to parent unless it's 1",
            ));
        }

        Ok(Self {
            taxid,
            parent,
            rank,
            depth,
            path_to_root: vec![],
        })
    }

    pub fn is_root(&self) -> bool {
        self.taxid == 1 || self.depth == 1
    }
}

impl Default for TaxonomyNode {
    fn default() -> Self {
        Self::new(1, 1, "N".into(), 0).unwrap()
    }
}

#[derive(Serialize, Deserialize)]
pub struct NCBITaxonomy {
    pub nodes: Vec<TaxonomyNode>,
    pub id_map: BiMap<u32>,
}

impl Default for NCBITaxonomy {
    fn default() -> Self {
        Self {
            nodes: Vec::new(),
            id_map: BiMap::new(),
        }
    }
}

impl NCBITaxonomy {
    pub fn save_to_file<P: AsRef<Path>>(&self, path: P) -> Result<(), TaxonomyError> {
        let mut file = File::create(path)?;
        let json = serde_json::to_string(self)?;
        writeln!(file, "{}", json)?;
        Ok(())
    }

    fn update_depth_path(&mut self) {
        // 首先，为每个节点计算深度，并存储在一个Vec中
        let depths_paths: Vec<(u32, Vec<u32>)> = self
            .nodes
            .iter()
            .map(|node| {
                let mut depth = 1;
                let mut current_taxid = node.taxid;
                let mut path = vec![];
                while current_taxid != 1 {
                    if let Some(&parent_index) = self.id_map.get_by_key(&current_taxid) {
                        current_taxid = self.nodes[parent_index as usize].parent;
                        path.push(current_taxid);
                        depth += 1;
                    } else {
                        break;
                    }
                }
                (depth, path)
            })
            .collect();

        // 然后，使用收集到的深度值更新每个节点
        for (node, depth_path) in self.nodes.iter_mut().zip(depths_paths.iter()) {
            node.depth = depth_path.0;
            node.path_to_root = depth_path.1.clone();
            node.path_to_root.reverse();
        }
    }

    pub fn get_parent(&self, taxid: &u32) -> Option<&TaxonomyNode> {
        self.get_node(taxid)
            .and_then(|node| self.get_node(&node.parent))
    }

    pub fn get_node_parent(&self, node: &TaxonomyNode) -> Option<&TaxonomyNode> {
        self.get_node(&node.parent)
    }

    pub fn get_node(&self, taxid: &u32) -> Option<&TaxonomyNode> {
        self.id_map
            .get_by_key(taxid)
            .and_then(|&nodeid| self.nodes.get(nodeid as usize))
    }

    pub fn lca(&self, a: u32, b: u32) -> u32 {
        if a == 0 || b == 0 || a == b {
            return if a != 0 { a } else { b };
        }

        let na = self.get_node(&a).unwrap();
        let nb = self.get_node(&b).unwrap();

        let path_a = &na.path_to_root;
        let path_b = &nb.path_to_root;

        let mut i = 0;
        while i < path_a.len() && i < path_b.len() && path_a[i] == path_b[i] {
            i += 1;
        }

        if i == 0 {
            return 0;
        }

        // 返回最后一个共同的祖先
        *path_a.get(i - 1).unwrap_or(&0)
    }

    pub fn load_ncbi_dmp<P: AsRef<Path>>(node_file: P) -> Result<NCBITaxonomy, TaxonomyError> {
        let nodes_file = std::fs::File::open(node_file)?;

        let mut ncbi_taxo = NCBITaxonomy::default();

        for (ix, line) in BufReader::new(nodes_file).lines().enumerate() {
            let line = line?;
            let fields: Vec<_> = line.split("\t|\t").collect();
            if fields.len() < 10 {
                // should be at least 14
                return Err(TaxonomyError::new(
                    "Not enough fields in nodes.dmp; bad line?",
                ));
            }
            let taxid = fields[0].trim().parse::<u32>()?;
            let parent = fields[1].trim().parse::<u32>()?;
            let rank = fields[2].trim().to_string();

            let depth = if taxid == 1 { 1 } else { 0 };
            let node = TaxonomyNode::new(taxid, parent, rank, depth)?;
            ncbi_taxo.nodes.push(node);
            ncbi_taxo.id_map.insert(taxid, ix as u32);
        }

        ncbi_taxo.update_depth_path();

        Ok(ncbi_taxo)
    }

    pub fn load_from_json<P: AsRef<Path>>(path: P) -> Result<Self, TaxonomyError> {
        let file = File::open(path).map_err(|e| TaxonomyError::new(&e.to_string()))?;
        // let reader = BufReader::new(file);

        let mmap = unsafe { Mmap::map(&file) }.map_err(|e| TaxonomyError::new(&e.to_string()))?;

        let taxo =
            serde_json::from_slice(&mmap[..]).map_err(|e| TaxonomyError::new(&e.to_string()))?;

        // 直接从 Buffered Reader 反序列化 JSON 数据
        // let taxo =
        //     serde_json::from_reader(reader).map_err(|e| TaxonomyError::new(&e.to_string()))?;
        Ok(taxo)
    }

    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, TaxonomyError> {
        match path.as_ref().extension().and_then(|s| s.to_str()) {
            Some("json") => Self::load_from_json(path),
            Some("dmp") => Self::load_ncbi_dmp(path),
            _ => Err(TaxonomyError::new("Unsupported file format")),
        }
    }
}

pub fn load_taxonomy(taxonomy_dir: PathBuf) -> Result<NCBITaxonomy, TaxonomyError> {
    let json_file = taxonomy_dir.join("nodes.json");
    let dmp_file = taxonomy_dir.join("nodes.dmp");

    print!("\t>>STEP 2: READING NODES.DMP FILE\n");
    let taxo = if json_file.exists() {
        NCBITaxonomy::load(json_file)?
    } else if dmp_file.exists() {
        let taxo = NCBITaxonomy::load(dmp_file)?;
        taxo.save_to_file(&json_file)?;
        taxo
    } else {
        return Err(TaxonomyError::new("No suitable taxonomy nodes file found."));
    };

    print!("\t\t{:?} total nodes read\n", taxo.nodes.len());
    Ok(taxo)
}
