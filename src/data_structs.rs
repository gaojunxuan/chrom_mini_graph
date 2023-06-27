use std::collections::{HashSet, BTreeSet, BTreeMap, HashMap};
use std::vec::Vec;
use debruijn::kmer::Kmer16;
use smallvec::SmallVec;
use serde::{Serialize, Deserialize};
use block_aligner::cigar::*;

//First is ref, second is query
pub type Anchors = Vec<(u32, u32)>;
pub type Color = u128;

//Use the SmallVec implementation to save lots of memory 
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerNode{
    pub kmer: Kmer16,
    pub order: u32,
    pub order_val: u32,
    pub color: Color,
    pub child_nodes: SmallVec<[u32;1]>,
    pub child_edge_distance: SmallVec<[(u16,(Color,u8));1]>,
    pub parent_nodes: SmallVec<[u32;1]>,
    pub id: u32,
    pub canonical: bool,
    pub actual_ref_positions: SmallVec<[usize;0]>,
    pub repetitive: bool,
    pub primary_base: Option<u32>,
    pub closest_ref: u32,
    pub dist_to_closest_ref: u16,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Bubble {
    pub id: u32,
    pub kmers: Vec<u32>,
    pub start: u32,
    pub end: u32,
    pub longest_path_length: Option<u32>,
    pub shortest_path_length: Option<u32>,
    pub parent_bubble: Option<u32>
}

pub struct BamInfo{
    pub cigar: Vec<OpLen>,
    pub sequence: String,
    pub quals: Vec<u8>,
    pub qname: String,
    pub strand: bool,
    pub ref_name: String,
    pub map_pos: i64,
    pub mapq: u8,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseMatrix {
    pub shape: (usize, usize),
    pub data: HashMap<(usize, usize), u32>,
    pub num_nonzero: usize,
}

impl SparseMatrix {
    pub fn new(shape: (usize, usize)) -> Self {
        Self {
            shape,
            data: HashMap::new(),
            num_nonzero: 0,
        }
    }

    pub fn get(&self, i: usize, j: usize) -> u32 {
        if self.data.contains_key(&(i,j)) {
            return self.data.get(&(i, j)).unwrap().clone();
        }
        u32::MAX
    }

    pub fn set(&mut self, i: usize, j: usize, val: u32) {
        self.data.insert((i, j), val);
        self.num_nonzero += 1;
    }
}