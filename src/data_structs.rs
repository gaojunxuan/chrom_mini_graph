use debruijn::kmer::Kmer16;
use smallvec::SmallVec;
use serde::{Serialize, Deserialize};
use block_aligner::cigar::*;

//First is ref, second is query
pub type Anchors = Vec<(u32, u32)>;
pub type Color = u128;

//Use the SmallVec impelementation to save lots of memory 
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
    pub primary_base: Option<u32>
}
pub struct Bubble {
    pub id: u32,
    pub kmers: Vec<u32>,
    pub start: u32,
    pub end: u32,
    pub longest_path_length: u32,
    pub shortest_path_length: u32,
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
