use debruijn::kmer::Kmer16;
use smallvec::SmallVec;
use serde::{Serialize, Deserialize};


//Use the SmallVec impelementation to save lots of memory 
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerNode{
    pub kmer: Kmer16,
    pub order: u32,
    pub color: u64,
    pub child_nodes: SmallVec<[u32;1]>,
    pub child_edge_distance: SmallVec<[(u8,(u64,u8));1]>,
//    pub child_nodes: Vec<u32>,
    pub id: u32,
    pub canonical: bool,
    pub actual_ref_positions: SmallVec<[usize;0]>,
}
