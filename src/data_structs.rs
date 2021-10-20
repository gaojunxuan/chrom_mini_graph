use debruijn::kmer::Kmer16;
use smallvec::SmallVec;
use serde::{Serialize, Deserialize};


//Use the SmallVec impelementation to save lots of memory 
//if we don't care about serializing. (SmallVec can't be serialized). 
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerNode{
    pub kmer: Kmer16,
    pub order: u32,
    pub color: u64,
    pub child_nodes: SmallVec<[u32;1]>,
    pub child_edge_distance: SmallVec<[u8;1]>,
//    pub child_nodes: Vec<u32>,
    pub id: u32,
}
