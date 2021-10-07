use debruijn::kmer::Kmer16;
use smallvec::SmallVec;
use serde::{Serialize};


//Use the SmallVec impelementation to save lots of memory 
//if we don't care about serializing. (SmallVec can't be serialized). 
#[derive(Debug, Clone, Serialize)]
pub struct KmerNode{
    pub kmer: Kmer16,
    pub order: u32,
    pub color: u64,
//    pub child_nodes: SmallVec<[u32;1]>,
    pub child_nodes: Vec<u32>,
    pub id: u32, // may not need the id?
}
