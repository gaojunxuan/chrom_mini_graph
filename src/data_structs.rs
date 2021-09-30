use debruijn::kmer::Kmer16;
use smallvec::SmallVec;


#[derive(Debug)]
pub struct KmerNode{
    pub kmer: Kmer16,
    pub order: u32,
    pub color: u64,
    pub child_nodes: SmallVec<[u32;1]>,
    pub id: u32, // may not need the id?
}
