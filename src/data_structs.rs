use debruijn::kmer::Kmer16;
#[derive(Debug)]
pub struct KmerNode{
    pub kmer: Kmer16,
    pub id: usize,
    pub color: u64,
    pub child_nodes: Vec<usize>
}
