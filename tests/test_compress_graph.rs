use chrom_mini_graph::graph_utils::compress_graph;
use cmg_shared::data_structs::{KmerNode,Bubble};
use chrom_mini_graph::data_structs::SimplifiedKmerNode;

#[cfg(test)]
mod tests {
    use debruijn::kmer::IntKmer;

    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_compress_graph() {
        let node1 = KmerNode{
            kmer: IntKmer { storage: 0 },
            order: 0,
            order_val: 0,
            color: 0,
            child_nodes: smallvec::smallvec![1],
            child_edge_distance: smallvec::smallvec![(0,(0,0))],
            parent_nodes: smallvec::smallvec![],
            id: 0,
            canonical: false,
            actual_ref_positions: smallvec::smallvec![],
            repetitive: false,
            primary_base: None,
            closest_ref: 0,
            dist_to_closest_ref: 0,
        };
        let node2 = KmerNode{
            kmer: IntKmer { storage: 0 },
            order: 1,
            order_val: 1,
            color: 0,
            child_nodes: smallvec::smallvec![2],
            child_edge_distance: smallvec::smallvec![(0,(0,0))],
            parent_nodes: smallvec::smallvec![0],
            id: 1,
            canonical: false,
            actual_ref_positions: smallvec::smallvec![],
            repetitive: false,
            primary_base: None,
            closest_ref: 0,
            dist_to_closest_ref: 0,
        };
        let node3 = KmerNode{
            kmer: IntKmer { storage: 0 },
            order: 2,
            order_val: 2,
            color: 0,
            child_nodes: smallvec::smallvec![3, 4],
            child_edge_distance: smallvec::smallvec![(0,(0,0)), (0,(0,0))],
            parent_nodes: smallvec::smallvec![2],
            id: 2,
            canonical: false,
            actual_ref_positions: smallvec::smallvec![],
            repetitive: false,
            primary_base: None,
            closest_ref: 0,
            dist_to_closest_ref: 0,
        };
        let node4 = KmerNode{
            kmer: IntKmer { storage: 0 },
            order: 3,
            order_val: 3,
            color: 0,
            child_nodes: smallvec::smallvec![],
            child_edge_distance: smallvec::smallvec![],
            parent_nodes: smallvec::smallvec![2],
            id: 3,
            canonical: false,
            actual_ref_positions: smallvec::smallvec![],
            repetitive: false,
            primary_base: None,
            closest_ref: 0,
            dist_to_closest_ref: 0,
        };
        let node5 = KmerNode{
            kmer: IntKmer { storage: 0 },
            order: 4,
            order_val: 4,
            color: 0,
            child_nodes: smallvec::smallvec![],
            child_edge_distance: smallvec::smallvec![],
            parent_nodes: smallvec::smallvec![2],
            id: 4,
            canonical: false,
            actual_ref_positions: smallvec::smallvec![],
            repetitive: false,
            primary_base: None,
            closest_ref: 0,
            dist_to_closest_ref: 0,
        };
        let mut test_graph = vec![node1, node2, node3, node4, node5];
        let mut cmg = compress_graph(&mut test_graph);
        print!("{:?}", cmg);
    }
}