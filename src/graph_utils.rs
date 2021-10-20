use crate::data_structs::KmerNode;
use std::cmp::Ordering;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use smallvec::SmallVec;

pub fn concat_graph<'a>(head: &KmerNode, ref_nodes: &Vec<KmerNode>) -> (Vec<(u32,u32,usize)>, FxHashMap<u32,Vec<u32>>) {

    let mut to_return = vec![];
    let mut bubble_nodes = FxHashSet::default();
    let mut in_edges = FxHashMap::default();
    let mut vertex_to_kmers_map = FxHashMap::default();

    bubble_nodes.insert(head.id);
    for node in ref_nodes.iter() {
        if node.child_nodes.len() > 1 {
            bubble_nodes.insert(node.id);
        }
        for child_id in node.child_nodes.iter() {
            let edges = in_edges.entry(child_id).or_insert(vec![]);
            edges.push(node.id);
        }
    }

    for (key, value) in in_edges {
        if value.len() > 1 {
            bubble_nodes.insert(*key);
        }
    }

    for start_node in bubble_nodes.iter(){
        let mut visited = FxHashSet::default();
        let mut nodes_to_visit = SmallVec::<[u32; 20]>::new();
        let mut len_edge = 0;
        for child_id in ref_nodes[*start_node as usize].child_nodes.iter() {
            nodes_to_visit.push(*child_id)
        }
        let mut intermediate_vertices = vec![];
        while nodes_to_visit.len() != 0 {
            len_edge += 1;
            let node = nodes_to_visit.pop().unwrap();
            if !bubble_nodes.contains(&node){
                intermediate_vertices.push(node);
                visited.insert(node);
                for child_id in ref_nodes[node as usize].child_nodes.iter() {
                    if visited.contains(&child_id) {
                        continue;
                    } else {
                        nodes_to_visit.push(*child_id);
                    }
                }
            }
            else{
                if intermediate_vertices.len() == 0{
                    to_return.push((ref_nodes[*start_node as usize].order, ref_nodes[node as usize].order, len_edge));
                }
                else{
                    to_return.push((ref_nodes[*start_node as usize].order, ref_nodes[intermediate_vertices[0] as usize].order, len_edge));
                    to_return.push((ref_nodes[intermediate_vertices[0] as usize].order, ref_nodes[node as usize].order, len_edge));
                    vertex_to_kmers_map.insert(intermediate_vertices[0],intermediate_vertices.clone());
                    intermediate_vertices.clear();
                }
                len_edge = 0;
            }
        }
    }

    return (to_return, vertex_to_kmers_map);
}

#[derive(Copy, Clone, Eq, PartialEq)]
struct State {
    inval: [u32;2],
    node_id: u32,
}

// The priority queue depends on `Ord`.
// Explicitly implement the trait so the queue becomes a min-heap
// instead of a max-heap.
impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        // Notice that the we flip the ordering on costs.
        // In case of a tie we compare positions - this step is necessary
        // to make implementations of `PartialEq` and `Ord` consistent.
        other.inval[1].cmp(&self.inval[1])
            .then_with(|| self.node_id.cmp(&other.node_id))
    }
}

// `PartialOrd` needs to be implemented as well.
impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
