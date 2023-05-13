use crate::{data_structs::{ Bubble, KmerNode }, graph_utils};
use fxhash::FxHashSet;

/// Find bubble given a source node
/// 
/// # Arguments
/// * `ref_nodes` - Vector of KmerNodes representing the reference graph
/// * `src` - Source node to start the search from
/// 
/// # Returns
/// * `Option<Bubble>` - `Bubble` struct containing the nodes inside the bubble if a bubble is found
pub fn find_bubbles_given_src(ref_nodes: &Vec<KmerNode>, src: &KmerNode) -> Option<Bubble> {
    let mut seen = FxHashSet::default();
    let mut visited = FxHashSet::default();
    let mut nodes_inside_bubble: Vec<u32> = Vec::new();
    let mut stack: Vec<u32> = Vec::new();
    stack.push(src.id);
    seen.insert(src.id);
    while !stack.is_empty() {
        let node_id = stack.pop();
        if node_id.is_none() {
            break;
        }
        let node_id = node_id.unwrap();
        visited.insert(node_id);
        nodes_inside_bubble.push(node_id);
        // remove from seen
        seen.remove(&node_id);
        let node = &ref_nodes[node_id as usize];
        if node.child_nodes.len() == 0 {
            // tip
            break;
        }
        for child in node.child_nodes.iter() {
            let child_node = &ref_nodes[*child as usize];
            if child_node.id == src.id {
                // loop
                break;
            }
            if !seen.contains(&child_node.id) {
                seen.insert(child_node.id);
            }
            let parents = &child_node.parent_nodes;
            // push u to stack if all parents are visited
            if parents.iter().all(|p| visited.contains(p)) {
                stack.push(child_node.id);
            }
        }
        if stack.len() == 1 && seen.len() == 1 {
            let t = stack.pop().unwrap();
            nodes_inside_bubble.push(t);
            // if no edge between t and s
            if !ref_nodes[t as usize].child_nodes.contains(&src.id) && !ref_nodes[src.id as usize].child_nodes.contains(&t) {
                let bubble = Bubble {
                    id: 0,
                    start: src.id,
                    end: t,
                    kmers: nodes_inside_bubble,
                    longest_path_length: None,
                    shortest_path_length: None,
                    parent_bubble: None
                };
                // println!("Found bubble with src: {} and sink: {}", src.id, t);
                return Some(bubble);
            } else {
                break;
            }
        }
    }
    return None;
}

/// Find bubbles in the reference graph
/// 
/// This algorithm is due to Onodera et al. (2016) and is described in the paper
/// "Detecting superbubbles in assembly graphs".
/// 
/// The algorithm runs in O(mn) time.
/// 
/// # Arguments
/// * `ref_nodes` - Vector of KmerNodes representing the reference graph
/// 
/// # Returns
/// * `Vec<Bubble>` - Vector of `Bubble` structs containing the nodes inside the bubble
pub fn find_bubbles(ref_nodes: &Vec<KmerNode>) -> Vec<Bubble> {
    let mut bubbles: Vec<Bubble> = Vec::new();
    let mut bubble_id = 0;
    for node in ref_nodes.iter() {
        if node.child_nodes.len() > 1 {
            let bubble = find_bubbles_given_src(ref_nodes, node);
            if !bubble.is_none() {
                let mut bubble = bubble.unwrap();
                bubble.id = bubble_id;
                bubble_id += 1;
                bubble.longest_path_length = Some(graph_utils::longest_path_length(&ref_nodes, &bubble.kmers, &bubble.start, &bubble.end));
                bubble.shortest_path_length = Some(graph_utils::shortest_path_length(&ref_nodes, &bubble.kmers, &bubble.start, &bubble.end));
                // clear kmers to save memory
                bubble.kmers = Vec::new();
                bubbles.push(bubble);
            }
        }
    }
    return bubbles;
}