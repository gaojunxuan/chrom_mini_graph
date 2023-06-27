use std::{collections::HashMap, iter::FromIterator};
use crate::{data_structs::{ Bubble }, graph_utils};
use cmg_shared::data_structs::{KmerNode};
use fxhash::FxHashSet;
use range_minimum_query::Rmq;

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
                // bubble.longest_path_length = Some(graph_utils::longest_path_length(&ref_nodes, &bubble.kmers, &bubble.start, &bubble.end));
                // bubble.shortest_path_length = Some(graph_utils::shortest_path_length(&ref_nodes, &bubble.kmers, &bubble.start, &bubble.end));
                // clear kmers to save memory
                bubble.kmers = Vec::new();
                bubbles.push(bubble);
            }
        }
    }
    return bubbles;
}

fn min_out_parent(ref_nodes: &Vec<KmerNode>, v: &KmerNode) -> usize {
    let mut minimum = usize::MAX;
    for parent in v.parent_nodes.iter() {
        minimum = minimum.min(ref_nodes[*parent as usize].order as usize);
    }
    return minimum;
}

fn max_out_child(ref_nodes: &Vec<KmerNode>, v: &KmerNode) -> usize {
    let mut maximum = usize::MIN;
    for child in v.child_nodes.iter() {
        maximum = maximum.max(ref_nodes[*child as usize].order as usize);
    }
    return maximum;
}

fn is_entrance(ref_nodes: &Vec<KmerNode>, pos: usize, top_sort: &Vec<u32>) -> bool {
    if pos + 1 == top_sort.len() {
        return false;
    }
    return pos == min_out_parent(ref_nodes, &ref_nodes[top_sort[pos + 1] as usize]);
}

fn is_exit(ref_nodes: &Vec<KmerNode>, pos: usize, top_sort: &Vec<u32>) -> bool {
    if pos == 0 {
        return false;
    }
    return pos == max_out_child(ref_nodes, &ref_nodes[top_sort[pos - 1] as usize]);
}

fn insert_exit(v: &KmerNode, candidates: &mut Vec<(u32, bool)>) {
    candidates.push((v.id, false));
}

fn insert_entrance(v: &KmerNode, candidates: &mut Vec<(u32, bool)>, backcand: &mut HashMap<u32, u32>) {
    backcand.entry(v.id).and_modify(|e| {*e = candidates.len() as u32}).or_insert(candidates.len() as u32);
    candidates.push((v.id, true));
}

fn delete_tail(candidates: &mut Vec<(u32, bool)>) {
    candidates.pop();
}

fn next_candidate(v: &KmerNode, candidates: &mut Vec<(u32, bool)>, backcand: &mut HashMap<u32, u32>) -> (u32, bool) {
    return candidates[backcand[&v.id] as usize + 1];
}

fn validate_superbubble(ref_nodes: &Vec<KmerNode>, top_sort: &Vec<u32>, start: &KmerNode, end: &KmerNode, out_child_rmq: &Rmq,
                        out_parent_rmq: &Rmq, out_children: &Vec<u32>, out_parents: &Vec<u32>, prev_entrance: &Vec<Option<u32>>) -> Option<u32> {
    let start_order = start.order;
    let end_order = end.order;
    let outchild_idx = out_child_rmq.range_minimum(start_order as usize..end_order as usize);
    let outparent_idx = out_parent_rmq.range_minimum((start_order + 1) as usize..(end_order + 1) as usize);
    if outchild_idx.is_none() || out_children[outchild_idx.unwrap()] != end_order || outparent_idx.is_none() {
        return None;
    }
    let out_parent = out_parents[outparent_idx.unwrap()];
    if out_parent == start_order {
        return Some(start.id);
    } else if is_entrance(ref_nodes, out_parent as usize, top_sort) {
        return Some(top_sort[out_parent as usize]);
    } else {
        return prev_entrance[out_parent as usize];
    }
}

fn report_superbubble(ref_nodes: &Vec<KmerNode>, top_sort: &Vec<u32>, start: Option<&KmerNode>, end: Option<&KmerNode>,
                      candidates: &mut Vec<(u32, bool)>, backcand: &mut HashMap<u32, u32>, 
                      prev_entrance: &Vec<Option<u32>>, alt_entrance: &mut Vec<Option<u32>>,
                      out_children: &Vec<u32>, out_parents: &Vec<u32>, out_child_rmq: &Rmq, out_parent_rmq: &Rmq,
                      bubbles: &mut Vec<Bubble>) {
    if start.is_none() || end.is_none() || start.unwrap().order >= end.unwrap().order {
        delete_tail(candidates);
        return;
    }
    let mut s: Option<u32> = prev_entrance[end.unwrap().order as usize];
    let mut valid = None;
    while s.is_some() && ref_nodes[s.unwrap() as usize].order >= start.unwrap().order {
        valid = validate_superbubble(ref_nodes, top_sort, start.unwrap(), end.unwrap(), out_child_rmq, out_parent_rmq,
                                     out_children, out_parents, prev_entrance);
        if (valid.is_some() && valid.unwrap() == s.unwrap()) 
             || (valid.is_some() && alt_entrance[s.unwrap() as usize].is_some() && valid.unwrap() == alt_entrance[s.unwrap() as usize].unwrap())
             || valid.is_none() {
            break;
        }
        alt_entrance[s.unwrap() as usize] = valid;
        s = valid;
    }
    delete_tail(candidates);
    if s.is_some() && (valid.is_some() && valid.unwrap() == s.unwrap()) {
        let s_node = &ref_nodes[s.unwrap() as usize];
        let t_node = &ref_nodes[end.unwrap().id as usize];
        // report (s, end) as a superbubble
        let mut bubble: Bubble = Bubble {
            id: bubbles.len() as u32,
            start: s.unwrap(),
            end: end.unwrap().id,
            kmers: top_sort[s_node.order as usize..=t_node.order as usize].to_vec(),
            longest_path_length: Some(graph_utils::longest_path_length(&ref_nodes, top_sort, &s.unwrap(), &end.unwrap().id)),
            shortest_path_length: Some(graph_utils::shortest_path_length(&ref_nodes, top_sort, &s.unwrap(), &end.unwrap().id)),
            parent_bubble: None
        };
        // if bubble.kmers.len() <= 10 {
        //     return;
        // }
        let kmer_count = bubble.kmers.len();
        bubble.kmers = Vec::new();
        // println!("Found bubble with src: {} and sink: {}; longest = {}, shortest = {}, s.ord = {}, t.ord = {}", s.unwrap(), end.unwrap().id, bubble.longest_path_length.unwrap(), bubble.shortest_path_length.unwrap(),
        //          s_node.order, t_node.order);
        if kmer_count > 10 {
            bubbles.push(bubble);
        }
        // then
        while candidates.len() > 0 && candidates[candidates.len() - 1].0 != s.unwrap() {
            if candidates[candidates.len() - 1].1 {
                delete_tail(candidates);
            } else {
                report_superbubble(ref_nodes, 
                                   top_sort, 
                                   Some(&ref_nodes[next_candidate(&ref_nodes[s.unwrap() as usize], candidates, backcand).0 as usize]), 
                                   Some(&ref_nodes[candidates[candidates.len() - 1].0 as usize]), 
                                   candidates, 
                                   backcand, 
                                   prev_entrance, 
                                   alt_entrance, 
                                   out_children, 
                                   out_parents, 
                                   out_child_rmq, 
                                   out_parent_rmq, 
                                   bubbles);
            }
        }

    }
}

pub fn find_superbubbles(ref_nodes: &Vec<KmerNode>, top_sort: &Vec<u32>, bubbles: &mut Vec<Bubble>) {
    let mut out_parents = Vec::new();
    let mut out_children = Vec::new();
    for k in 0..top_sort.len() {
        out_parents.push(min_out_parent(ref_nodes, &ref_nodes[top_sort[k] as usize]) as u32);
        out_children.push(max_out_child(ref_nodes, &ref_nodes[top_sort[k] as usize]) as u32);
    }
    // out_parent range minimum query
    let out_parent_rmq = Rmq::from_iter(&out_parents);
    // out_child range maximum query
    let out_child_rmq = Rmq::from_iter(out_children.iter().map(|x| -(*x as i32)).collect::<Vec<i32>>());

    let mut prev_ent: Option<u32> = None;
    
    let mut alt_entrance: Vec<Option<u32>> = vec![None; ref_nodes.len()];
    let mut prev_entrance: Vec<Option<u32>> = Vec::new();
    let mut candidates: Vec<(u32, bool)> = Vec::new();
    let mut backcand: HashMap<u32, u32> = HashMap::new();

    for j in 0..top_sort.len() {
        prev_entrance.push(prev_ent);
        let v = &ref_nodes[top_sort[j] as usize];
        alt_entrance[v.id as usize] = None;
        if is_exit(ref_nodes, j, top_sort) {
            insert_exit(v, &mut candidates);
        }
        if is_entrance(ref_nodes, j, top_sort) {
            insert_entrance(v, &mut candidates, &mut backcand);
            prev_ent = Some(v.id);
        }
    }
    while candidates.len() != 0 {
        if candidates[candidates.len() - 1].1 {
            delete_tail(&mut candidates);
        } else {
            report_superbubble(ref_nodes, top_sort, 
                               Some(&ref_nodes[candidates[0].0 as usize]), Some(&ref_nodes[candidates[candidates.len() - 1].0 as usize]),
                               &mut candidates, &mut backcand, &prev_entrance, 
                               &mut alt_entrance, &out_children, &out_parents, 
                               &out_child_rmq, &out_parent_rmq, bubbles);
        }
    }

}