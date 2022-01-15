use crate::data_structs::KmerNode;
use std::hash::Hash;
use std::collections::HashSet;
use debruijn::Kmer;
use debruijn::Mer;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use smallvec::SmallVec;
use std::cmp::Ordering;

fn has_unique_elements<T>(iter: T) -> bool
where
    T: IntoIterator,
    T::Item: Eq + Hash,
{
    let mut uniq = HashSet::new();
    iter.into_iter().all(move |x| uniq.insert(x))
}

pub fn concat_graph<'a>(
    head: &KmerNode,
    ref_nodes: &Vec<KmerNode>,
) -> (Vec<(u32, u32, usize)>, FxHashMap<u32, Vec<u32>>) {
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

    for start_node in bubble_nodes.iter() {
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
            if !bubble_nodes.contains(&node) {
                intermediate_vertices.push(node);
                visited.insert(node);
                for child_id in ref_nodes[node as usize].child_nodes.iter() {
                    if visited.contains(&child_id) {
                        continue;
                    } else {
                        nodes_to_visit.push(*child_id);
                    }
                }
            } else {
                if intermediate_vertices.len() == 0 {
                    to_return.push((
                        ref_nodes[*start_node as usize].order,
                        ref_nodes[node as usize].order,
                        len_edge,
                    ));
                } else {
                    to_return.push((
                        ref_nodes[*start_node as usize].order,
                        ref_nodes[intermediate_vertices[0] as usize].order,
                        len_edge,
                    ));
                    to_return.push((
                        ref_nodes[intermediate_vertices[0] as usize].order,
                        ref_nodes[node as usize].order,
                        len_edge,
                    ));
                    vertex_to_kmers_map
                        .insert(intermediate_vertices[0], intermediate_vertices.clone());
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
    inval: [u32; 2],
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
        other.inval[1]
            .cmp(&self.inval[1])
            .then_with(|| self.node_id.cmp(&other.node_id))
    }
}

// `PartialOrd` needs to be implemented as well.
impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub fn top_sort(ref_nodes: &mut Vec<KmerNode>) -> Vec<u32> {
    let mut nodes_to_visit = Vec::new();
    let mut stack_of_visited = Vec::new();
    stack_of_visited.push(0 as u32);
    nodes_to_visit.push(0 as u32);
    let mut visited = FxHashSet::default();
    let mut rev_sort_list = vec![];
    let mut order_to_id = vec![];
    let mut already_seen = FxHashSet::default();

    while nodes_to_visit.len() != 0 {
        let node = nodes_to_visit.pop().unwrap();
        visited.insert(node);
        let mut no_further = true;
        for child_id in ref_nodes[node as usize].child_nodes.iter() {
            //Circular cutoff
            if ref_nodes[(*child_id) as usize].order < ref_nodes[node as usize].order {
                continue;
            }
            if visited.contains(&child_id) {
                continue;
            } else {
                nodes_to_visit.push(*child_id);
                stack_of_visited.push(*child_id);
                no_further = false;
            }
        }
        if no_further {
            loop {
                //last node
                if nodes_to_visit.len() == 0 {
                    for _i in 0..stack_of_visited.len() {
                        let sorted_node = stack_of_visited.pop().unwrap();
                        if !already_seen.contains(&sorted_node) {
                            already_seen.insert(sorted_node);
                            rev_sort_list.push(sorted_node);
                        }
                    }
                    break;
                }

                if *stack_of_visited.last().unwrap() == *nodes_to_visit.last().unwrap() {
                    break;
                }
                let sorted_node = stack_of_visited.pop().unwrap();
                if !already_seen.contains(&sorted_node) {
                    rev_sort_list.push(sorted_node);
                    already_seen.insert(sorted_node);
                }
            }
        }
    }

    for (i, id) in rev_sort_list.iter().rev().enumerate() {
        ref_nodes[(*id) as usize].order = i as u32;
        order_to_id.push(*id);
    }

    return order_to_id;
}

pub fn add_align_to_graph(
    ref_nodes: &mut Vec<KmerNode>,
    aln_nodes: Vec<KmerNode>,
    mut anchors: Vec<(u32, u32)>,
    forward_strand: bool,
    samp_freq: usize,
    circular: bool
) {
    let mut strand_aln_nodes;
    //Reverse the array so that the orders correspond to the indices
    //(orders are swapped during chaining procedure for rev comps)
    if forward_strand == false {
        strand_aln_nodes = vec![];
        for node in aln_nodes.into_iter().rev() {
            strand_aln_nodes.push(node);
        }
    } else {
        strand_aln_nodes = aln_nodes;
    }
    //    let clone = ref_nodes.clone();
    let mut new_nodes = vec![];
    for node in ref_nodes.iter_mut() {
        node.color = node.color << 1;
        for edge in node.child_edge_distance.iter_mut() {
            edge.1 .0 = edge.1 .0 << 1;
        }
    }
    // TODO do this later, but we'll forget about circular rn because it's messy
//    if circular{
//        anchors.push(anchors[0]);
//    }
    let mut num_contains_dist_not = 0;
    for i in 0..anchors.len() - 1 {
        let ref_node_len = ref_nodes.len();

        let kmer1r;
        let kmer2r;
        let largest_anchor_id;

        let left;
        let right;
        if anchors[i].0 > anchors[i + 1].0 {
            largest_anchor_id = anchors[i].0;
            let (_left, _right) = ref_nodes.split_at_mut(largest_anchor_id as usize);
            left = _left;
            right = _right;
            kmer1r = &mut right[0];
            kmer2r = &mut left[anchors[i + 1].0 as usize];
        } else {
            largest_anchor_id = anchors[i + 1].0;
            let (_left, _right) = ref_nodes.split_at_mut(largest_anchor_id as usize);
            left = _left;
            right = _right;
            //            dbg!(anchors[i+1],anchors[i],i,i+1);
            kmer1r = &mut left[anchors[i].0 as usize];
            kmer2r = &mut right[0];
        }

        let ind_q_i;
        let ind_q_ip1;
        if forward_strand {
            ind_q_i = anchors[i].1 as usize;
            ind_q_ip1 = anchors[i + 1].1 as usize;
        } else {
            ind_q_i = strand_aln_nodes.len() - anchors[i].1 as usize - 1;
            ind_q_ip1 = strand_aln_nodes.len() - anchors[i + 1].1 as usize - 1;
        }

        let kmer1q = &strand_aln_nodes[ind_q_i];
        let kmer2q = &strand_aln_nodes[ind_q_ip1];

        //Wrapped around 
        let q_adjacent;
        if kmer2q.order < kmer1q.order{
            q_adjacent = (kmer2q.order == 0) && (kmer1q.order as usize == strand_aln_nodes.len() - 1)
        }
        else{
            q_adjacent = kmer2q.order == kmer1q.order + 1;
        }
        let r_adjacent = kmer1r.child_nodes.contains(&kmer2r.id);

        let kmer1rorder = kmer1r.order;
        kmer1r.color |= 1;
        kmer2r.color |= 1;

//        if kmer1r.id == 274529 || kmer2r.id == 274529{
//            dbg!(&kmer1r,&kmer2r,&kmer1q,&kmer2q);
//        }

        //        dbg!(anchors[i], anchors[i+1], &kmer1r, &kmer2r, &kmer1q, &kmer2q);

        if kmer1r.actual_ref_positions.len() > 0 && i == 0{
            kmer1r
                .actual_ref_positions
                .push(kmer1q.actual_ref_positions[0]);
        }
        if kmer2r.actual_ref_positions.len() > 0{
            kmer2r
                .actual_ref_positions
                .push(kmer2q.actual_ref_positions[0]);
        }

        if q_adjacent && r_adjacent {
            let genome_dist_query_adj;
            if forward_strand {
                genome_dist_query_adj = kmer1q.child_edge_distance[0];
            } else {
                genome_dist_query_adj = kmer2q.child_edge_distance[0];
                //                dbg!(kmer2q.child_edge_distance[0], kmer1q.child_edge_distance[0]);
                //                dbg!(&kmer1q,&kmer2q);
                //                dbg!(&kmer1r,&kmer2r);
                //                dbg!(kmer1q.kmer.to_string(),kmer2q.kmer.to_string(), kmer1q.kmer.rc().to_string(), kmer2q.kmer.rc().to_string());
                //                dbg!(&kmer1r.child_edge_distance);
                //                panic!();
            }
            let mut contains_dist = false;
            let mut edge_id = u8::MAX;
            for edge in kmer1r.child_edge_distance.iter_mut() {
                if kmer1r.child_nodes[edge.1 .1 as usize] == kmer2r.id {
                    edge_id = edge.1 .1;
                    if edge.0 == genome_dist_query_adj.0 {
                        edge.1 .0 |= 1;
                        contains_dist = true;
                        break;
                    }
                }
            }

            if !contains_dist {
                num_contains_dist_not += 1;
//                println!(
//                    "Unequal distance between adj minimizers: {:?}, {:?}, {}, {}, {}, {}, {}",
//                    genome_dist_query_adj,
//                    kmer1r.child_edge_distance,
//                    kmer1r.kmer.rc().to_string(),
//                    kmer2r.kmer.to_string(),
//                    kmer1q.kmer.to_string(),
//                    kmer2q.kmer.rc().to_string(),
//                    num_contains_dist_not
//                );
                if edge_id == u8::MAX {}
                kmer1r
                    .child_edge_distance
                    .push((genome_dist_query_adj.0, (1, edge_id)));
            }
            continue;
        } else {
            //            dbg!(&kmer2q, &kmer1q);
            //            dbg!(&kmer2r, &kmer1r);
            let mut parent_node = kmer1r;
            let mut nn_len = new_nodes.len();
            let mut range = vec![];
            if kmer2q.order < kmer1q.order {
                for node_id in kmer1q.order+1 .. strand_aln_nodes.len() as u32{
                    range.push(node_id);
                }
                for node_id in 0..kmer2q.order{
                    range.push(node_id);
                }
            }
            else{
                for node_id in kmer1q.order+1 .. kmer2q.order{
                    range.push(node_id);
                }
            }
            for i in range{
                //                dbg!(&strand_aln_nodes[i as usize]);
                let new_id = ref_node_len + nn_len;
                let mut new_kmer_node = KmerNode {
                    id: new_id as u32,
                    order: kmer1rorder,
                    kmer: strand_aln_nodes[i as usize].kmer,
                    child_nodes: SmallVec::<[u32; 1]>::new(),
                    child_edge_distance: SmallVec::<[(u8, (u64, u8)); 1]>::new(),
                    color: 1,
                    //xnor hack. truth table is
                    //11 1
                    //10 0
                    //01 0
                    //00 1
                    canonical: strand_aln_nodes[i as usize].canonical == forward_strand,
                    actual_ref_positions: SmallVec::<[usize; 0]>::new(),
                };

                if i as usize % samp_freq == 0 {
                    new_kmer_node
                        .actual_ref_positions
                        .push(strand_aln_nodes[i as usize].actual_ref_positions[0]);
                }

                

                let genome_dist_query;
                if forward_strand {
                    if i != 0{
                        genome_dist_query = strand_aln_nodes[(i - 1) as usize].child_edge_distance[0];
                    }
                    else{
                        genome_dist_query = strand_aln_nodes[strand_aln_nodes.len()-1].child_edge_distance[0];
                    }
                } else {
                    genome_dist_query = strand_aln_nodes[i as usize].child_edge_distance[0];
                }
                parent_node.child_nodes.push(new_id as u32);
                parent_node.child_edge_distance.push((
                    genome_dist_query.0,
                    (1, (parent_node.child_nodes.len() - 1) as u8),
                ));
                new_nodes.push(new_kmer_node);
                nn_len += 1;
                parent_node = &mut new_nodes[nn_len - 1];
            }
            parent_node.child_nodes.push(kmer2r.id);
            let genome_dist_query;
            if forward_strand {
                if kmer2q.order == 0 {
                    genome_dist_query =
                        strand_aln_nodes[strand_aln_nodes.len() - 1].child_edge_distance[0];
                } else {
                    genome_dist_query =
                        strand_aln_nodes[(kmer2q.order - 1) as usize].child_edge_distance[0];
                }
            } else {
                genome_dist_query =
                    strand_aln_nodes[(kmer2q.order) as usize].child_edge_distance[0];
            }
            parent_node.child_edge_distance.push((
                genome_dist_query.0,
                (1, (parent_node.child_nodes.len() - 1) as u8),
            ));
            if !has_unique_elements(&parent_node.child_nodes){
                dbg!(&parent_node,&kmer2r, &kmer1q, &kmer2q);
                dbg!(&q_adjacent,&r_adjacent);
                panic!()
            }
        }
    }

    for node in new_nodes {
        ref_nodes.push(node);
    }
}
