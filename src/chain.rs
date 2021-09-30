use crate::avl_tree::SearchTree;
use crate::data_structs::KmerNode;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use smallvec::SmallVec;
use std::mem;
use std::time::Instant;

pub fn position_max_f64(slice: &[f64]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value0.partial_cmp(value1).unwrap())
        .map(|(idx, _)| idx)
}

fn score(ref_order_dist: f64, query_order_dist: f64) -> f64 {
    let dist_multiplier = 3.0;
    let dist_ref = 16.0 - dist_multiplier * ref_order_dist;
    let dist_query = 16.0 - dist_multiplier * query_order_dist;
    let mut score = f64::min(dist_ref, dist_query);
    //try max chain length score
    score = 1.0;

    //    if score < 0.0{
    //        score = 0.05 * score - score.abs().log2();
    //    }
    //    let num_bases_query = f64::max(dist_query, 0.0);
    //    let num_bases_ref = f64::max(dist_ref, 0.0);
    //    let mut overlap = f64::min(num_bases_query, num_bases_ref);

    //    let score;
    //    if overlap < 0.001 {
    //        overlap = f64::max(dist_ref, dist_query);
    //        //        score = (f64::min(dist_query, dist_ref) / 10.0).exp();
    //        //        score = f64::min(dist_query, dist_ref).exp();
    //        score = 0.0;
    //    } else {
    //        score = overlap;
    //    }
    //
    return score;
}

fn alpha(j: usize, i: usize, anchors: &Vec<(&KmerNode, &KmerNode)>) -> f64 {
    let dist_multiplier = 1.0;
    let num_bases_query = f64::max(
        16.0 - dist_multiplier * (anchors[i].0.order as f64 - anchors[j].0.order as f64),
        0.0,
    );
    let num_bases_ref = f64::max(
        16.0 - dist_multiplier * (anchors[i].1.order as f64 - anchors[j].1.order as f64),
        0.0,
    );
    return f64::min(num_bases_query, num_bases_ref);
}

fn beta(ref_order_dist: f64, query_order_dist: f64) -> f64 {
    let dist_multiplier = 5.0;

    //    if f64::max(ref_diff, query_diff) > g {
    //        return f64::MAX;
    //    }

    let diff = (query_order_dist - ref_order_dist) * dist_multiplier;
    if diff as usize > 1000000 {
        return f64::MAX;
    }
    if diff == 0.0 {
        return 0.0;
    } else {
        //        return (0.01 * 16.0 as f64 * diff) + 0.5 * diff.log2();
        //        return diff.log2()/100.0;
        return 0.0;
    }
}

pub fn chain_seeds<'a>(
    seeds_ref: &'a mut Vec<KmerNode>,
    seeds_q: &'a Vec<KmerNode>,
    h: usize,
    chain_heuristic: bool
) -> Vec<(u32, u32)> {
    let n = 1;

    let now = Instant::now();
    let mut mini_hash_map1 = FxHashMap::default();
    for (_i, kmer_node) in seeds_ref.iter().enumerate() {
        let kmer = &kmer_node.kmer;
        let pos_vec = mini_hash_map1.entry(kmer).or_insert(vec![]);
        pos_vec.push(kmer_node.id);
    }
    let mut mini_hash_map2 = FxHashMap::default();
    for (_i, kmer_node) in seeds_q.iter().enumerate() {
        let kmer = &kmer_node.kmer;
        let pos_vec = mini_hash_map2.entry(kmer).or_insert(vec![]);
        pos_vec.push(kmer_node.id);
    }

    let mut anchors = vec![];
    let mut avl_tree: SearchTree<[usize; 2]> = SearchTree::new();

    let mut most_repet_kmer = -1;
    for kmer in mini_hash_map1.keys() {
        if mini_hash_map2.contains_key(kmer) {
            let mut count = 0;
            let ref_positions = mini_hash_map1.get(kmer).unwrap();
            let query_positions = mini_hash_map2.get(kmer).unwrap();
            for p1 in ref_positions {
                for p2 in query_positions {
                    anchors.push((&seeds_ref[*p1 as usize], &seeds_q[*p2 as usize]));
                    count += 1
                }
            }

            if count > most_repet_kmer {
                most_repet_kmer = count.clone();
            }
        }
    }
    dbg!(most_repet_kmer);
    println!(
        "Got {} anchors in {}",
        anchors.len(),
        now.elapsed().as_secs_f32()
    );

    let now = Instant::now();
    anchors.sort_by(|a, b| a.0.order.cmp(&b.0.order));
    for (i, anchor) in anchors.iter().enumerate() {
        avl_tree.insert([anchor.1.order as usize, anchor.1.id as usize]);
    }
    println!("Sorting anchors {}", now.elapsed().as_secs_f32());
    //dbg!(anchors[1],anchors[2],anchors[3]);
    //dbg!(alpha(1,2,&anchors,k),alpha(2,3,&anchors,k));
    //dbg!(beta(1,2,&anchors,k,g),beta(2,3,&anchors,k,g));
    //chaining
    let mut initial_set = FxHashSet::default();
    initial_set.insert(0);

    let mut f = vec![0.0 as f64];
    let mut pointer_array = vec![];
    for i in 0..anchors.len() {
        pointer_array.push(i);
    }

    let mut last_best_j = usize::MAX;
    avl_tree.update_query_info([anchors[0].1.order as usize, anchors[0].1.id as usize],0.0, 0, anchors[0].1.id as usize);
    for i in 1..anchors.len() {
        let mut best_f_i = usize::MIN as f64;
        let mut best_j = usize::MAX;
        let mut start = 0;
        if chain_heuristic {
            if i > h {
                start = i - h
            }
            for mut j in start..i {
                if start != 0 && j == start && last_best_j != usize::MAX {
                    j = last_best_j;
                }
                let mut incompat_chain = false;

                if anchors[j].0.order >= anchors[i].0.order {
                    incompat_chain = true;
                }

                if anchors[j].1.order >= anchors[i].1.order {
                    incompat_chain = true;
                }

                let f_cand_i;

                if incompat_chain {
                    f_cand_i = f64::MIN;
                } else {
                    let ref_order_dist = (anchors[i].0.order - anchors[j].0.order) as f64;
                    let query_order_dist = (anchors[i].1.order - anchors[j].1.order) as f64;
                    if ref_order_dist > n as f64 {
                        f_cand_i = f[j] + score(ref_order_dist, query_order_dist)
                            - beta(ref_order_dist, query_order_dist);
                    } else {
                        let graph_dist = graph_dist(&anchors[i].0, &anchors[j].0, n, &seeds_ref);
                        f_cand_i = f[j] + score(graph_dist, query_order_dist)
                            - beta(graph_dist, query_order_dist);
                    }
                }

                if f_cand_i > best_f_i {
                    best_f_i = f_cand_i;
                    best_j = j;
                }
            }
        }
        else{
            let gap_start;
            if anchors[i].1.order > 100000{
                gap_start = anchors[i].1.order - 100000;
            }
            else{
                gap_start = 0;
            }
            let (best_score, best_id) = avl_tree.mrq([gap_start as usize,0], [anchors[i].1.order as usize, anchors[i].1.id as usize], anchors[i].0.id as usize);

            if best_score == i64::MIN{
                best_f_i = 0.0;
                best_j = i;
            }
            else{
                best_f_i = best_score as f64 + 1.0;
    //            dbg!(best_id,anchors[i].1.order);
                best_j = best_id;
            }
            if anchors[i].1.order < anchors[best_j].1.order{
                dbg!(anchors[i], anchors[best_j]);
                panic!()
            }

        }

        last_best_j = best_j;
        f.push(best_f_i);
        avl_tree.update_query_info([anchors[i].1.order as usize, anchors[i].1.id as usize], best_f_i, i, anchors[i].0.id as usize);
        if best_j != usize::MAX {
            pointer_array[i] = best_j;
        }
    }

    let best_i = position_max_f64(&f).unwrap();
    let mut chain_sequence = vec![];
    let mut curr_i = best_i;
    let mut prev_i = pointer_array[curr_i];
    chain_sequence.push(curr_i);
    while curr_i != prev_i {
        chain_sequence.push(prev_i);
        curr_i = prev_i;
        prev_i = pointer_array[curr_i];
    }

    let mut best_seq_anchors = vec![];
    for i in (0..chain_sequence.len()).rev() {
        //       dbg!(anchors[chain_sequence[i]], pos1[anchors[chain_sequence[i]].0.order]);
        best_seq_anchors.push((
            anchors[chain_sequence[i]].0.id,
            anchors[chain_sequence[i]].1.id,
        ));
    }

    dbg!(f[best_i]);

    dbg!(
        (anchors[chain_sequence[0]].1.order)
            - (anchors[chain_sequence[chain_sequence.len() - 1]].1.order),
    );
    dbg!(
        (anchors[chain_sequence[0]].0.order)
            - (anchors[chain_sequence[chain_sequence.len() - 1]].0.order),
    );
    return best_seq_anchors;
}

pub fn add_align_to_graph(
    ref_nodes: &mut Vec<KmerNode>,
    aln_nodes: Vec<KmerNode>,
    anchors: Vec<(u32, u32)>,
) {
    let mut new_nodes = vec![];
    for node in ref_nodes.iter_mut() {
        node.color = node.color << 1;
    }
    for i in 0..anchors.len() - 1 {
        let ref_node_len = ref_nodes.len();

        let kmer1r;
        let kmer2r;
        let largest_anchor_id;

        if anchors[i].0 > anchors[i + 1].0 {
            largest_anchor_id = anchors[i].0;
            let (left, right) = ref_nodes.split_at_mut(largest_anchor_id as usize);
            kmer1r = &mut right[0];
            kmer2r = &mut left[anchors[i + 1].0 as usize];
        } else {
            largest_anchor_id = anchors[i + 1].0;
            let (left, right) = ref_nodes.split_at_mut(largest_anchor_id as usize);
            kmer1r = &mut left[anchors[i].0 as usize];
            kmer2r = &mut right[0];
        }

        let kmer1q = &aln_nodes[anchors[i].1 as usize];
        let kmer2q = &aln_nodes[anchors[i + 1].1 as usize];

        let q_adjacent = kmer2q.order == kmer1q.order + 1;
        let r_adjacent = kmer1r.child_nodes.contains(&kmer2r.id);

        kmer1r.color |= 1;
        kmer2r.color |= 1;

        if q_adjacent && r_adjacent {
            continue;
        } else {
            let mut parent_node = kmer1r;
            let mut nn_len = new_nodes.len();
            if kmer2q.order < kmer1q.order {
                dbg!(&anchors[i], &anchors[i + 1]);
                dbg!(&kmer2q,&kmer1q);
                panic!();
            }
            for i in kmer1q.order + 1..kmer2q.order {
                let new_id = ref_node_len + nn_len;
                let new_kmer_node = KmerNode {
                    id: new_id as u32,
                    order: i,
                    kmer: aln_nodes[i as usize].kmer,
                    child_nodes: SmallVec::<[u32; 1]>::new(),
                    color: 1,
                };

                parent_node.child_nodes.push(new_id as u32);
                new_nodes.push(new_kmer_node);
                nn_len += 1;
                parent_node = &mut new_nodes[nn_len - 1];
            }
            parent_node.child_nodes.push(kmer2r.id);
        }
    }

    for node in new_nodes {
        ref_nodes.push(node);
    }
}

//L ← Empty list that will contain the sorted nodes
//while exists nodes without a permanent mark do
//    select an unmarked node n
//    visit(n)
//
//function visit(node n)
//    if n has a permanent mark then
//        return
//    if n has a temporary mark then
//        stop   (not a DAG)
//
//    mark n with a temporary mark
//
//    for each node m with an edge from n to m do
//        visit(m)
//
//    remove temporary mark from n
//    mark n with a permanent mark
//    add n to head of L

pub fn top_sort(ref_nodes: &mut Vec<KmerNode>) {
    let mut nodes_to_visit = Vec::new();
    let mut stack_of_visited = Vec::new();
    stack_of_visited.push(0 as u32);
    nodes_to_visit.push(0 as u32);
    let mut visited = FxHashSet::default();
    let mut rev_sort_list = vec![];

    while nodes_to_visit.len() != 0 {
        let node = nodes_to_visit.pop().unwrap();
        visited.insert(node);
        let mut no_further = true;
        for child_id in ref_nodes[node as usize].child_nodes.iter() {
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
                        rev_sort_list.push(stack_of_visited.pop().unwrap());
                    }
                    break;
                }

                if *stack_of_visited.last().unwrap() == *nodes_to_visit.last().unwrap() {
                    break;
                }
                let sorted_node = stack_of_visited.pop().unwrap();
                rev_sort_list.push(sorted_node);
            }
        }
    }

    for (i, id) in rev_sort_list.iter().rev().enumerate() {
        ref_nodes[(*id) as usize].order = i as u32;
    }
}

pub fn graph_dist(
    node_forward: &KmerNode,
    node_back: &KmerNode,
    n: usize,
    all_nodes: &Vec<KmerNode>,
) -> f64 {
    let order_dist = node_forward.order - node_back.order;
    let mut visited = FxHashSet::default();

    let mut nodes_to_visit = SmallVec::<[u32; 20]>::new();
    let mut nodes_visited = SmallVec::<[u32; 20]>::new();
    nodes_to_visit.push(node_back.id);
    for i in 0..usize::min(order_dist as usize, n) {
        for node in nodes_to_visit.iter() {
            for child in all_nodes[*node as usize].child_nodes.iter() {
                nodes_visited.push(*child);
                visited.insert(child);
            }

            if visited.contains(&node_forward.id) {
                return (i + 1) as f64;
            }
        }
        mem::swap(&mut nodes_to_visit, &mut nodes_visited);
        nodes_visited.clear();
    }

    return f64::MAX;
}
