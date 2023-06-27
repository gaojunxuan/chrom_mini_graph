use crate::align;
use crate::graph_utils::get_closest_node;
use simple_logger::SimpleLogger;
use crate::avl_tree::SearchTree;
use crate::constants;
use crate::data_structs::{SparseMatrix};
use cmg_shared::data_structs::{KmerNode, Bubble};
use crate::data_structs::{Anchors, Color};
use debruijn::kmer::Kmer16;
use debruijn::Kmer;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use smallvec::SmallVec;
use core::num;
use std::cmp::min;
use std::i32::MAX;
use std::mem::{self, swap};
use std::time::Instant;

pub fn print_as_binary(color: Color, string: String) {
    let binary_color = format!("{:#08b}", color);
    println!("{},{}", string, binary_color);
}

pub fn get_kmer_dict_mut(seeds: &mut Vec<KmerNode>) -> FxHashMap<Kmer16, Vec<u32>> {
    let mut mini_hash_map = FxHashMap::default();
    for (_i, kmer_node) in seeds.iter().enumerate() {
        let kmer = &kmer_node.kmer;
        let pos_vec = mini_hash_map.entry(*kmer).or_insert(vec![]);
        pos_vec.push(kmer_node.id);
    }

    return mini_hash_map;
}

pub fn get_kmer_dict(seeds: &Vec<KmerNode>) -> FxHashMap<Kmer16, Vec<u32>> {
    let mut mini_hash_map = FxHashMap::default();
    for (_i, kmer_node) in seeds.iter().enumerate() {
        let kmer = &kmer_node.kmer;
        let pos_vec = mini_hash_map.entry(*kmer).or_insert(vec![]);
        pos_vec.push(kmer_node.id);
    }

    return mini_hash_map;
}

pub fn position_max_f64(slice: &[f64]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value0.partial_cmp(value1).unwrap())
        .map(|(idx, _)| idx)
}

pub fn position_max_f64_nth(slice: &[f64], n: usize) -> usize {
    let mut vec: Vec<_> = slice.iter().enumerate().collect();
    vec.sort_by(|(_, v0), (_, v1)| v1.partial_cmp(v0).unwrap());
    return vec[n].0;
}

#[inline]
fn heuristic_score(ref_order_dist: f64, query_order_dist: f64) -> f64 {
    //    let dist_ref = -2.0 * f64::sqrt(ref_order_dist);
    //    let max_dist = f64::max(ref_order_dist,query_order_dist);
    //    let dist_ref = -1.0 * ref_order_dist;
    //    let dist_query = 30.0 - f64::sqrt(max_dist);
    let mut gap_cost = (ref_order_dist - query_order_dist).abs();
    //    if query_order_dist * query_order_dist + 50. < ref_order_dist {
    //    } else {
    //        gap_cost = f64::min((ref_order_dist - query_order_dist).sqrt(), gap_cost);
    //    }
    //    let fuzzy_gap_cost = 500.;
    //    let fuzzy_gap_cost = gap_cost.powi(2) / query_order_dist;
    //    let linear_cost = f64::sqrt(ref_order_dist) + query_order_dist;
    //    let score = constants::BANDED_CHAINING_BASE_SCORE - f64::min(gap_cost, fuzzy_gap_cost);
    let score = constants::BANDED_CHAINING_BASE_SCORE - gap_cost;
    //    let score = constants::BANDED_CHAINING_BASE_SCORE - linear_cost;
    return score;
}

fn get_chains<'a>(
    seeds_ref: &'a Vec<KmerNode>,
    seeds_query: &'a mut Vec<KmerNode>,
    h: usize,
    chain_heuristic: bool,
    chain_reads: bool,
    circular: bool,
    forward_strand: bool,
    anchors: &mut Vec<(u32, u32)>,
    closest_bubble_source: Option<&Vec<Option<(u32,u32,u16)>>>,
    bubbles: Option<&Vec<Bubble>>,
    ref_nodes: Option<&Vec<KmerNode>>
) -> Vec<(Anchors, f64)> {
    let q_len = seeds_query.len();
    let order_val_last = seeds_query[q_len - 1].order_val;

    if forward_strand == false {
        for node in seeds_query.iter_mut() {
            node.order = q_len as u32 - node.order - 1;
            node.order_val = order_val_last - node.order_val;
            //            for child_id in node.child_nodes.iter_mut(){
            //                *child_id = q_len as u32 - *child_id - 1;
            //            }
        }
    }

    if anchors.len() == 0 {
        return vec![(vec![], 0.0)];
    }

    //dbg!(anchors[1],anchors[2],anchors[3]);
    //dbg!(alpha(1,2,&anchors,k),alpha(2,3,&anchors,k));
    //dbg!(beta(1,2,&anchors,k,g),beta(2,3,&anchors,k,g));
    //chaining

    let mut f = vec![0.0 as f64];
    f.reserve(anchors.len());

    let mut pointer_array = vec![];
    for i in 0..anchors.len() {
        pointer_array.push(i);
    }

    anchors.sort_by(|a, b| {
        seeds_ref[a.0 as usize]
            .order
            .cmp(&seeds_ref[b.0 as usize].order)
    });

    score_anchors(
        &mut f,
        &mut pointer_array,
        anchors,
        chain_heuristic,
        chain_reads,
        h,
        &seeds_ref,
        &seeds_query,
        (0, 0),
        closest_bubble_source,
        bubbles,
        ref_nodes,
    );

    //    let (mut best_seq_anchors_1, range_ref, range_query) =
    let mut chain_range_vec = get_best_chains(
        f,
        pointer_array,
        &anchors,
        &seeds_ref,
        &seeds_query,
        chain_reads,
    );

    //Return all secondary chains for read chaining
    if chain_reads {
        return chain_range_vec
            .iter_mut()
            .map(|x| (mem::take(&mut x.0), x.3))
            .collect();
    }
    //Return best chain for genome chaining
    else if !circular {
        return vec![(mem::take(&mut chain_range_vec[0].0), chain_range_vec[0].3)];
    } else {
        //For second round of chaining, include anchors outside the range
        //of previous chain
        let mut best_seq_anchors_1 = mem::take(&mut chain_range_vec[0].0);
        let range_ref = chain_range_vec[0].1;
        let range_query = chain_range_vec[0].2;
        let best_aln_score = chain_range_vec[0].3;
        let mut second_round_anchors = vec![];
        for anchor in anchors {
            if seeds_ref[anchor.0 as usize].order < range_ref.0
                || seeds_ref[anchor.0 as usize].order > range_ref.1
            {
                if seeds_query[anchor.1 as usize].order < range_query.0
                    || seeds_query[anchor.1 as usize].order > range_query.1
                {
                    second_round_anchors.push(*anchor);
                }
            }
        }

        let mut f = vec![0.0 as f64];
        let mut pointer_array = vec![];
        for i in 0..second_round_anchors.len() {
            pointer_array.push(i);
        }

        if second_round_anchors.len() == 0 || chain_reads || !circular {
            return vec![(best_seq_anchors_1, best_aln_score)];
        }

        score_anchors(
            &mut f,
            &mut pointer_array,
            &mut second_round_anchors,
            chain_heuristic,
            chain_reads,
            h,
            &seeds_ref,
            &seeds_query,
            (range_ref.1, range_query.1),
            closest_bubble_source,
            bubbles,
            ref_nodes,
        );

        let (mut best_seq_anchors_2, _range_ref, _range_query, second_best_aln_score) = mem::take(
            &mut get_best_chains(
                f,
                pointer_array,
                &second_round_anchors,
                &seeds_ref,
                &seeds_query,
                chain_reads,
            )[0],
        );

        best_seq_anchors_1.append(&mut best_seq_anchors_2);
        //Reverse back because we may not want q reversed if the best alignment was actually forward.
        if forward_strand == false {
            for node in seeds_query.iter_mut() {
                node.order = q_len as u32 - node.order - 1;
                //            for child_id in node.child_nodes.iter_mut(){
                //                *child_id = q_len as u32 - *child_id - 1;
                //            }
            }
        }
        //    best_seq_anchors_1.sort_by(|a, b| a.0.cmp(&b.0));
        return vec![(best_seq_anchors_1, second_best_aln_score + best_aln_score)];
    }
}

pub fn anchors_from_seeds(
    seeds_ref: &Vec<KmerNode>,
    seeds_query: &Vec<KmerNode>,
    ref_hash_map: &FxHashMap<Kmer16, Vec<u32>>,
    q_hash_map: &FxHashMap<Kmer16, Vec<u32>>,
    not_used_kmers: &FxHashSet<Kmer16>,
    only_primary: bool,
) -> (Anchors, Anchors, usize, usize) {
    let mut forward_anchors = vec![];
    let mut forward_hits_set = FxHashSet::default();
    let mut backward_anchors = vec![];
    let mut backward_hits_set = FxHashSet::default();

    for kmer in q_hash_map.keys() {
        if not_used_kmers.contains(kmer) {
            continue;
        }
        if ref_hash_map.contains_key(kmer) {
            let mut count = 0;
            let ref_positions = ref_hash_map.get(kmer).unwrap();
            let query_positions = q_hash_map.get(kmer).unwrap();
            for p1 in ref_positions {
                for p2 in query_positions {
                    let node_r = &seeds_ref[*p1 as usize];
                    let node_q = &seeds_query[*p2 as usize];
                    if node_r.primary_base.is_none() && only_primary {
                        continue;
                    }
                    if node_r.canonical == node_q.canonical {
                        //                        forward_anchors.push((node_r,node_q));
                        forward_anchors.push((node_r.id, node_q.id));
                        forward_hits_set.insert(&node_r.kmer);
                    } else {
                        //                        backward_anchors.push((node_r,node_q));
                        backward_anchors.push((node_r.id, node_q.id));
                        backward_hits_set.insert(&node_r.kmer);
                    }
                    count += 1
                }
            }
        }
    }

    //    let num_forward_anchors = forward_anchors.len();
    //    let num_backward_anchors = backward_anchors.len();
    let num_forward_anchors = forward_hits_set.len();
    let num_backward_anchors = backward_hits_set.len();

    return (
        forward_anchors,
        backward_anchors,
        num_forward_anchors,
        num_backward_anchors,
    );
}

pub fn chain_seeds<'a>(
    seeds_ref: &'a Vec<KmerNode>,
    seeds_query: &'a mut Vec<KmerNode>,
    ref_hash_map: &'a FxHashMap<Kmer16, Vec<u32>>,
    q_hash_map: &'a FxHashMap<Kmer16, Vec<u32>>,
    h: usize,
    chain_heuristic: bool,
    chain_reads: bool,
    not_used_kmers: &FxHashSet<Kmer16>,
    circular: bool,
    closest_bubble_source: Option<&Vec<Option<(u32,u32,u16)>>>,
    bubbles: Option<&Vec<Bubble>>,
    ref_nodes: Option<&Vec<KmerNode>>
) -> Vec<(Anchors, f64, bool)> {
    let (mut forward_anchors, mut backward_anchors, num_forward_anchors, num_backward_anchors) =
        anchors_from_seeds(
            seeds_ref,
            seeds_query,
            ref_hash_map,
            q_hash_map,
            not_used_kmers,
            false,
        );

    let mut forward_strand;
    let ambig;
    if num_forward_anchors as f64 > num_backward_anchors as f64 * constants::AMBIGUOUS_FRACTION {
        forward_strand = true;
        ambig = false;
    } else if num_backward_anchors as f64
        > num_forward_anchors as f64 * constants::AMBIGUOUS_FRACTION
    {
        forward_strand = false;
        ambig = false;
    } else {
        ambig = true;
        forward_strand = true;
    }

    log::trace!(
        "Num forward/back anchors {},{}, Query length {}",
        forward_anchors.len(),
        backward_anchors.len(),
        seeds_query.len()
    );

    if !ambig && !chain_reads {
        let mut anchors;
        if forward_strand {
            anchors = forward_anchors;
        } else {
            anchors = backward_anchors
        }
        let mut chains_scores = get_chains(
            seeds_ref,
            seeds_query,
            h,
            chain_heuristic,
            chain_reads,
            circular,
            forward_strand,
            &mut anchors,
            closest_bubble_source,
            bubbles,
            ref_nodes,
        );
        return chains_scores
            .iter_mut()
            .map(|x| (mem::take(&mut x.0), x.1, forward_strand))
            .collect();
    } else {
        // //        let chain_heuristic;
        // if forward_anchors.len() > 100000 || backward_anchors.len() > 100000 {
        //     //            chain_heuristic = true;
        // } else {
        //     //            chain_heuristic = false;
        // }
        forward_strand = true;
        let chains_scores_forward = get_chains(
            seeds_ref,
            seeds_query,
            h,
            chain_heuristic,
            chain_reads,
            circular,
            forward_strand,
            &mut forward_anchors,
            closest_bubble_source,
            bubbles,
            ref_nodes,
        );
        forward_strand = false;
        let chains_scores_backward = get_chains(
            seeds_ref,
            seeds_query,
            h,
            chain_heuristic,
            chain_reads,
            circular,
            forward_strand,
            &mut backward_anchors,
            closest_bubble_source,
            bubbles,
            ref_nodes,
        );
        let mut return_chains = vec![];
        for (chain, score) in chains_scores_forward.into_iter() {
            return_chains.push((chain, score, true));
        }
        for (chain, score) in chains_scores_backward.into_iter() {
            return_chains.push((chain, score, false));
        }
        return return_chains;
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

pub fn get_best_path_from_chain_rewrite(
    anchors: &Anchors,
    ref_nodes: &Vec<KmerNode>,
    order_to_id: &Vec<u32>,
    query_nodes: &Vec<KmerNode>,
    read_length: usize,
    is_primary_ref: bool,
) -> (Vec<Color>, Vec<(Anchors, f64)>) {
    let mut colour_paths: FxHashMap<usize, (f64, usize, usize, usize)> = FxHashMap::default();

    if anchors.len() == 0 {
        log::trace!("No anchors found");
        return (vec![], vec![]);
    }
    let last_node = &ref_nodes[anchors.last().unwrap().0 as usize];
    let first_node = &ref_nodes[anchors[0].0 as usize];
    let mut current_anchor_id = 0;

    for i in first_node.order..last_node.order + 1 {
        let id_of_node = order_to_id[i as usize];
        let intermediate_node = &ref_nodes[id_of_node as usize];
        let anchor_hit;
        if anchors[current_anchor_id].0 == intermediate_node.id {
            anchor_hit = true;
        } else {
            anchor_hit = false;
        }
        let colours = align::get_nonzero_bits_fast(intermediate_node.color);
        for colour in colours {
            if colour_paths.contains_key(&colour) {
                if anchor_hit {
                    let mut tup = colour_paths.get_mut(&colour).unwrap();
                    let last_anchor_seen = tup.2;
                    let query_dist;
                    query_dist = (query_nodes[anchors[current_anchor_id].1 as usize]
                        .actual_ref_positions[0] as i64
                        - query_nodes[anchors[last_anchor_seen].1 as usize].actual_ref_positions[0]
                            as i64)
                        .abs();
                    let ref_dist = tup.3;
                    let gap_cost = ((query_dist as i64).abs() - ref_dist as i64).abs() as f64;
                    let new_score_add = (constants::PATH_CHAIN_BASE_SCORE - gap_cost) as f64;
                    tup.0 += new_score_add;
                    tup.2 = current_anchor_id;
                    tup.3 = 0;
                }
            } else {
                colour_paths.insert(colour, (0.0, current_anchor_id, current_anchor_id, 0));
            }
        }

        for edge in intermediate_node.child_edge_distance.iter() {
            let edge_dist = edge.0;
            let edge_colour = edge.1 .0;
            for colour in align::get_nonzero_bits_fast(edge_colour) {
                let tup = colour_paths.get_mut(&colour);
                if tup.is_none() {
                    dbg!(
                        &colour,
                        &colour_paths,
                        &intermediate_node,
                        align::get_nonzero_bits_fast(intermediate_node.color)
                    );
                    panic!()
                }
                let tup = tup.unwrap();
                tup.3 += edge_dist as usize;
            }
        }
        if anchor_hit {
            current_anchor_id += 1;
        }
    }

    let mut best_path_colors = vec![];
    let mut best_path_scores = vec![];
    let mut best_path_start_anchors = vec![];

    if colour_paths.is_empty() {
        log::trace!("Colour paths is empty");
        return (vec![], vec![]);
    }

    let read_length = (query_nodes.iter().last().unwrap().actual_ref_positions[0] as i64
        - query_nodes[0].actual_ref_positions[0] as i64)
        .abs();

    for (_colour, tup) in colour_paths.iter_mut() {
        let last_node_q = &query_nodes[anchors.last().unwrap().1 as usize];
        let first_node_q = &query_nodes[anchors[tup.1].1 as usize];
        let total_query_dist = (last_node_q.actual_ref_positions[0] as i64
            - first_node_q.actual_ref_positions[0] as i64)
            .abs();

        let length_diff = read_length - total_query_dist;
        tup.0 -= length_diff as f64;
    }

    let mut colour_paths_vec: Vec<(usize, (f64, usize, usize, usize))> =
        colour_paths.into_iter().collect();
    colour_paths_vec.sort_by(|x, y| y.1 .0.partial_cmp(&x.1 .0).unwrap());
    let best_path_score = colour_paths_vec[0].1 .0;

    if best_path_score < (anchors.len() as f64) * constants::PATH_THRESHOLD_FRACTION
        && !is_primary_ref
    {
        log::trace!(
            "Best path {} < {} is bad",
            best_path_score,
            anchors.len() as f64 * constants::PATH_THRESHOLD_FRACTION
        );
        return (vec![], vec![]);
    } else if best_path_score < -(query_nodes.len() as f64) * constants::PATH_THRESHOLD_FRACTION {
        log::trace!(
            "Best path {} < {} is bad",
            best_path_score,
            anchors.len() as f64 * constants::PATH_THRESHOLD_FRACTION
        );
        return (vec![], vec![]);
    }

    for i in 0..usize::min(5, colour_paths_vec.len()) {
        let colour_index = colour_paths_vec[i].0;
        let tup = colour_paths_vec[i].1;
        best_path_colors.push(1 << colour_index);
        best_path_scores.push(tup.0);
        best_path_start_anchors.push(tup.1);
    }

    let mut consistent_color_anchors = vec![];
    for j in 0..best_path_start_anchors.len() {
        let mut consistent_anchors = vec![];
        let start_anchor = best_path_start_anchors[j];
        for i in start_anchor..anchors.len() {
            let anchor = anchors[i];
            if &ref_nodes[anchor.0 as usize].color & best_path_colors[j] == best_path_colors[j] {
                consistent_anchors.push(anchor);
            }
        }
        consistent_color_anchors.push((consistent_anchors, best_path_scores[j]));
    }

    for tup in colour_paths_vec.iter() {
        log::trace!("{:?}", tup);
    }

    return (best_path_colors, consistent_color_anchors);
}

fn modulo_n(value: u32, n: u32, modulo_position: u32) -> usize {
    if value >= modulo_position {
        return (value - modulo_position) as usize;
    } else {
        if modulo_position > value + n {
            dbg!(value, n, modulo_position);
            panic!()
        }
        return (value + n - modulo_position) as usize;
    }
}

/// Score anchors using the banded dynamic programming algorithm
/// 
/// Results are stored in the f array; tracebacks are stored in the pointer array
/// 
pub fn score_anchors(
    f: &mut Vec<f64>,
    pointer_array: &mut [usize],
    anchors: &Anchors,
    chain_heuristic: bool,
    chain_reads: bool,
    h: usize,
    seeds_ref: &Vec<KmerNode>,
    seeds_query: &Vec<KmerNode>,
    modulo_positions: (u32, u32),
    closest_bubble_source: Option<&Vec<Option<(u32,u32,u16)>>>,
    bubbles: Option<&Vec<Bubble>>,
    ref_nodes: Option<&Vec<KmerNode>>,
) {
    let q_len = seeds_query.len() as u32;
    let r_len = seeds_ref.len() as u32;
    let mut ref_color_vec = vec![0; anchors.len()];
    if chain_reads && chain_heuristic {
        for (i, anchor) in anchors.iter().enumerate() {
            ref_color_vec[i] = seeds_ref[anchor.0 as usize].color;
        }
    }
    //    let n = 1;
    let w;
    if chain_reads {
        w = constants::READ_CHAIN_BASE_SCORE;
    } else {
        w = constants::CONTIG_CHAIN_BASE_SCORE;
    }
    let c1 = 1.0;

    let mut avl_tree: SearchTree<[usize; 2]> = SearchTree::new();
    if !chain_reads || !chain_heuristic {
        for (i, anchor) in anchors.iter().enumerate() {
            avl_tree.insert([
                modulo_n(
                    seeds_query[anchor.1 as usize].order_val,
                    q_len,
                    modulo_positions.1,
                ),
                i,
            ]);
        }
        avl_tree.update_query_info(
            [
                modulo_n(
                    seeds_query[anchors[0].1 as usize].order_val,
                    q_len,
                    modulo_positions.1,
                ),
                0,
            ],
            c1 * (modulo_n(
                seeds_query[anchors[0].1 as usize].order_val,
                q_len,
                modulo_positions.1,
            ) + modulo_n(
                seeds_ref[anchors[0].0 as usize].order_val,
                r_len,
                modulo_positions.0,
            )) as f64,
            0,
            anchors[0].0 as usize,
            anchors[0].1 as usize,
        );
    }
    let mandatory_color = Color::MAX;
    let adj_h = (h * anchors.len()) as f64 / 12500 as f64 + 20.;
    let adj_h = usize::min(adj_h as usize, h * 2);

    log::trace!("Starting chaining");

    for i in 1..anchors.len() {
        let q_oval;
        let r_oval;
        if chain_heuristic && chain_reads {
            q_oval = seeds_query[anchors[i].1 as usize].order_val;
            r_oval = seeds_ref[anchors[i].0 as usize].order_val;
        } else {
            q_oval = seeds_query[anchors[i].1 as usize].order_val;
            r_oval = seeds_ref[anchors[i].0 as usize].order_val;
        }
        let anchoriq = modulo_n(q_oval, q_len, modulo_positions.1);
        let anchorir = modulo_n(r_oval, r_len, modulo_positions.0);
        //        let mut best_f_i = usize::MIN as f64;
        let mut best_f_i = 0. as f64;
        let mut best_j = usize::MAX;
        let mut num_iter = 0;
        let mut max_num_iter = 0;
        
        if chain_heuristic && chain_reads {
            for j in (0..i).rev() {
                if num_iter == adj_h || max_num_iter > 3 * adj_h {
                    // log::trace!("Breaking at i = {}", i);
                    break;
                }

                let anchorjq = seeds_query[anchors[j].1 as usize].order_val as usize;
                let anchorjr = seeds_ref[anchors[j].0 as usize].order_val as usize;
                let anchoriq = seeds_query[anchors[i].1 as usize].order_val as usize;

                if anchoriq < anchorjq {
                    //panic!("Don't deal with circular mappings right now");
                    //                    anchoriq = q_len as usize + anchoriq;
                }
                let anchorir = seeds_ref[anchors[i].0 as usize].order_val as usize;

                let color_history_jr = ref_color_vec[j];
                let color_ir = seeds_ref[anchors[i].0 as usize].color;
                //                if anchorir >= anchorjr{
                //                    anchorir = r_len as usize + anchorir;
                //                }

                //                println!("anchors {}, {}, {}, {}",anchorjq,anchorjr,anchoriq,anchorir);

                //                if start != 0 && j == start && last_best_j != usize::MAX {
                //                    j = last_best_j;
                //                }

                let mut incompat_chain = false;

                if anchorjr >= anchorir {
                    incompat_chain = true;
                }

                //This forces the chain to be a walkable path in the DAG
                if color_history_jr & color_ir == 0 {
                    incompat_chain = true;
                }

                if anchorjq >= anchoriq {
                    incompat_chain = true;
                }

                let mut f_cand_i;
                max_num_iter += 1;

                if incompat_chain {
                    f_cand_i = f64::MIN;
                } else {
                    num_iter += 1;
                    let ref_order_dist = (anchorir - anchorjr) as f64;
                    let query_order_dist = (anchoriq - anchorjq) as f64;
                    f_cand_i = f[j] + heuristic_score(ref_order_dist, query_order_dist);
                    if heuristic_score(ref_order_dist, query_order_dist) < 50.0 {
                        let mut dist = ref_order_dist;
                        if !closest_bubble_source.is_none() && !bubbles.is_none() {
                            //     let bubble_pos = bubble_pos.unwrap();
                            let bubbles = bubbles.unwrap();
                            let closest_bubble_source = closest_bubble_source.unwrap();
                            // negative score, check if i and j cross a superbubble
                            // log::trace!("Anchors: {}-to-{} ({},{}), ({},{}), score = {}, ref_dis = {}, query_dist = {}", 
                            //     seeds_ref[anchors[i].0 as usize].id,
                            //     seeds_ref[anchors[j].0 as usize].id,
                            //     seeds_ref[anchors[i].0 as usize].order_val,
                            //     seeds_query[anchors[i].1 as usize].order_val, 
                            //     seeds_ref[anchors[j].0 as usize].order_val,
                            //     seeds_query[anchors[j].1 as usize].order_val,
                            //     heuristic_score(ref_order_dist, query_order_dist),
                            //     ref_order_dist,
                            //     query_order_dist);

                            let i_id = seeds_ref[anchors[i].0 as usize].id;
                            let j_id = seeds_ref[anchors[j].0 as usize].id;
                            let i_pos = seeds_ref[anchors[i].0 as usize].order_val;
                            let j_pos = seeds_ref[anchors[j].0 as usize].order_val;

                            let left_most_pos = u32::min(i_pos, j_pos);
                            let right_most_pos = u32::max(i_pos, j_pos);
                            let left_most_id = if left_most_pos == i_pos { i_id } else { j_id };
                    
                            let closest_bubble = closest_bubble_source[left_most_id as usize];
                            let v_id = bubbles[closest_bubble.unwrap().1 as usize].end;
                            let u_id = bubbles[closest_bubble.unwrap().1 as usize].start;
                            
                            let longest_path_on_bubble = bubbles[closest_bubble.unwrap().1 as usize].longest_path_length.unwrap();
                            let shortest_path_on_bubble = bubbles[closest_bubble.unwrap().1 as usize].shortest_path_length.unwrap();
                            let is_unbalanced = longest_path_on_bubble > 2 * shortest_path_on_bubble;
                    
                            let u_pos = seeds_ref[u_id as usize].order_val;
                            let v_pos = seeds_ref[v_id as usize].order_val;
                    
                            // check if (i,j) contains the bubble interval (u,v)
                            if left_most_pos <= u32::min(u_pos, v_pos) && right_most_pos >= u32::max(u_pos, v_pos) && is_unbalanced {
                                // log::trace!("Anchor pair {}-{} crosses a superbubble", i_id, j_id);
                                // use alternative scoring matrix to assign a score to the anchor pair
                                let mut i_to_u = u32::MAX;
                                let mut v_to_j = u32::MAX;
                                let mut v = &seeds_ref[v_id as usize];
                                let mut u = &seeds_ref[u_id as usize];
                                // if u.order_val > v.order_val {
                                //     mem::swap(&mut u, &mut v);
                                // }
                                if u.parent_nodes.len() > 0 {
                                    // find the closest color consistent parent
                                    // for parent in u.parent_nodes.iter() {
                                    //     // check if parent is color consistent with i
                                    //     if seeds_ref[*parent as usize].color & seeds_ref[anchors[i].0 as usize].color > 0 {
                                    //         // update i_to_u
                                    //         let dist_to_parent = u.order_val - seeds_ref[*parent as usize].order_val;
                                    //         if dist_to_parent < i_to_u {
                                    //             i_to_u = dist_to_parent;
                                    //         }
                                    //     }
                                    // }
                                    let parent = u.parent_nodes[0];
                                    if seeds_ref[parent as usize].color & seeds_ref[anchors[i].0 as usize].color > 0 {
                                        // update i_to_u
                                        let dist_to_parent = u.order_val - seeds_ref[parent as usize].order_val;
                                        if dist_to_parent < i_to_u {
                                            i_to_u = dist_to_parent;
                                        }
                                    }
                                }
                                // similarly, find v_to_j
                                if v.child_nodes.len() > 0 {
                                    // for child in v.child_nodes.iter() {
                                    //     if seeds_ref[*child as usize].color & seeds_ref[anchors[j].0 as usize].color > 0 {
                                    //         let dist_to_child = seeds_ref[*child as usize].order_val - v.order_val;
                                    //         // update v_to_j
                                    //         if dist_to_child < v_to_j {
                                    //             v_to_j = dist_to_child;
                                    //         }
                                    //     }
                                    // }
                                    let child = v.child_nodes[0];
                                    if seeds_ref[child as usize].color & seeds_ref[anchors[j].0 as usize].color > 0 {
                                        let dist_to_child = seeds_ref[child as usize].order_val - v.order_val;
                                        // update v_to_j
                                        if dist_to_child < v_to_j {
                                            v_to_j = dist_to_child;
                                        }
                                    }
                                }
                                // update dist
                                if i_to_u != u32::MAX && v_to_j != u32::MAX {
                                    if !ref_nodes.is_none() {
                                        let ref_nodes = ref_nodes.unwrap();
                                        let u_pos_est = ref_nodes[u.closest_ref as usize].actual_ref_positions[0];
                                        let v_pos_est = ref_nodes[v.closest_ref as usize].actual_ref_positions[0];
                                        let u_to_v = u_pos_est.abs_diff(v_pos_est) as u16 + u.dist_to_closest_ref + v.dist_to_closest_ref;
                                        // dist = f64::min((i_to_u + u_to_v + v_to_j) as f64, dist) as f64;
                                        dist = f64::min(u_to_v as f64, dist) as f64;
                                    }
                                    
                                }
                            }
                            f_cand_i = f[j] + heuristic_score(dist, query_order_dist);
                    
                            //     // binary search for the closest bubble source u using i_pos as the search key
                            //     let mut l: i64 = 0;
                            //     let mut r: i64 = (bubble_pos.len() - 1) as i64;
                            //     while l <= r {
                            //         let m: i64 = l + (r - l) / 2;
                            //         if bubble_pos[m as usize].0 >= left_most_pos {
                            //             r = m - 1;
                            //         } else if bubble_pos[m as usize].0 < left_most_pos {
                            //             l = m + 1;
                            //         } else {
                            //             break;
                            //         }
                            //     }
                            //     // l is the index of the closest bubble source u
                            //     // third position of the tuple is the bubble id
                            //     if l < bubble_pos.len() as i64 {
                            //         let l = (usize::max(l as usize, 0));
                            //         let u = bubble_pos[l].2;
                            //         let closest_bubble = &bubbles[u];
                            //         let u_pos = seeds_ref[closest_bubble.start as usize].order_val;
                            //         let v_pos = seeds_ref[closest_bubble.end as usize].order_val;
                            //         // check if (i,j) contains the bubble interval (u,v)
                            //         if left_most_pos <= u32::min(u_pos, v_pos) && right_most_pos >= u32::max(u_pos, v_pos) {
                            //             // log::trace!("Anchor pair {}-{} crosses a superbubble", i_id, j_id);
                            //             // use alternative scoring matrix to assign a score to the anchor pair
                            //             let mut i_to_u = u32::MAX;
                            //             let mut v_to_j = u32::MAX;
                            //             let mut v = &seeds_ref[bubbles[u as usize].end as usize];
                            //             let mut u = &seeds_ref[closest_bubble.start as usize];
                            //             if u.order_val > v.order_val {
                            //                 mem::swap(&mut u, &mut v);
                            //             }
                            //             if u.parent_nodes.len() > 0 {
                            //                 // find the closest color consistent parent
                            //                 for parent in u.parent_nodes.iter() {
                            //                     // check if parent is color consistent with i
                            //                     if seeds_ref[*parent as usize].color & seeds_ref[anchors[i].0 as usize].color > 0 {
                            //                         // update i_to_u
                            //                         let dist_to_parent = u.order_val - seeds_ref[*parent as usize].order_val;
                            //                         if dist_to_parent < i_to_u {
                            //                             i_to_u = dist_to_parent;
                            //                         }
                            //                     }
                            //                 }
                            //             }
                            //             // similarly, find v_to_j
                            //             if v.child_nodes.len() > 0 {
                            //                 for child in v.child_nodes.iter() {
                            //                     if seeds_ref[*child as usize].color & seeds_ref[anchors[j].0 as usize].color > 0 {
                            //                         let dist_to_child = seeds_ref[*child as usize].order_val - v.order_val;
                            //                         if dist_to_child < v_to_j {
                            //                             v_to_j = dist_to_child;
                            //                         }
                            //                     }
                            //                 }
                            //             }
                            //             // update dist
                            //             if i_to_u != u32::MAX && v_to_j != u32::MAX {
                            //                 let mut u_to_v = min(dist_mat.get(u.id as usize, v.id as usize), dist_mat.get(v.id as usize, u.id as usize));
                            //                 if u_to_v == u32::MAX {
                            //                     u_to_v = 0;
                            //                 }
                            //                 dist = f64::min((i_to_u + u_to_v + v_to_j) as f64, dist) as f64;
                            //             }
                            //         }
                            //     }
                            // }
                            // f_cand_i = f[j] + heuristic_score(dist, query_order_dist);
                        }
                    }
                }
                
                if f_cand_i > best_f_i {
                    best_f_i = f_cand_i;
                    best_j = j;
                }
            }
    
        } else {
            let gap_start = 0;
            let (best_score, best_id) = avl_tree.mrq(
                [gap_start as usize, 0],
                [anchoriq, i],
                anchors[i].0 as usize,
                anchors[i].1 as usize,
            );

            if best_score == i64::MIN {
                best_f_i = 0.0;
                best_j = i;
            } else {
                best_j = best_id;
                if chain_reads {
                    best_f_i = best_score as f64 - c1 * (anchorir + anchoriq) as f64 + w;
                } else {
                    best_f_i = best_score as f64 - c1 * (anchorir + anchoriq) as f64 + w;
                }

                if best_f_i < 0.0 {
                    best_f_i = 0.0;
                    best_j = i;
                }
            }
            if anchoriq
                < modulo_n(
                    seeds_query[anchors[best_j].1 as usize].order_val,
                    q_len,
                    modulo_positions.1,
                )
            {
                dbg!(&anchors[i], &anchors[best_j]);
                panic!()
            }
        }

        //        last_best_j = best_j;
        if best_f_i <= 0.0 {
            best_j = i
        }
        f.push(best_f_i);
        if !chain_heuristic || !chain_reads {
            avl_tree.update_query_info(
                [anchoriq, i],
                best_f_i + w + c1 * (anchorir + anchoriq) as f64,
                i,
                anchors[i].0 as usize,
                anchors[i].1 as usize,
            );
        }
        if best_j != usize::MAX {
            pointer_array[i] = best_j;
            //TODO
            //            ref_color_vec[i] &= ref_color_vec[best_j];
        }
    }
}

/// Get the best color coherent chains from the f and pointer array
/// 
/// # Returns
/// A vector of tuples of (Anchors, (u32, u32), (u32, u32), f64). The first element is the
/// anchors in the chain; the second element is the range of the reference anchors; the third
/// element is the range of the query anchors; the fourth element is the score of the chain.
pub fn get_best_chains(
    f: Vec<f64>,
    pointer_array: Vec<usize>,
    anchors: &Anchors,
    seeds_ref: &Vec<KmerNode>,
    seeds_query: &Vec<KmerNode>,
    chain_reads: bool,
) -> Vec<(Anchors, (u32, u32), (u32, u32), f64)> {
    //    for anchor in anchors.iter(){
    //        dbg!(align::get_nonzero_bits(seeds_ref[anchor.0 as usize].color));
    //        dbg!(&seeds_ref[anchor.0 as usize]);
    //    }
    let mut vec: Vec<_> = f.iter().enumerate().collect();
    vec.sort_by(|(_, v0), (_, v1)| v1.partial_cmp(v0).unwrap());
    let mut return_chains: Vec<(Anchors, (u32, u32), (u32, u32), f64)> = vec![];
    let mut already_used_anchors = FxHashSet::default();
    let best_chain_score = f[vec[0].0];
    let cutoff_percent = constants::SECONDARY_CHAIN_CUTOFF_PERCENT;
    let cutoff_score = constants::CHAIN_CUTOFF_SCORE;
    if best_chain_score < cutoff_score {
        log::trace!("Poor chaining score {}", best_chain_score);
        return vec![];
    }
    for i in 0..vec.len() {
        let mut coherent_color = Color::MAX;
        // Only get the best chain for genome chaining
        if i > 0 && !chain_reads {
            break;
        }

        if already_used_anchors.contains(&i) {
            continue;
        }

        let mut best_seq_anchors = vec![];
        let mut chain_sequence = vec![];
        let mut curr_i = vec[i].0;
        let ith_score = f[curr_i];

        if ith_score < cutoff_percent * best_chain_score {
            break;
        }

        let mut prev_i = pointer_array[curr_i];
        chain_sequence.push(curr_i);

        while curr_i != prev_i {
            chain_sequence.push(prev_i);
            already_used_anchors.insert(curr_i);
            curr_i = prev_i;
            prev_i = pointer_array[curr_i];
        }
        already_used_anchors.insert(curr_i);

        let chain_end_ref = seeds_ref[anchors[chain_sequence[0]].0 as usize].order;
        let chain_start_ref =
            seeds_ref[anchors[chain_sequence[chain_sequence.len() - 1]].0 as usize].order;

        let mut secondary = false;
        for tup in return_chains.iter() {
            let range_ref = tup.1;
            if chain_start_ref >= range_ref.0 && chain_start_ref <= range_ref.1 {
                secondary = true
            }
            if range_ref.0 >= chain_start_ref && range_ref.0 <= chain_end_ref {
                secondary = true
            }
        }

        //TODO Don't output secondary chains for now.
        if secondary {
            //            println!("Secondary chain!")
            continue;
        }

        for i in (0..chain_sequence.len()).rev() {
            //       dbg!(anchors[chain_sequence[i]], pos1[anchors[chain_sequence[i]].0.order]);
            best_seq_anchors.push((anchors[chain_sequence[i]].0, anchors[chain_sequence[i]].1));
            coherent_color &= seeds_ref[anchors[chain_sequence[i]].0 as usize].color;
        }

        let color;
        if chain_reads && coherent_color != 0 {
            color = coherent_color
        } else {
            color = 1;
        }

        log::trace!(
            "Chain number {}, Chain score {}, Color {:?}, Length {}",
            i,
            ith_score,
            align::get_nonzero_bits(color),
            best_seq_anchors.len()
        );
        let ref_end = seeds_ref[anchors[chain_sequence[0]].0 as usize].order;
        let ref_start =
            seeds_ref[anchors[chain_sequence[chain_sequence.len() - 1]].0 as usize].order;

        log::trace!("End/start of ref anchors: {},{}", ref_end, ref_start);
        let q_end = seeds_query[anchors[chain_sequence[0]].1 as usize].order;
        let q_start = seeds_query[anchors[chain_sequence[chain_sequence.len() - 1]].1 as usize].order;

        log::trace!("End/start of query anchors: {},{}", q_end, q_start);

        if q_start >= q_end || ref_start >= ref_end {
            //            best_seq_anchors = best_seq_anchors.into_iter().rev().collect();
        }
        //    println!(
        //        "Ref/Query order anchor dist {},{}",
        //        (anchors[chain_sequence[0]].0.order)
        //            - (anchors[chain_sequence[chain_sequence.len() - 1]].0.order),
        //        (anchors[chain_sequence[0]].1.order)
        //            - (anchors[chain_sequence[chain_sequence.len() - 1]].1.order),
        //    );
        let first_anchor = &anchors[chain_sequence[chain_sequence.len() - 1]];
        let last_anchor = &anchors[chain_sequence[0]];
        let range_ref = (
            seeds_ref[first_anchor.0 as usize].order,
            seeds_ref[last_anchor.0 as usize].order,
        );
        let range_query = (
            seeds_query[first_anchor.1 as usize].order,
            seeds_query[last_anchor.1 as usize].order,
        );
        return_chains.push((best_seq_anchors, range_ref, range_query, ith_score));
    }

    return return_chains;
}
