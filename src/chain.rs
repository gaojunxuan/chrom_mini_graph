use crate::avl_tree::SearchTree;
use crate::data_structs::KmerNode;
use debruijn::kmer::Kmer16;
use debruijn::Kmer;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use smallvec::SmallVec;
use std::mem;
use std::time::Instant;

pub fn print_as_binary(color: u64, string: String) {
    let binary_color = format!("{:#08b}", color);
    println!("{},{}", string, binary_color);
}

//                    let parent_color= format!("{:#08b}", parent_path_color);

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

fn score(ref_order_dist: f64, query_order_dist: f64) -> f64 {
    let dist_ref = 100.0 - 1.0 * ref_order_dist;
    let dist_query = 100.0 - 3.0 * query_order_dist;
    let score = dist_ref + dist_query;
    //try max chain length score

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

fn _alpha(j: usize, i: usize, anchors: &Vec<(&KmerNode, &KmerNode)>) -> f64 {
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
    seeds_q: &'a mut Vec<KmerNode>,
    ref_hash_map: &'a FxHashMap<Kmer16, Vec<u32>>,
    q_hash_map: &'a FxHashMap<Kmer16, Vec<u32>>,
    h: usize,
    chain_heuristic: bool,
    chain_reads: bool,
    not_used_kmers: &FxHashSet<&Kmer16>,
    circular: bool
) -> (Vec<(u32, u32)>, f64, bool) {
    let q_len = seeds_q.len();
    let now = Instant::now();

    let mut forward_anchors = vec![];
    let mut backward_anchors = vec![];

    let mut most_repet_kmer = -1;
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
                    let node_q = &seeds_q[*p2 as usize];
                    if node_r.canonical == node_q.canonical {
                        //                        forward_anchors.push((node_r,node_q));
                        forward_anchors.push((node_r.id, node_q.id));
                    } else {
                        //                        backward_anchors.push((node_r,node_q));
                        backward_anchors.push((node_r.id, node_q.id));
                    }
                    count += 1
                }
            }

            if count > most_repet_kmer {
                most_repet_kmer = count.clone();
                if !chain_reads {
                    dbg!(&kmer);
                }
            }
        }
    }

    let num_forward_anchors = forward_anchors.len();
    let num_backward_anchors = backward_anchors.len();

    let forward_strand;
    let mut anchors;
    if num_forward_anchors > num_backward_anchors {
        forward_strand = true;
        anchors = forward_anchors;
    } else {
        forward_strand = false;
        anchors = backward_anchors;
    }

    if forward_strand == false {
        for node in seeds_q.iter_mut() {
            node.order = q_len as u32 - node.order - 1;
//            for child_id in node.child_nodes.iter_mut(){
//                *child_id = q_len as u32 - *child_id - 1;
//            }
        }
    }

    println!(
        "Num forward/back anchors {},{}",
        num_forward_anchors, num_backward_anchors
    );

    if !chain_reads {
        dbg!(most_repet_kmer);
    }
    if !chain_reads {
        println!(
            "Got {} anchors in {}",
            anchors.len(),
            now.elapsed().as_secs_f32()
        );
    }

    if anchors.len() == 0 {
        return (vec![], 0.0, true);
    }

    //dbg!(anchors[1],anchors[2],anchors[3]);
    //dbg!(alpha(1,2,&anchors,k),alpha(2,3,&anchors,k));
    //dbg!(beta(1,2,&anchors,k,g),beta(2,3,&anchors,k,g));
    //chaining

    let mut f = vec![0.0 as f64];
    let mut pointer_array = vec![];
    for i in 0..anchors.len() {
        pointer_array.push(i);
    }

    anchors.sort_by(|a, b| {
        seeds_ref[a.0 as usize]
            .order
            .cmp(&seeds_ref[b.0 as usize].order)
    });

    get_chain_from_anchors(
        &mut f,
        &mut pointer_array,
        &mut anchors,
        chain_heuristic,
        chain_reads,
        h,
        &seeds_ref,
        &seeds_q,
        (0, 0),
    );
    let mut aln_score = 0.0;
    aln_score += f[position_max_f64(&f).unwrap()];

    let (mut best_seq_anchors_1, range_ref, range_query) =
        get_best_chain(f, pointer_array, &anchors, &seeds_ref, &seeds_q);

    let mut second_round_anchors = vec![];
    for anchor in anchors {
        if seeds_ref[anchor.0 as usize].order < range_ref.0
            || seeds_ref[anchor.0 as usize].order > range_ref.1
        {
            if seeds_q[anchor.1 as usize].order < range_query.0
                || seeds_q[anchor.1 as usize].order > range_query.1
            {
                second_round_anchors.push(anchor);
            }
        }
    }

    let mut f = vec![0.0 as f64];
    let mut pointer_array = vec![];
    for i in 0..second_round_anchors.len() {
        pointer_array.push(i);
    }

    if second_round_anchors.len() == 0 || chain_reads || !circular{
        return (best_seq_anchors_1, aln_score, forward_strand);
    }

    get_chain_from_anchors(
        &mut f,
        &mut pointer_array,
        &mut second_round_anchors,
        chain_heuristic,
        chain_reads,
        h,
        &seeds_ref,
        &seeds_q,
        (range_ref.1, range_query.1),
    );

    aln_score += f[position_max_f64(&f).unwrap()];
    let (mut best_seq_anchors_2, _range_ref, _range_query) = get_best_chain(
        f,
        pointer_array,
        &second_round_anchors,
        &seeds_ref,
        &seeds_q,
    );

    best_seq_anchors_1.append(&mut best_seq_anchors_2);
    //    best_seq_anchors_1.sort_by(|a, b| a.0.cmp(&b.0));
    return (best_seq_anchors_1, aln_score, forward_strand);
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

pub fn get_best_path_from_chain2(
    anchors: Vec<(u32, u32)>,
    ref_nodes: &Vec<KmerNode>,
    order_to_id: &Vec<u32>,
    qlen: u32
) -> (u64, Vec<(u32, u32)>) {
    let mut in_edges_dict: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    let mut best_paths: FxHashMap<u32, Vec<(u64, f64, usize, usize)>> = FxHashMap::default();
    let mut current_anchor_id = 0;

    if anchors.len() == 0 {
        return (1, vec![]);
    }

    let last_node = &ref_nodes[anchors.last().unwrap().0 as usize];
    let first_node = &ref_nodes[anchors[0].0 as usize];

    //Create a starting path for the first node
    for child_id in first_node.child_nodes.iter() {
        let parent_vec = in_edges_dict.entry(*child_id).or_insert(vec![]);
        parent_vec.push(first_node.id);
    }

    best_paths.insert(first_node.id, vec![(first_node.color, 10.0, 0, 0)]);
    current_anchor_id += 1;

    for i in first_node.order + 1..last_node.order + 1 {
        let id_of_node = order_to_id[i as usize];
        //        println!("Id/order of node {},{}", id_of_node, i);

        let intermediate_node = &ref_nodes[id_of_node as usize];

        for child_id in intermediate_node.child_nodes.iter() {
            let parent_vec = in_edges_dict.entry(*child_id).or_insert(vec![]);
            parent_vec.push(intermediate_node.id);
        }

        //Reevaluate paths
        let mut best_node_paths = vec![];
        best_node_paths.reserve(20);

        let parent_vec;
        if !in_edges_dict.contains_key(&intermediate_node.id) {
            best_paths.insert(
                intermediate_node.id,
                vec![(intermediate_node.color, 0.0, 0, current_anchor_id)],
            );
            parent_vec = vec![];
        } else {
            parent_vec = in_edges_dict.get(&intermediate_node.id).unwrap().to_vec();
        }

        let anchor_hit;
        let mut unreachable_past = false;

        if anchors[current_anchor_id].0 == intermediate_node.id {
            anchor_hit = true;
        } else {
            anchor_hit = false;
            if ref_nodes[anchors[current_anchor_id].0 as usize].order < intermediate_node.order {
                unreachable_past = true;
            }
        }

        //Weird stuff happens w.r.t coloring if there is an insertion in the graph
        //Consider    o
        //          /   \
        //        o   -  o
        //
        //        if the bottom-left is 111, top is 110, bottom-right is 111, then the path
        //        from BL to BR should be 001.

        let mut all_colors = 0;
        for parent_id in parent_vec.iter() {
            

            //Get the color of the edge from the parent to the current node
            let parent_node = &ref_nodes[*parent_id as usize];
            let mut edge_color = 0;
            for edge in parent_node.child_edge_distance.iter(){
                let ind = edge.1.1;
                if parent_node.child_nodes[ind as usize] == intermediate_node.id{
                    edge_color |= edge.1.0;
                }
            }

            

            let parent_paths = best_paths.get(parent_id).unwrap();
            for parent_path in parent_paths {
                let parent_path_color = parent_path.0;
                all_colors |= parent_path_color;
                //Color coherence is needed for new paths w.r.t intermediate node
                if parent_path_color & edge_color != 0 {
                    let new_color = parent_path_color & edge_color;
                    //If the intermediate node is one of the anchor nodes
                    if anchor_hit {
                        //Calculate gap cost
                        //circular stuff
                        //
                        let query_dist;
                        //                        if anchors[current_anchor_id].1 < anchors[current_anchor_id - 1].1 {
                        //                            query_dist = anchors[current_anchor_id].1 + q_len
                        //                                - anchors[current_anchor_id - 1].1;
                        //                        } else {
                        query_dist = anchors[current_anchor_id].1 as i64
                            - anchors[current_anchor_id - 1].1 as i64;
                        //                        }
                        let ref_dist = parent_path.2 + 1;
                        let gap_cost =  ((query_dist as i64).abs() - ref_dist as i64).abs().pow(2)/2 + (query_dist as i64 + ref_dist as i64)/3;

                        let new_score_add = (10 - gap_cost) as f64;
                        let updated_path = (new_color, parent_path.1 + new_score_add, 0, parent_path.3);
                        best_node_paths.push(updated_path);
                    }
                    //Just update the paths by increasing the distance by 1
                    else {
                        let updated_path = (new_color, parent_path.1, parent_path.2 + 1, parent_path.3);
                        best_node_paths.push(updated_path);
                    }
                }
            }
        }

        if intermediate_node.color & all_colors != intermediate_node.color && parent_vec.len() > 0 {
            let remaining_colors = intermediate_node.color ^ all_colors;
            let mut lowest_score = 0.0;
            for path in best_node_paths.iter(){
                if lowest_score < path.1{
                    lowest_score = path.1;
                }
            }
            if remaining_colors > 0 {
                //                println!("{:?}", &best_node_paths);
                //                print_as_binary(remaining_colors,"".to_string());
                //                print_as_binary(intermediate_node.color,"".to_string());
                //                print_as_binary(all_colors,"".to_string());
                best_node_paths.push((remaining_colors, lowest_score, 0, current_anchor_id));
            }
        }
//        if intermediate_node.actual_ref_positions.len() > 0{
//            dbg!(&best_node_paths, &intermediate_node);
//        }

        let mut color_set = FxHashSet::default();
        for path in best_node_paths.iter(){
            if color_set.contains(&path.0){
            dbg!(&best_node_paths, &intermediate_node);
            for id in in_edges_dict.get(&intermediate_node.id).unwrap(){
                dbg!(&ref_nodes[*id as usize]);
            }
            panic!();
            }
            else{
                color_set.insert(path.0);
            }
        }
//        if intermediate_node.id == 274837{
//            dbg!(&best_node_paths, &intermediate_node);
//        }
        best_paths.insert(intermediate_node.id, best_node_paths);
        if anchor_hit || unreachable_past {
            current_anchor_id += 1;
        }

        
    }

    let mut best_path_color = u64::MAX;
    let mut best_path_score = f64::MIN;
    let mut best_path_start_anchor = 0;
    let best_path = best_paths.get(&last_node.id);
    if let None = best_path {
        return (u64::MAX, vec![(0, 0)]);
    }
    for path in best_path.unwrap() {
        if path.1 > best_path_score {
            best_path_color = path.0;
            best_path_score = path.1;
            best_path_start_anchor = path.3;
        }
    }

    if best_path_score < -5000.0{
        return (u64::MAX, vec![(0, 0)]);
    }

    let mut consistent_color_anchors = vec![];
    for i in best_path_start_anchor..anchors.len(){
//    for i in 0..anchors.len(){
        let anchor = anchors[i];
        if &ref_nodes[anchor.0 as usize].color & best_path_color == best_path_color {
            consistent_color_anchors.push(anchor);
        }
    }

    println!("{:?}", best_paths.get(&last_node.id));
    return (best_path_color, consistent_color_anchors);

}

pub fn get_best_path_from_chain(
    anchors: Vec<(u32, u32)>,
    ref_nodes: &Vec<KmerNode>,
    order_to_id: &Vec<u32>,
    q_len: u32,
) -> (u64, Vec<(u32, u32)>) {
    let mut in_edges_dict: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    let mut nodes_to_search = FxHashSet::default();
    let mut best_paths: FxHashMap<u32, Vec<(u64, f64, usize)>> = FxHashMap::default();
    let n = 25;
    let mut current_anchor_id = 0;

    if anchors.len() == 0 {
        return (1, vec![]);
    }

    let last_node = &ref_nodes[anchors.last().unwrap().0 as usize];
    let first_node = &ref_nodes[anchors[0].0 as usize];

    //Create a starting path for the first node
    for child_id in first_node.child_nodes.iter() {
        nodes_to_search.insert(child_id);
        let parent_vec = in_edges_dict.entry(*child_id).or_insert(vec![]);
        parent_vec.push(first_node.id);
    }

    best_paths.insert(first_node.id, vec![(first_node.color, 10.0, 0)]);
    current_anchor_id += 1;

    for i in first_node.order + 1..last_node.order + 1 {
        let id_of_node = order_to_id[i as usize];
        //        println!("Id/order of node {},{}", id_of_node, i);

        //Don't iterate over nodes which do not have an ancestor
        if !nodes_to_search.contains(&id_of_node) {
            if i == last_node.order {
                //Edge case needs to be sorted out
                //                dbg!(&first_node,&last_node);
                //                dbg!(&ref_nodes[first_node.child_nodes[0] as usize]);
                //                dbg!(&anchors);
            }
            continue;
        }

        let intermediate_node = &ref_nodes[id_of_node as usize];

        for child_id in intermediate_node.child_nodes.iter() {
            nodes_to_search.insert(child_id);
            let parent_vec = in_edges_dict.entry(*child_id).or_insert(vec![]);
            parent_vec.push(intermediate_node.id);
        }

        //Reevaluate paths
        let mut best_node_paths = vec![];
        best_node_paths.reserve(20);
        let parent_vec = in_edges_dict.get(&intermediate_node.id).unwrap();
        let anchor_hit;
        let mut unreachable_past = false;

        if anchors[current_anchor_id].0 == intermediate_node.id {
            anchor_hit = true;
        } else {
            anchor_hit = false;
            if ref_nodes[anchors[current_anchor_id].0 as usize].order < intermediate_node.order {
                unreachable_past = true;
            }
        }

        //Weird stuff happens w.r.t coloring if there is an insertion in the graph
        //Consider    o
        //          /   \
        //        o   -  o
        //
        //        if the bottom-left is 111, top is 110, bottom-right is 111, then the path
        //        from BL to BR should be 001.

        //TODO change this when using more references.

        let mut all_colors = 0;
        for parent_id in parent_vec.iter() {
            let parent_paths = best_paths.get(parent_id).unwrap();
            for parent_path in parent_paths {
                let parent_path_color = parent_path.0;
                all_colors |= parent_path_color;
                //Color coherence is needed for new paths w.r.t intermediate node
                if parent_path_color & intermediate_node.color != 0 {
                    let new_color = parent_path_color & intermediate_node.color;
                    //If the intermediate node is one of the anchor nodes
                    if anchor_hit {
                        //Calculate gap cost
                        //circular stuff
                        //
                        let query_dist;
                        //                        if anchors[current_anchor_id].1 < anchors[current_anchor_id - 1].1 {
                        //                            query_dist = anchors[current_anchor_id].1 + q_len
                        //                                - anchors[current_anchor_id - 1].1;
                        //                        } else {
                        query_dist = anchors[current_anchor_id].1 as i64
                            - anchors[current_anchor_id - 1].1 as i64;
                        //                        }
                        let ref_dist = parent_path.2 + 1;
                        let gap_cost = 1 * ((query_dist as i64).abs() - ref_dist as i64).abs();

                        let new_score_add = (10 - gap_cost) as f64;
                        let updated_path = (new_color, parent_path.1 + new_score_add, 0);
                        best_node_paths.push(updated_path);
                    }
                    //Just update the paths by increasing the distance by 1
                    else {
                        let updated_path = (new_color, parent_path.1, parent_path.2 + 1);
                        best_node_paths.push(updated_path);
                    }
                }
            }
        }

        if intermediate_node.color & all_colors != intermediate_node.color {
            let remaining_colors = intermediate_node.color ^ all_colors;
            if remaining_colors > 0 {
                //                println!("{:?}", &best_node_paths);
                //                print_as_binary(remaining_colors,"".to_string());
                //                print_as_binary(intermediate_node.color,"".to_string());
                //                print_as_binary(all_colors,"".to_string());
                best_node_paths.push((remaining_colors, 0.0, 0));
            }
        }

        //Update the path
        //the xor operation is for the weird deletion case mentioned above.
        if parent_vec.len() > 1 {
            let mut test_set = FxHashSet::default();
            let mut path_split = true;
            best_node_paths.sort_by(|a, b| b.partial_cmp(&a).unwrap());
            //            dbg!(&best_node_paths);

            for path in best_node_paths.iter_mut() {
                let x = path.0;
                if !((x & (x - 1)) == 0) {
                    //                    dbg!(x);
                    path_split = false;
                }
                //collapse paths
                if test_set.contains(&path.0) {
                    path.0 = 0;
                } else {
                    test_set.insert(path.0);
                }
            }

            if !path_split {
                let mut changes = vec![];
                for i in 0..best_node_paths.len() {
                    let path_i = &best_node_paths[i];
                    for j in i + 1..best_node_paths.len() {
                        let path_j = &best_node_paths[j];
                        if path_i.0 & path_j.0 == path_j.0 {
                            changes.push((i, path_i.0 ^ path_j.0));
                            break;
                        }
                    }
                }

                for (i, color) in changes {
                    best_node_paths[i].0 = color;
                }
                //            if best_node_paths.len() > n{
                //                best_node_paths.drain(n..);
                //                    let parent_color= format!("{:#08b}", parent_path_color);
                //                    let int_color= format!("{:#08b}", intermediate_node.color);
                //                    dbg!(parent_color, int_color);
                //            }
                //            dbg!(&best_node_paths);
            }
        }
        best_paths.insert(intermediate_node.id, best_node_paths);
        if anchor_hit || unreachable_past {
            current_anchor_id += 1;
        }
    }

    let mut best_path_color = u64::MAX;
    let mut best_path_score = f64::MIN;
    let best_path = best_paths.get(&last_node.id);
    if let None = best_path {
        return (u64::MAX, vec![(0, 0)]);
    }
    for path in best_path.unwrap() {
        if path.1 > best_path_score {
            best_path_color = path.0;
            best_path_score = path.1;
        }
    }

    let mut consistent_color_anchors = vec![];
    for anchor in anchors {
        if &ref_nodes[anchor.0 as usize].color & best_path_color == best_path_color {
            consistent_color_anchors.push(anchor);
        }
    }

    println!("{:?}", best_paths.get(&last_node.id));
    return (best_path_color, consistent_color_anchors);
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

fn get_chain_from_anchors(
    f: &mut Vec<f64>,
    pointer_array: &mut [usize],
    anchors: &Vec<(u32, u32)>,
    chain_heuristic: bool,
    chain_reads: bool,
    h: usize,
    seeds_ref: &Vec<KmerNode>,
    seeds_q: &Vec<KmerNode>,
    modulo_positions: (u32, u32),
) {
    let q_len = seeds_q.len() as u32;
    let r_len = seeds_ref.len() as u32;
    let n = 1;
    let mut avl_tree: SearchTree<[usize; 2]> = SearchTree::new();
    let now = Instant::now();
    for (i, anchor) in anchors.iter().enumerate() {
        avl_tree.insert([
            modulo_n(seeds_q[anchor.1 as usize].order, q_len, modulo_positions.1),
            i,
        ]);
    }
    if !chain_reads {
        println!("Sorting anchors {}", now.elapsed().as_secs_f32());
    }
    let w = 100.0;
    let c1 = 1.0;

    let mut last_best_j = usize::MAX;
    avl_tree.update_query_info(
        [
            modulo_n(
                seeds_q[anchors[0].1 as usize].order,
                q_len,
                modulo_positions.1,
            ),
            0,
        ],
        c1 * (modulo_n(
            seeds_q[anchors[0].1 as usize].order,
            q_len,
            modulo_positions.1,
        ) + modulo_n(
            seeds_ref[anchors[0].0 as usize].order,
            r_len,
            modulo_positions.0,
        )) as f64,
        0,
        anchors[0].0 as usize,
        anchors[0].1 as usize,
    );

    for i in 1..anchors.len() {
        let anchoriq = modulo_n(
            seeds_q[anchors[i].1 as usize].order,
            q_len,
            modulo_positions.1,
        );
        let anchorir = modulo_n(
            seeds_ref[anchors[i].0 as usize].order,
            r_len,
            modulo_positions.0,
        );
        let mut best_f_i = usize::MIN as f64;
        let mut best_j = usize::MAX;
        let mut start = 0;
        if chain_heuristic {
            if i > h {
                start = i - h
            }
            for mut j in start..i {
                let anchorjq = seeds_q[anchors[j].1 as usize].order as usize;
                let anchorjr = seeds_ref[anchors[j].0 as usize].order as usize;
                let mut anchoriq = seeds_q[anchors[i].1 as usize].order as usize;

                if anchoriq < anchorjq {
                    anchoriq = q_len as usize + anchoriq;
                }
                let anchorir = seeds_ref[anchors[i].0 as usize].order as usize;
                //                if anchorir >= anchorjr{
                //                    anchorir = r_len as usize + anchorir;
                //                }

                //                println!("anchors {}, {}, {}, {}",anchorjq,anchorjr,anchoriq,anchorir);

                if start != 0 && j == start && last_best_j != usize::MAX {
                    j = last_best_j;
                }

                let mut incompat_chain = false;

                if anchorjr >= anchorir {
                    incompat_chain = true;
                }

                if anchorjq >= anchoriq {
                    incompat_chain = true;
                }

                let f_cand_i;

                if incompat_chain {
                    f_cand_i = f64::MIN;
                } else {
                    let ref_order_dist = (anchorir - anchorjr) as f64;
                    let query_order_dist = (anchoriq - anchorjq) as f64;
                    f_cand_i = f[j] + score(ref_order_dist, query_order_dist)
                        - beta(ref_order_dist, query_order_dist);

                    //                    if ref_order_dist > n as f64 {
                    //                        f_cand_i = f[j] + score(ref_order_dist, query_order_dist)
                    //                            - beta(ref_order_dist, query_order_dist);
                    //                    } else {
                    //                        let graph_dist = graph_dist(&anchors[i].0, &anchors[j].0, n, &seeds_ref);
                    //                        f_cand_i = f[j] + score(graph_dist, query_order_dist)
                    //                            - beta(graph_dist, query_order_dist);
                    //                    }
                }
                if f_cand_i > best_f_i {
                    best_f_i = f_cand_i;
                    best_j = j;
                }
            }
        } else {
            let gap_start = 0;
            //            if anchors[i].1.order > 10000 {
            //                gap_start = anchors[i].1.order - 10000;
            //            } else {
            //                gap_start = 0;
            //            }
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
                //                best_f_i = best_score as f64 + 1.0;
                //            dbg!(best_id,anchors[i].1.order);
                best_j = best_id;
                if chain_reads {
                    best_f_i = best_score as f64 - c1 * (anchorir + anchoriq) as f64 + w;
                //                    let diff_anc_ref = anchors[i].0.order - anchors[best_j].0.order;
                //                    let diff_anc_q = anchors[i].1.order - anchors[best_j].1.order;
                ////                    best_f_i = best_score as f64 + 150.0
                ////                        - f64::max(diff_anc_ref as f64, diff_anc_q as f64);
                //                    best_f_i = best_score as f64 + 300.0 - (diff_anc_ref as f64 - diff_anc_q as f64).abs();
                } else {
                    //                    best_f_i = best_score as f64
                    //                        + f64::max(
                    //                            100.0 - (anchors[i].1.order - anchors[best_j].1.order) as f64,
                    //                            -0.01,
                    //                        );
                    best_f_i = best_score as f64 - c1 * (anchorir + anchoriq) as f64 + w;
                }

                if best_f_i < 0.0 {
                    best_f_i = 0.0;
                    best_j = i;
                }
            }
            if anchoriq
                < modulo_n(
                    seeds_q[anchors[best_j].1 as usize].order,
                    q_len,
                    modulo_positions.1,
                )
            {
                dbg!(&anchors[i], &anchors[best_j]);
                panic!()
            }
        }

        last_best_j = best_j;
        if best_f_i < 0.0 {
            last_best_j = i
        }
        f.push(best_f_i);
        avl_tree.update_query_info(
            [anchoriq, i],
            best_f_i + w + c1 * (anchorir + anchoriq) as f64,
            i,
            anchors[i].0 as usize,
            anchors[i].1 as usize,
        );
        if best_j != usize::MAX {
            pointer_array[i] = best_j;
        }
    }
}

fn get_best_chain(
    f: Vec<f64>,
    pointer_array: Vec<usize>,
    anchors: &Vec<(u32, u32)>,
    seeds_ref: &Vec<KmerNode>,
    seeds_q: &Vec<KmerNode>,
) -> (Vec<(u32, u32)>, (u32, u32), (u32, u32)) {
    let best_i = position_max_f64(&f).unwrap();
    let mut best_seq_anchors = vec![];
    let mut chain_sequence = vec![];
    let mut curr_i = best_i;
    let mut prev_i = pointer_array[curr_i];
    chain_sequence.push(curr_i);
    while curr_i != prev_i {
        chain_sequence.push(prev_i);
        curr_i = prev_i;
        prev_i = pointer_array[curr_i];
    }

    for i in (0..chain_sequence.len()).rev() {
        //       dbg!(anchors[chain_sequence[i]], pos1[anchors[chain_sequence[i]].0.order]);
        best_seq_anchors.push((anchors[chain_sequence[i]].0, anchors[chain_sequence[i]].1));
    }

    //    dbg!(f[best_i]);

    println!(
        "End/start of ref anchors: {},{}",
        seeds_ref[anchors[chain_sequence[0]].0 as usize].order,
        seeds_ref[anchors[chain_sequence[chain_sequence.len() - 1]].0 as usize].order
    );
    println!(
        "End/start of query anchors: {},{}",
        seeds_q[anchors[chain_sequence[0]].1 as usize].order,
        seeds_q[anchors[chain_sequence[chain_sequence.len() - 1]].1 as usize].order
    );
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
        seeds_q[first_anchor.1 as usize].order,
        seeds_q[last_anchor.1 as usize].order,
    );

    return (best_seq_anchors, range_ref, range_query);
}
