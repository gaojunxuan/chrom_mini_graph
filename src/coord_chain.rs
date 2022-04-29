use crate::align;
use crate::avl_tree::SearchTree;
use crate::chain;
use crate::constants;
use crate::data_structs::KmerNode;
use crate::data_structs::{Anchors, Color};
use debruijn::kmer::Kmer16;
use debruijn::Kmer;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use smallvec::SmallVec;
use std::mem;
use std::time::Instant;

#[inline]
fn min_dist_ref_nodes(ref_node1: &Vec<(usize, usize)>, ref_node2: &Vec<(usize, usize)>) -> f64 {
    let mut used_colors = vec![];
    let mut color_map = vec![usize::MAX; 128];
    for (col, pos) in ref_node2 {
        used_colors.push(col);
        color_map[*col] = *pos;
    }
    for (col, pos) in ref_node1 {
        if color_map[*col] < *pos {
            color_map[*col] = *pos - color_map[*col];
        } else {
            color_map[*col] -= *pos;
        }
    }

    let use_min = true;
    if use_min {
        let mut best_dist = usize::MAX;
        for col in used_colors {
            if color_map[*col] < best_dist {
                best_dist = color_map[*col];
            }
        }
        return best_dist as f64;
    } else {
        let mut best_dist = usize::MIN;
        for col in used_colors {
            if color_map[*col] > best_dist {
                best_dist = color_map[*col];
            }
        }
        return best_dist as f64;
    }
}

#[inline]
fn min_dist_ref_nodes_lazy(
    ref_node1: &Vec<(usize, usize)>,
    ref_node2: &Vec<(usize, usize)>,
    shared_col: Color,
) -> f64 {
    let shared_col_ind = shared_col.trailing_zeros() as usize;
    let mut shared_pos = 0;
    for (col, pos) in ref_node2 {
        if *col == shared_col_ind {
            shared_pos = *pos;
            break;
        }
    }
    for (col, pos) in ref_node1 {
        if *col == shared_col_ind {
            if shared_pos < *pos {
                return (*pos - shared_pos) as f64;
            } else {
                return (shared_pos - *pos) as f64;
            }
        }
    }

    return f64::MAX;
}

#[inline]
fn heuristic_super_score(ref_order_dist: f64, query_order_dist: f64, samp_freq: usize) -> f64 {
    let gap = f64::abs(ref_order_dist - query_order_dist);
    let gap_cost;
    let dist_between_nodes_guess = samp_freq as f64 * constants::BASE_SCORE;
    let base_score = dist_between_nodes_guess * constants::ERROR.powf(constants::K as f64);
    if gap < (samp_freq as f64 * dist_between_nodes_guess) {
        gap_cost = gap;
    } else {
        gap_cost = gap;
        //        gap_cost = (gap - (samp_freq * dist_between_nodes_guess) as f64) / samp_freq as f64;
    }
    let score = base_score as f64 - gap_cost;
    //    let score = dist_between_nodes_guess as f64 - gap_cost;

    return score;
}

pub fn get_super_chains(
    seeds_ref: &Vec<KmerNode>,
    seeds_q: &Vec<KmerNode>,
    ref_hash_map: &FxHashMap<Kmer16, Vec<u32>>,
    q_hash_map: &FxHashMap<Kmer16, Vec<u32>>,
    h: usize,
    not_used_kmers: &FxHashSet<Kmer16>,
    closest_kmer_vec: &Vec<Option<u32>>,
    samp_freq: usize,
    read_length: usize,
) -> Vec<(Anchors, f64, bool)> {
    let now = Instant::now();
    let (forward_anchors, backward_anchors, num_forward_anchors, num_backward_anchors) =
        chain::anchors_from_seeds(seeds_ref, seeds_q, ref_hash_map, q_hash_map, not_used_kmers, false);

    let forward_strand;
    if num_forward_anchors > num_backward_anchors {
        forward_strand = true;
    } else if num_backward_anchors > num_forward_anchors {
        forward_strand = false;
    } else {
        forward_strand = true;
    }

    let normal_anchors;
    if forward_strand {
        normal_anchors = forward_anchors;
    } else {
        normal_anchors = backward_anchors;
    }

    let mut super_nodes_map = FxHashMap::default();
    for (ref_id, q_id) in normal_anchors.iter() {
        let opt = closest_kmer_vec[(*ref_id) as usize];
        let closest_coord_node;
        if !opt.is_none() {
            closest_coord_node = opt.unwrap();
        } else {
            continue;
        }
        //TODO can RC the vector inside the container for performance bonus
        let list_of_hits = super_nodes_map.entry(closest_coord_node).or_insert(vec![]);
        list_of_hits.push((
            seeds_q[(*q_id) as usize].actual_ref_positions[0],
            *q_id,
            ref_id,
        ));
    }

    let mut debug = vec![];
    let mut all_colors = FxHashSet::default();
    for (coord_node, list_of_hits) in super_nodes_map.iter() {
        let cutoff = (0.35 * 0.9f64.powf(16.) * samp_freq as f64) as usize;
        if list_of_hits.len() > cutoff || true {
            let ref_pos = &seeds_ref[(*coord_node) as usize].actual_ref_positions;
            let order = &seeds_ref[(*coord_node) as usize].order;
            let color = seeds_ref[(*coord_node) as usize].color;
            let list_color = align::get_nonzero_bits(color);
            let list_color: Vec<usize> = list_color.into_iter().rev().collect();
            for color_ind in list_color.iter() {
                all_colors.insert(*color_ind);
            }
            let col_ref_zip: Vec<(usize, &usize)> =
                list_color.into_iter().zip(ref_pos.into_iter()).collect();
            debug.push((order, col_ref_zip, coord_node))
        }
    }
    debug.sort_by(|x, y| x.0.cmp(&y.0));

//    dbg!(&debug);
    let mut super_anchors = vec![];
    let mut col_list = FxHashMap::default();
    for (_i, (_order, col_pos, coord_node)) in debug.iter().enumerate() {
        for (_pos, q_id, _ref_id) in super_nodes_map[coord_node].iter() {
            super_anchors.push((**coord_node, *q_id));
        }
        let mut list = vec![];
        for (col, pos) in col_pos {
            list.push((*col, **pos));
        }
        col_list.insert(*coord_node, list);
    }

    println!(
        "Super chain preprocess time {}",
        now.elapsed().as_secs_f32()
    );
    println!("Number of super anchors {}", super_anchors.len());
    let now = Instant::now();

    super_anchors.sort_by(|x, y| x.0.cmp(&y.0));

    let mut pointer_array = vec![0; super_anchors.len()];
    let mut f = vec![0.];

    score_anchors_coords(
        &mut f,
        &mut pointer_array,
        &super_anchors,
        h,
        seeds_ref,
        seeds_q,
        forward_strand,
        read_length,
        &col_list,
        samp_freq,
    );

    println!("Super chain chaining time {}", now.elapsed().as_secs_f32());

    let mut best_superchain =
        chain::get_best_chains(f, pointer_array, &super_anchors, seeds_ref, seeds_q, true);

    if best_superchain.is_empty() {
        println!("no best superchain");
        return vec![];
    } else {
        for iter in best_superchain.iter() {
            //            let (best_colors, best_list_anchors) =
            //                chain::get_best_path_from_chain2(&iter.0, &seeds_ref, &order_to_id, &seeds_q);

            let mut ref_interval = [
                seeds_ref[iter.0[0].0 as usize].order,
                seeds_ref[iter.0.last().unwrap().0 as usize].order,
            ];
            ref_interval.sort();
            for anchor in iter.0.iter() {
                eprintln!(
                    "{:?}\n{:?}",
                    col_list[&anchor.0],
                    //                seeds_ref[anchor.0 as usize].actual_ref_positions,
                    seeds_q[anchor.1 as usize].actual_ref_positions
                );
            }
            //            eprintln!(
            //                "{:?}\n{:?}",
            //                seeds_q[iter.0[0].1 as usize].a,
            //                seeds_q[iter.0.last().unwrap().1 as usize]

            //TODO
            //            let mut retrieved_anchors = vec![];
            //            let mut num_anchors_in_superchain = 0;
            //            let mut set = FxHashSet::default();
            //            for (coord_node, hits) in super_nodes_map.iter() {
            //                if seeds_ref[*coord_node as usize].order <= ref_interval[1]
            //                    && seeds_ref[*coord_node as usize].order >= ref_interval[0]
            //                {
            //                    for (_, q_id, ref_id) in hits {
            //                        retrieved_anchors.push((**ref_id, *q_id));
            //                    }
            //                    set.insert(hits);
            //                }
            //            }
            //            //            dbg!(&set);
            //            for s in set {
            //                num_anchors_in_superchain += s.len();
            //            }
            //            println!(
            //                "Number of anchors in superchain is {}",
            //                num_anchors_in_superchain
            //            );
            //            retrieved_anchors.sort();
            //            println!("Chaining induced anchors");
            //            let mut pointer_array = vec![0; super_anchors.len()];
            //            let mut f = vec![0.];
            //            chain::score_anchors(
            //                &mut f,
            //                &mut pointer_array,
            //                &retrieved_anchors,
            //                chain_heuristic,
            //                true,
            //                h,
            //                &seeds_ref,
            //                &seeds_q,
            //                (0, 0),
            //            );
            //            let mut chain_range_vec = chain::get_best_chains(
            //                f,
            //                pointer_array,
            //                &retrieved_anchors,
            //                &seeds_ref,
            //                &seeds_q,
            //                false,
            //            );
            //            if chain_range_vec.len() > 0 {}
        }

        println!(
            "Number of chains:{}, strand {}",
            best_superchain.len(),
            forward_strand
        );
        let mut best_superchain_ret = vec![];
        for i in 0..best_superchain.len() {
            let mut best_superchain_localized = vec![];
            for anchor in best_superchain[i].0.iter() {
                let candidate_anchors = &super_nodes_map[&anchor.0];
                let mut dist_ref = vec![];
                for c_anchor in candidate_anchors {
                    if c_anchor.1 == anchor.1 {
                        dist_ref.push((
                            (*c_anchor.2 as f64 - anchor.0 as f64).abs(),
                            (*c_anchor.2, c_anchor.1),
                        ));
                        //                        best_superchain_localized.push((*c_anchor.2, c_anchor.1));
                    }
                }
                best_superchain_localized.push(
                    dist_ref
                        .iter()
                        .min_by(|x, y| x.0.partial_cmp(&y.0).unwrap())
                        .unwrap()
                        .1,
                );
            }
            best_superchain_ret.push((
                best_superchain_localized,
                best_superchain[i].3,
                forward_strand,
            ));
        }
        return best_superchain_ret;
        return best_superchain
            .iter_mut()
            .map(|x| (mem::take(&mut x.0), x.3, forward_strand))
            .collect();
    }
}

pub fn score_anchors_coords(
    f: &mut Vec<f64>,
    pointer_array: &mut [usize],
    super_anchors: &Anchors,
    h: usize,
    seeds_ref: &Vec<KmerNode>,
    seeds_q: &Vec<KmerNode>,
    forward_strand: bool,
    read_length: usize,
    col_list: &FxHashMap<&u32, Vec<(usize, usize)>>,
    samp_freq: usize,
) {
    let mut interval_pointer_array: Vec<usize> = (0..pointer_array.len()).collect();
    let mut ref_color_vec = vec![0; super_anchors.len()];
    let use_interval_heuristic = true;
    for (i, anchor) in super_anchors.iter().enumerate() {
        ref_color_vec[i] = seeds_ref[(anchor.0) as usize].color;
    }
    for i in 1..super_anchors.len() {
        let mut best_f_i = 0. as f64;
        let mut best_j = usize::MAX;
        let mut num_iter = 0;
        let mut max_num_iter = 0;
        for j in (0..i).rev() {
            if num_iter == h || max_num_iter > 3 * h {
                break;
            }
            let mut incompat_chain = false;

            let anchorjr_o = seeds_ref[super_anchors[j].0 as usize].order;
            let anchorir_o = seeds_ref[super_anchors[i].0 as usize].order;
            if anchorjr_o >= anchorir_o {
                incompat_chain = true;
            }

            let color_history_jr = ref_color_vec[j];
            let color_ir = seeds_ref[super_anchors[i].0 as usize].color;

            //This forces the chain to be a walkable path in the DAG
            if color_history_jr & color_ir == 0 {
                incompat_chain = true;
            }

            let anchorjq;
            let anchoriq;
            if forward_strand {
                anchorjq = seeds_q[super_anchors[j].1 as usize].actual_ref_positions[0];
                anchoriq = seeds_q[super_anchors[i].1 as usize].actual_ref_positions[0];
            } else {
                anchorjq =
                    read_length - seeds_q[super_anchors[j].1 as usize].actual_ref_positions[0];
                anchoriq =
                    read_length - seeds_q[super_anchors[i].1 as usize].actual_ref_positions[0];
            }

            //            let anchorjr_n = &seeds_ref[super_anchors[j].0 as usize];
            //            let anchorir_n = &seeds_ref[super_anchors[i].0 as usize];

            if anchorjq >= anchoriq {
                incompat_chain = true;
            }

            let mut f_cand_i;
            let use_lazy = false;
            max_num_iter += 1;

            if incompat_chain {
                f_cand_i = f64::MIN;
            } else {
                num_iter += 1;
                //                let ref_order_dist = min_dist_ref_nodes(
                //                    &col_list[&super_anchors[j].0],
                //                    &col_list[&super_anchors[i].0],
                //                );
                let ref_order_dist;
                if use_lazy {
                    ref_order_dist = min_dist_ref_nodes_lazy(
                        &col_list[&super_anchors[j].0],
                        &col_list[&super_anchors[i].0],
                        color_history_jr & color_ir,
                    );
                } else {
                    ref_order_dist = min_dist_ref_nodes(
                        &col_list[&super_anchors[j].0],
                        &col_list[&super_anchors[i].0],
                    );
                }

                let query_order_dist = (anchoriq - anchorjq) as f64;
                let heur_score = heuristic_super_score(ref_order_dist, query_order_dist, samp_freq);

                f_cand_i = f[j] + heur_score;
                if use_interval_heuristic {
                    if interval_pointer_array[j] != j {
                        let s = interval_pointer_array[j];
                        let anchorsq;
                        if forward_strand {
                            anchorsq = seeds_q[super_anchors[s].1 as usize].actual_ref_positions[0];
                        } else {
                            anchorsq = read_length
                                - seeds_q[super_anchors[s].1 as usize].actual_ref_positions[0];
                        }
                        let ref_order_dist_s = min_dist_ref_nodes(
                            &col_list[&super_anchors[s].0],
                            &col_list[&super_anchors[i].0],
                        );
                        let query_order_dist_s = (anchoriq - anchorsq) as f64;
                        f_cand_i +=
                            heuristic_super_score(ref_order_dist_s, query_order_dist_s, samp_freq);
                    }
                }
            }
            if f_cand_i > best_f_i {
                best_f_i = f_cand_i;
                best_j = j;
            }
        }
        if best_f_i <= 0.0 {
            best_j = i;
        }
        f.push(best_f_i);
        if best_j != usize::MAX {
            pointer_array[i] = best_j;
            interval_pointer_array[i] = interval_pointer_array[best_j];
        }
    }
}

pub fn get_base_chains(
    seeds_ref: &Vec<KmerNode>,
    seeds_q: &Vec<KmerNode>,
    ref_hash_map: &FxHashMap<Kmer16, Vec<u32>>,
    q_hash_map: &FxHashMap<Kmer16, Vec<u32>>,
    h: usize,
    not_used_kmers: &FxHashSet<Kmer16>,
    read_length: usize,
) -> Vec<(Anchors, f64, bool)> {
    let now = Instant::now();
    let (forward_anchors, backward_anchors, num_forward_anchors, num_backward_anchors) =
        chain::anchors_from_seeds(seeds_ref, seeds_q, ref_hash_map, q_hash_map, not_used_kmers, true);

    let forward_strand;
    if num_forward_anchors > num_backward_anchors {
        forward_strand = true;
    } else if num_backward_anchors > num_forward_anchors {
        forward_strand = false;
    } else {
        forward_strand = true;
    }

    let mut normal_anchors;
    if forward_strand {
        normal_anchors = forward_anchors;
    } else {
        normal_anchors = backward_anchors;
    }

    normal_anchors.sort_by(|x, y| x.0.cmp(&y.0));
    let mut pointer_array = vec![0; normal_anchors.len()];
    let mut f = vec![0.];
    score_primary_anchors_coords(
        &mut f,
        &mut pointer_array,
        &normal_anchors,
        h,
        seeds_ref,
        seeds_q,
        forward_strand,
        read_length,
    );

    println!("Primary ref chain chaining time {}", now.elapsed().as_secs_f32());

    let mut best_primary_ref_chain =
        chain::get_best_chains(f, pointer_array, &normal_anchors, seeds_ref, seeds_q, true);

    if best_primary_ref_chain.is_empty() {
        println!("no best primary ref chain");
        return vec![];
    } else {
        return best_primary_ref_chain
            .iter_mut()
            .map(|x| (mem::take(&mut x.0), x.3, forward_strand))
            .collect();
    }
}

pub fn score_primary_anchors_coords(
    f: &mut Vec<f64>,
    pointer_array: &mut [usize],
    primary_ref_anchors: &Anchors,
    h: usize,
    seeds_ref: &Vec<KmerNode>,
    seeds_q: &Vec<KmerNode>,
    forward_strand: bool,
    read_length: usize,
) {
    let mut interval_pointer_array: Vec<usize> = (0..pointer_array.len()).collect();
    let adj_h = (h * primary_ref_anchors.len()) as f64 / 12500 as f64 + 20.;
    let adj_h = adj_h as usize;
    let adj_h = usize::max(adj_h as usize, h * 3);
    for i in 1..primary_ref_anchors.len() {
        let mut best_f_i = 0. as f64;
        let mut best_j = usize::MAX;
        let mut num_iter = 0;
        let mut max_num_iter = 0;
        for j in (0..i).rev() {
            if num_iter == adj_h || max_num_iter > 3 *  adj_h {
                break;
            }
            let mut incompat_chain = false;

            let anchorjr = seeds_ref[primary_ref_anchors[j].0 as usize].primary_base.unwrap();
            let anchorir = seeds_ref[primary_ref_anchors[i].0 as usize].primary_base.unwrap();
            if anchorjr >= anchorir {
                incompat_chain = true;
            }

            let anchorjq;
            let anchoriq;
            if forward_strand {
                anchorjq = seeds_q[primary_ref_anchors[j].1 as usize].actual_ref_positions[0];
                anchoriq = seeds_q[primary_ref_anchors[i].1 as usize].actual_ref_positions[0];
            } else {
                anchorjq =
                    read_length - seeds_q[primary_ref_anchors[j].1 as usize].actual_ref_positions[0];
                anchoriq =
                    read_length - seeds_q[primary_ref_anchors[i].1 as usize].actual_ref_positions[0];
            }

            //            let anchorjr_n = &seeds_ref[super_anchors[j].0 as usize];
            //            let anchorir_n = &seeds_ref[super_anchors[i].0 as usize];

            if anchorjq >= anchoriq {
                incompat_chain = true;
            }

            let f_cand_i;
            max_num_iter += 1;

            if incompat_chain {
                f_cand_i = f64::MIN;
            } else {
                num_iter += 1;
                //                let ref_order_dist = min_dist_ref_nodes(
                //                    &col_list[&super_anchors[j].0],
                //                    &col_list[&super_anchors[i].0],
                //                );
                let ref_order_dist = (anchorir - anchorjr) as f64;
                let query_order_dist = (anchoriq - anchorjq) as f64;
                let heur_score = heuristic_score_coord(ref_order_dist, query_order_dist);

                f_cand_i = f[j] + heur_score;
            }
            if f_cand_i > best_f_i {
                best_f_i = f_cand_i;
                best_j = j;
            }
        }
        if best_f_i <= 0.0 {
            best_j = i;
        }
        f.push(best_f_i);
        if best_j != usize::MAX {
            pointer_array[i] = best_j;
            interval_pointer_array[i] = interval_pointer_array[best_j];
        }
    }
}

#[inline]
pub fn heuristic_score_coord(ref_order_dist: f64, query_order_dist: f64) -> f64 {
    //    let dist_ref = -2.0 * f64::sqrt(ref_order_dist);
    //    let max_dist = f64::max(ref_order_dist,query_order_dist);
    //    let dist_ref = -1.0 * ref_order_dist;
    //    let dist_query = 30.0 - f64::sqrt(max_dist);
    let mut gap_cost = (ref_order_dist - query_order_dist).abs().powf(1.25);
    if query_order_dist > ref_order_dist {
    } else {
        //        gap_cost = (ref_order_dist - query_order_dist) * 0.2;
    }
    //    let fuzzy_gap_cost = 500.;
    //    let fuzzy_gap_cost = gap_cost.powi(2) / query_order_dist;
    //    let linear_cost = f64::sqrt(ref_order_dist) + query_order_dist;
    //    let score = constants::BANDED_CHAINING_BASE_SCORE - f64::min(gap_cost, fuzzy_gap_cost);
    let score = constants::COORD_CHAIN_BASE_SCORE - gap_cost;
    //    let score = constants::BANDED_CHAINING_BASE_SCORE - linear_cost;
    return score;
}


