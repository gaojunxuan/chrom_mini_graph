use crate::align;
use crate::data_structs::{Anchors, Color, KmerNode};
use fxhash::{FxHashMap, FxHashSet};
use linfa::prelude::*;
use linfa_elasticnet::{ElasticNet, Result};
use linfa_elasticnet::{ElasticNetError, ElasticNetParams};
use ndarray::{array, s, Array, Array2};
use probability::prelude::*;
use rand::distributions::{Bernoulli, Distribution};
use std::mem;
use std::time::Instant;

//pub struct RefGraphReadHits {
//    pub subsampled_graph_from_index: FxHashMap<u32, f64>,
//}
//
//impl RefGraphReadHits {
//    pub fn new() -> RefGraphReadHits {
//        return RefGraphReadHits {
//            subsampled_graph_from_index: FxHashMap::default(),
//        };
//    }
//}

pub fn subsampled_ref_graph(
    ref_graph: &Vec<KmerNode>,
    samp_freq: usize,
    num_genomes: usize,
) -> FxHashSet<u32> {
    let mut return_set = FxHashSet::default();
    for node in ref_graph.iter() {
        if node.id as usize % samp_freq == 0 {
            let num_colors = align::get_nonzero_bits(node.color);
            if num_colors.len() > num_genomes {
                //            if num_colors.len() > num_genomes * 9 / 10 {
                continue;
            }
            let d = Bernoulli::new(num_colors.len() as f64 / num_genomes as f64).unwrap();
            let v = d.sample(&mut rand::thread_rng());
            if v == true {
                return_set.insert(node.id);
            }
        }
    }
    return return_set;
}

pub fn collect_good_chains<'a>(
    successful_chains: &mut Vec<(Anchors, String)>,
    read_chains: &mut Vec<(Anchors, f64, bool)>,
    short_reads: bool,
    read_id: String,
) {
    let cutoff;
    if short_reads {
        cutoff = 0.98;
    } else {
        cutoff = 0.9;
    }
    if read_chains.len() == 0 {
        return;
    }
    let best_chain_score = read_chains
        .iter()
        .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
        .unwrap()
        .1;
    let mut bad_chains = FxHashSet::default();
    for (i, (_chain, score, _strand)) in read_chains.iter().enumerate() {
        if *score < cutoff * best_chain_score {
            println!("Best chain {}, bad chain {}", best_chain_score, score);
            bad_chains.insert(i);
        }
    }
    //No secondary chains TODO test this out
    if read_chains.len() - bad_chains.len() > 1 {
        return;
    }

    for (i, (chain, _score, _strand)) in read_chains.iter_mut().enumerate() {
        if bad_chains.contains(&i) || chain.len() < 5 {
            continue;
        }

        successful_chains.push((mem::take(chain), read_id.clone()));
    }
}

fn get_chain_colour(
    chain: &Anchors,
    ref_graph: &Vec<KmerNode>,
    num_genomes: usize,
    read: &String,
) -> Color {
    let now = Instant::now();
    let mut consensus_color: Color = 0;
    let mut color_vec = vec![0; num_genomes];
    for anchor in chain.iter() {
        for color in align::get_nonzero_bits(ref_graph[anchor.0 as usize].color) {
            color_vec[color] += 1;
        }
    }
    //    println!(
    //        "get_chain_colour time {}",
    //        now.elapsed().as_secs_f32()
    //    );
    let mut best_colors: Vec<(usize, usize)> = color_vec.into_iter().enumerate().collect();
    best_colors.sort_by(|x, y| y.1.cmp(&x.1));
    let mut colors_to_ret = vec![];
    let best_score = best_colors[0].1;
    for (color, score) in best_colors {
        if score < best_score * 97 / 100 {
            break;
        }
        colors_to_ret.push(color);
    }

    //    println!("{:?}, {}", colors_to_ret, read);

    for color in colors_to_ret {
        consensus_color += Color::pow(2, color as u32);
    }
    return consensus_color;
}

fn mcmc(colour_chains: &Vec<Color>, putative_strains: &Vec<usize>) {
    let mut str_in_read_vec = vec![];
    let strain_id_to_ind: FxHashMap<usize, usize> = putative_strains
        .iter()
        .enumerate()
        .map(|x| (*x.1, x.0))
        .collect();
    let mut s_mat = vec![vec![0.; putative_strains.len()]; putative_strains.len()];

    for color in colour_chains {
        let mut str_in_read = vec![];
        for id in putative_strains {
            let colpow = Color::pow(2, *id as u32);
            if color & colpow == colpow {
                str_in_read.push(strain_id_to_ind[id]);
            }
        }
        str_in_read_vec.push(str_in_read);
    }

    for str_vec in str_in_read_vec.iter() {
        let norm;
        if str_vec.len() < 2 {
            //Double degree for edge
            norm = 0.5;
        } else {
            norm = (str_vec.len() * (str_vec.len() - 1) / 2) as f64;
        }
        for id1 in str_vec.iter() {
            for id2 in str_vec.iter() {
                if id1 != id2 || str_vec.len() == 1 {
                    s_mat[*id1][*id2] += 1. / norm;
                }
            }
        }
    }

    println!("{:?}", s_mat);
    println!("{:?}", strain_id_to_ind);
    let pi: Vec<f64> = s_mat.iter().map(|x| x.iter().sum()).collect();
    let pi_sum: f64 = pi.iter().sum();
    let pi: Vec<f64> = pi.iter().map(|x| x / (pi_sum)).collect();
    for i in 0..s_mat.len() {
        s_mat[i][i] = 0.;
    }
    for vec in s_mat.iter_mut() {
        let vec_sum: f64 = vec.iter().sum();
        for elt in vec.iter_mut() {
            *elt = *elt / vec_sum;
        }
    }
    let k12 = 1200.;
    let ec = 2000.;
    let w = 400.;
    let h5 = 800.;
    let fift = 0.;
    let twosix = 0.;
    let total_reads = k12 + ec + w + h5 + fift + twosix;
    let true_vec = vec![
        (50, k12 / total_reads),
        (58, w / total_reads),
        (59, h5 / total_reads),
        (60, ec / total_reads),
//        (15, fift / total_reads),
//        (26, twosix / total_reads),
    ];
    let mut true_pi = vec![0.; pi.len()];
    for (ind, abund) in true_vec {
        true_pi[strain_id_to_ind[&ind]] = abund;
    }

    //Calculate likelihood
    let epsilon = 0.010;

    let mut source = source::default();
    let iterations = 1000;
    let mut best_dist = true_pi.clone();
    let mut dist = pi.clone();
    let mut prev_ll = 0.;
    let mut best_ll = f64::MIN;
    for iter in 0..iterations {
        let mut ll = 0.;
        let mut cand_dist = dist.clone();
        if iter != 0{
            for j in 0..dist.len(){
                let distribution = Gaussian::new(dist[j], f64::min(dist[j]/10., 0.1));
                let sampler = Independent(&distribution, &mut source);
                let sample = sampler.take(1).sum::<f64>();
                if sample < 0. || sample > 1.{
                    continue;
                }
                else{
                    cand_dist[j] = sample;
                }
            }
        }
        let cand_sum : f64 = cand_dist.iter().sum();
        let cand_dist : Vec<f64> = cand_dist.iter().map(|x| x/cand_sum).collect();
        for str_vec in str_in_read_vec.iter() {
            let mut prod_terms = vec![];
            let mut inter = 0.;
            for i in 0..putative_strains.len() {
                let gi = cand_dist[i] + epsilon;
                if str_vec.len() == 0 {
                    break;
                }
                let num_edges = str_vec.len() * (str_vec.len() - 1) / 2;
                let mut prod_term;

                if str_vec.contains(&i) {
                    prod_term = 1. - epsilon;
                    for j in 0..putative_strains.len() {
                        if j == i {
                            continue;
                        }
                        if str_vec.contains(&j) {
//                            prod_term *= (s_mat[i][j] + epsilon) * gi;
                        }
                    }
                } else {
//                    prod_term = f64::powi(epsilon, num_edges as i32);
                    prod_term = epsilon;
                }
                inter += prod_term * gi;
                prod_terms.push(prod_term);
            }
            if str_vec.len() == 0 {
                continue;
            }
//            ll += f64::ln(inter / prod_terms.iter().sum::<f64>());
            ll += f64::ln(inter);
        }
        dbg!(ll);
        if iter == 0 {
            prev_ll = ll;
        }
        else{
            if best_ll < ll{
                best_ll = ll;
                best_dist = cand_dist.clone();
            }
            let ratio = f64::exp(ll - prev_ll);
            let distribution = Uniform::new(0.,1.);
            let sampler = Independent(&distribution, &mut source);
            let sample = sampler.take(1).sum::<f64>();
            if sample < ratio{
                prev_ll = ll;
                dist = cand_dist;
            }
        }
    }
    dbg!(dist, true_pi, best_dist, best_ll);
}

pub fn solve_coeffs(
    successful_chains: &Vec<(Anchors, String)>,
    sampled_nodes: FxHashSet<u32>,
    ref_graph: &Vec<KmerNode>,
    num_genomes: usize,
    chroms: &Vec<String>,
    penalty: f64,
) {
    //Massage data structures to get it into matrix form
    let mut hit_counter = vec![0.; num_genomes];
    let mut nodes_with_hits = FxHashMap::default();
    let mut colour_node_hit_vec = vec![FxHashSet::default(); num_genomes];
    let mut num_without_cons = 0;
    let use_slack = true;
    let dag_chain = true;
    let mut success_chain_colours = vec![];
    let mut putative_strains = vec![];

    let now = Instant::now();
    //Look through chains to make color-coherent node-hit-matrix
    for (chain, read_id) in successful_chains {
        let slack_colours;
        if use_slack {
            slack_colours = get_chain_colour(chain, &ref_graph, num_genomes, read_id);
        } else {
            slack_colours = Color::MAX;
        }
        let mut consensus_color = Color::MAX;
        if !dag_chain {
            for anchor in chain.iter() {
                consensus_color &= ref_graph[anchor.0 as usize].color;
            }
            assert!(consensus_color != 0);
        } else {
            //Can do some rescue here probably TODO
            consensus_color = slack_colours
                & ref_graph[chain[0].0 as usize].color
                & ref_graph[chain[chain.len() - 1].0 as usize].color;
            if consensus_color == 0 {
                num_without_cons += 1;
                continue;
            }
            //            println!("{:?}", align::get_nonzero_bits(consensus_color));
        }
        //Need to pick a representative colour to traverse the path.
        let rep_color = Color::pow(2, align::get_first_nonzero_bit(consensus_color) as u32);
        let mut curr_node = &ref_graph[chain[0].0 as usize];

        if use_slack{
            success_chain_colours.push(slack_colours);
        }
        else{
            success_chain_colours.push(consensus_color);
        }
        while curr_node.id != chain[chain.len() - 1].0 {
            let primary_color;
            if use_slack {
                primary_color = slack_colours;
            } else {
                primary_color = consensus_color;
            }
            if primary_color & curr_node.color == primary_color {
                if sampled_nodes.contains(&curr_node.id) {
                    let num_hits = nodes_with_hits.entry(curr_node.id).or_insert(0.);
                    *num_hits += 1.;
                    for color in align::get_nonzero_bits(primary_color & curr_node.color) {
                        colour_node_hit_vec[color].insert(curr_node.id);
                    }
                }
            }
            let mut found = false;
            for (_distance, (edge_color, node_ind)) in curr_node.child_edge_distance.iter() {
                if *edge_color & rep_color == rep_color {
                    curr_node = &ref_graph[curr_node.child_nodes[*node_ind as usize] as usize];
                    found = true;
                }
            }
            if !found {
                //                println!("Circularity issue detected. Fix later");
                //                dbg!(curr_node, align::get_nonzero_bits(consensus_color), align::get_nonzero_bits(curr_node.child_edge_distance[0].1.0));
                break;
            }
        }
    }
    println!(
        "Processed colour-hit matrix {}",
        now.elapsed().as_secs_f32()
    );

    let now = Instant::now();
    //Greedy set cover
    let covering_colours = colour_cover(&nodes_with_hits, &colour_node_hit_vec);
    //    let covering_colours : Vec<usize> = (0..num_genomes).collect();
    dbg!(&covering_colours, covering_colours.len());
    dbg!(&num_without_cons);

    //Massage data structures into matrix form.
    let num_samples = nodes_with_hits.len();
    let mut values_mat = Array2::<f64>::zeros((num_samples, covering_colours.len()));
    let mut target_vec = vec![0.; num_samples];
    let mut id_to_pos = FxHashMap::default();
    let mut id_to_pos_vec: Vec<(u32, f64)> = nodes_with_hits.iter().map(|x| (*x.0, *x.1)).collect();
    id_to_pos_vec.sort_by(|x, y| x.0.cmp(&y.0));
    for (i, (id, _num_hits)) in id_to_pos_vec.iter().enumerate() {
        id_to_pos.insert(*id, i);
    }

    for (i, color) in covering_colours.iter().enumerate() {
        let colour_hits = &colour_node_hit_vec[*color];
        for id in colour_hits {
            values_mat[[id_to_pos[&id], i]] = 1.;
        }
    }
    for id in id_to_pos.keys() {
        target_vec[id_to_pos[&id]] = nodes_with_hits[&id];
    }

    for i in 0..covering_colours.len() {
        hit_counter[i] = values_mat.slice(s![.., i]).sum();
    }

    dbg!(&hit_counter);
    dbg!(penalty);
    let ds = Dataset::new(values_mat, Array::from_vec(target_vec));
    let model1 = ElasticNetParams::new()
        .penalty(penalty)
        .l1_ratio(1.0)
        .with_intercept(false);

    //    let mut models = vec![];
    //    for i in 0..50{
    //         let model = ElasticNetParams::new()
    //        .penalty(i as f64 / 250.)
    //        .l1_ratio(1.0)
    //        .with_intercept(false);
    //        models.push(model);
    //    }
    //    let r2_scores = ds.cross_validate(5,&models, |prediction, truth| prediction.r2(truth)).unwrap();
    //    dbg!(r2_scores);

    let result = model1.fit(&ds).unwrap();
    dbg!(result.hyperplane());

    for (i, value) in result.hyperplane().iter().enumerate() {
        if value > &0.0 {
            println!(
                "Chrom {} Colour {} Score {} Num hits {}",
                chroms[chroms.len() - covering_colours[i] - 1],
                covering_colours[i],
                value,
                hit_counter[i]
            );
            putative_strains.push(covering_colours[i]);
        }
    }
    println!("Solving LASSO {}", now.elapsed().as_secs_f32());
    putative_strains.sort();

    mcmc(&success_chain_colours, &putative_strains);

    //dbg!(&values_vec[3*num_samples..4*num_samples]);
}

pub fn colour_cover(
    hit_counts: &FxHashMap<u32, f64>,
    colour_hit_vec: &Vec<FxHashSet<u32>>,
) -> Vec<usize> {
    let mut uncovered: FxHashSet<u32> = hit_counts.iter().map(|x| *x.0).collect();
    let mut return_cover = FxHashSet::default();
    while !uncovered.is_empty() {
        let mut best_i = usize::MAX;
        let mut best_intersect = 0;
        for (i, set) in colour_hit_vec.iter().enumerate() {
            let mut int_count = 0;
            for elt in set.iter() {
                if uncovered.contains(elt) {
                    int_count += 1
                }
            }
            if int_count > best_intersect {
                best_i = i;
                best_intersect = int_count;
            }
        }
        if best_intersect < 10 {
            break;
        }
        for elt in colour_hit_vec[best_i].iter() {
            uncovered.remove(elt);
        }
        println!("Greediest color {} with {} hits", best_i, best_intersect);
        return_cover.insert(best_i);
    }

    let mut toret: Vec<usize> = return_cover.into_iter().collect();
    toret.sort();
    return toret;
}
