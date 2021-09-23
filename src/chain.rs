use debruijn::kmer::Kmer16;
use disjoint_sets::UnionFind;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use crate::data_structs::KmerNode;

pub fn position_max_f64(slice: &[f64]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value0.partial_cmp(value1).unwrap())
        .map(|(idx, _)| idx)
}

fn alpha(j: usize, i: usize, anchors: &Vec<(usize, usize, &KmerNode ,&KmerNode)>, w: usize) -> f64 {
    let dist_multiplier = 5.0;
//    let query_diff = anchors[i].0 as f64 - anchors[j].0 as f64;
//    let ref_diff = anchors[i].1 as f64 - anchors[j].1 as f64;
//    let min_diff = f64::min(ref_diff, query_diff);
//    return f64::min(min_diff, (w) as f64);
    let num_bases_query = f64::max(16.0 - dist_multiplier*(anchors[i].0 as f64 - anchors[j].0 as f64), 0.0);
    let num_bases_ref = f64::max(16.0 - dist_multiplier*(anchors[i].1 as f64 - anchors[j].1 as f64), 0.0);
    return f64::min(num_bases_query,num_bases_ref);
}

fn beta(
    j: usize,
    i: usize,
    anchors: &Vec<(usize, usize, &KmerNode,&KmerNode)>,
    w: usize,
    g: f64,
    best_f_i: f64,
) -> f64 {
    let dist_multiplier = 5.0;
    if anchors[j].0 > anchors[i].0 {
        return f64::MAX;
    }

    let query_diff = anchors[i].0 as f64 - anchors[j].0 as f64;
    let ref_diff = anchors[i].1 as f64 - anchors[j].1 as f64;

//    if f64::max(ref_diff, query_diff) > g {
//        return f64::MAX;
//    }

    let diff = (query_diff - ref_diff) * dist_multiplier;
    if diff == 0.0 {
        return 0.0;
    } else {
//        return (0.01 * 16.0 as f64 * diff) + 0.5 * diff.log2();
        return (diff)*0.01 + 0.5 * diff.log2();
    }
}

pub fn chain_seeds(
    seeds1: &Vec<KmerNode>,
    seeds2: &Vec<KmerNode>,
    pos1: &Vec<usize>,
    pos2: &Vec<usize>,
    w: usize,
    h: usize,
    k: usize,
) {
    let g = f64::MAX;

    let mut mini_hash_map1 = FxHashMap::default();
    let mut positions_to_index1 = FxHashMap::default();
    for (i, kmer_node) in seeds1.iter().enumerate() {
        let kmer = &kmer_node.kmer;
        let pos_vec = mini_hash_map1.entry(kmer).or_insert(vec![]);
        positions_to_index1.entry(pos1[i]).or_insert(kmer_node);
        pos_vec.push(kmer_node.id);
    }

    let mut mini_hash_map2 = FxHashMap::default();
    let mut positions_to_index2 = FxHashMap::default();
    for (i, kmer_node) in seeds2.iter().enumerate() {
        let kmer = &kmer_node.kmer;
        let pos_vec = mini_hash_map2.entry(kmer).or_insert(vec![]);
        positions_to_index2.entry(pos2[i]).or_insert(kmer_node);
        pos_vec.push(kmer_node.id);
    }

    let mut anchors = vec![];

    for kmer in mini_hash_map1.keys() {
        if mini_hash_map2.contains_key(kmer) {
            let query_positions = mini_hash_map1.get(kmer).unwrap();
            let ref_positions = mini_hash_map2.get(kmer).unwrap();
            for p1 in query_positions {
                for p2 in ref_positions {
                    anchors.push((*p1, *p2, &seeds1[*p1], &seeds2[*p2]));
                }
            }
        }
    }

    anchors.sort_by(|a, b| a.0.cmp(&b.0));
    //dbg!(anchors[1],anchors[2],anchors[3]);
    //dbg!(alpha(1,2,&anchors,k),alpha(2,3,&anchors,k));
    //dbg!(beta(1,2,&anchors,k,g),beta(2,3,&anchors,k,g));
    //chaining
    let mut initial_set = FxHashSet::default();
    initial_set.insert(0);

    let mut f = vec![w as f64];
    let mut uf = UnionFind::new(anchors.len());
    let mut pointer_array = vec![];
    for i in 0..anchors.len() {
        pointer_array.push(i);
    }

    for i in 1..anchors.len() {
        let mut best_f_i = w as f64;
        let mut best_j = usize::MAX;
        let mut start = 0;
        if i > h {
            start = i - h
        }
        let mut already_checked = FxHashSet::default();
        for j in start..i {
            let j_rep = uf.find(j);
            if already_checked.contains(&j_rep) {
                //continue;
            }
//            let f_cand_i = f[j] + alpha(j, i, &anchors, k) - beta(j, i, &anchors, k, g, best_f_i);
            let f_cand_i = f[j] + alpha(j, i, &anchors, k) - beta(j, i, &anchors, k, g, best_f_i);
            if f_cand_i > best_f_i {
                best_f_i = f_cand_i;
                best_j = j;
            }
            already_checked.insert(j_rep);
        }
        f.push(best_f_i);
        if best_j != usize::MAX {
            uf.union(i, best_j);
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

    for i in 0..chain_sequence.len(){
       dbg!(anchors[chain_sequence[i]], pos1[anchors[chain_sequence[i]].0]);
    }
    dbg!(f[best_i]);

//    dbg!(
//        (anchors[chain_sequence[0]].1) - (anchors[chain_sequence[chain_sequence.len() - 1]].1),
//        anchors[chain_sequence[0]].1
//    );
//    dbg!(
//        anchors[chain_sequence[0]].0,
//        anchors[chain_sequence[chain_sequence.len() - 1]].0
//    );
}
