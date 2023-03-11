use crate::data_structs::{Color, KmerNode};
use debruijn::dna_string::*;
use debruijn::kmer::Kmer10;
use debruijn::kmer::Kmer12;
use debruijn::kmer::Kmer16;
use debruijn::kmer::Kmer8;
use debruijn::Mer;
use debruijn::Vmer;
use debruijn::*;
use fnv::FnvHasher;
use fxhash::hash;
use fxhash::{FxHashMap, FxHashSet};
use smallvec::SmallVec;
use std::fs::File;
use std::hash::{Hash as _Hash, Hasher as _Hasher};
use std::io::{BufRead, BufReader};

pub fn get_masked_kmers(
    s: &DnaString,
    w: usize,
    k: usize,
    s_sync: usize,
    t: usize,
    fraction_mask_f64: f64,
    use_minimizers: bool,
    frequent_kmers: &FxHashMap<Kmer16, usize>,
) -> FxHashSet<Kmer16> {
    //Get the discarded k-mers here and don't use these k-mers when seeding
    let seeds1;
    if use_minimizers {
        let (seeds, _p1) =
            minimizer_seeds(s, w, k, 100, &FxHashSet::default(), frequent_kmers, false);
        seeds1 = seeds;
    } else {
        let (seeds, _p1) = 
            open_sync_seeds(s, k, t, s_sync, 100, &FxHashSet::default(), frequent_kmers, false);
        seeds1 = seeds;
    }
    let mut kmer_count_dict = FxHashMap::default();
    for node in seeds1.into_iter() {
        let count_num = kmer_count_dict.entry(node.kmer).or_insert(0);
        *count_num += 1;
    }

    let mut hash_vec: Vec<(Kmer16, usize)> = kmer_count_dict.into_iter().collect();
    hash_vec.sort_by(|a, b| b.1.cmp(&a.1));
    let mut dont_use_kmers = FxHashSet::default();

    for i in 0..(hash_vec.len() as f64 * fraction_mask_f64) as usize {
        dont_use_kmers.insert(hash_vec[i].0);
    }
    return dont_use_kmers;
}

fn position_min<T: Ord>(slice: &[T]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value1.cmp(value0))
        .map(|(idx, _)| idx)
}

pub fn minimizer_seeds(
    s: &DnaString,
    w: usize,
    k: usize,
    samp_freq: usize,
    dont_use_kmers: &FxHashSet<Kmer16>,
    frequent_kmers: &FxHashMap<Kmer16, usize>,
    primary_reference: bool,
) -> (Vec<KmerNode>, Vec<u32>) {
    let use_fnv = false;
    let mut minimizer_seeds: Vec<KmerNode> = vec![];
    let mut positions_selected: Vec<u32> = Vec::new();

    //look at windows
    let mut running_pos = 0;
    let mut min_running_pos = usize::MAX;
    let mut window_hashes = vec![0; w];
    if s.len() < k + 1 {
        return (vec![], vec![]);
    }
    let mut num_samp_coord = 0;
    for i in 0..s.len() - k + 1 {
        let kmer: Kmer16 = s.slice(i, i + k).get_kmer(0);
        let rc_kmer = kmer.rc();
        let hash_kmer;
        if kmer < rc_kmer {
            hash_kmer = kmer;
        } else {
            hash_kmer = rc_kmer;
        }
        if use_fnv {
            let mut state = FnvHasher::default();
            hash_kmer.hash(&mut state);
            window_hashes[running_pos] = state.finish() as usize;
        } else {
            window_hashes[running_pos] = hash(&hash_kmer);
        }
        //winnowmap weighting test TODO
        if frequent_kmers.contains_key(&hash_kmer) {
            let mut hash_val = window_hashes[running_pos];
            let num_iters = (frequent_kmers[&hash_kmer] as f64).log(2.0);
            for _i in 0..num_iters as usize {
                hash_val += (usize::MAX - hash_val) / 2;
            }
            window_hashes[running_pos] = hash_val;
        }
        if i < w - 1 {
            continue;
        }

        if min_running_pos == usize::MAX {
            min_running_pos = position_min(&window_hashes).unwrap();
        } else {
            if min_running_pos == running_pos {
                min_running_pos = position_min(&window_hashes).unwrap();
                if window_hashes[min_running_pos] == window_hashes[running_pos] {
                    min_running_pos = running_pos;
                }
            } else {
                if window_hashes[running_pos] < window_hashes[min_running_pos] {
                    min_running_pos = running_pos;
                } else {
                    running_pos += 1;
                    running_pos %= w;
                    continue;
                }
            }
        }

        let offset;
        if min_running_pos > running_pos {
            offset = w - (min_running_pos - running_pos);
        } else {
            offset = running_pos - min_running_pos;
        }

        let kmer: Kmer16 = s.slice(i - offset, i - offset + k).get_kmer(0);
        let canonical;
        let mut node_kmer = kmer.rc();
        if node_kmer < kmer {
            canonical = false;
        } else {
            canonical = true;
            node_kmer = kmer;
        }

        let mut distance_from_last = 0;
        let mut distance_from_start = 0;
        let mut sample_coord = false;
        if !positions_selected.last().is_none() {
            distance_from_last = (i - offset) as u32 - *positions_selected.last().unwrap() as u32;
            distance_from_start = i - offset;
        }

        if distance_from_last as usize > 500 {
            sample_coord = true;
            //            ignore_mask = true;
            //            dbg!(distance_from_last);
        }
        if !dont_use_kmers.contains(&node_kmer) || (sample_coord && primary_reference) {
            if sample_coord {
                num_samp_coord += 1;
            }
            positions_selected.push((i - offset) as u32);
            let mut kmer_node = KmerNode {
                kmer: node_kmer,
                id: positions_selected.len() as u32 - 1,
                order: positions_selected.len() as u32 - 1,
                order_val: distance_from_start as u32,
                color: 1,
                child_nodes: SmallVec::<[u32; 1]>::new(),
                child_edge_distance: SmallVec::<[(u16, (Color, u8)); 1]>::new(),
                //            child_nodes: vec![],
                canonical: canonical,
                actual_ref_positions: SmallVec::<[usize; 0]>::new(),
                repetitive: sample_coord,
                primary_base: Some(*positions_selected.last().unwrap() as u32),
            };
            if positions_selected.len() % samp_freq == 0 || sample_coord {
                kmer_node
                    .actual_ref_positions
                    .push(*positions_selected.last().unwrap() as usize);
            }

            //Repetitive causes far spaced k-mers, make sure to index pairs
            //that are far apart
            if sample_coord {
                minimizer_seeds.last_mut().unwrap().repetitive = true;
                if minimizer_seeds
                    .last()
                    .unwrap()
                    .actual_ref_positions
                    .is_empty()
                {
                    minimizer_seeds
                        .last_mut()
                        .unwrap()
                        .actual_ref_positions
                        .push(positions_selected[positions_selected.len() - 2] as usize);
                }
            }

            minimizer_seeds.push(kmer_node);
        }
        //        let pos_vec = minimizer_seeds
        //            .entry(kmer)
        //            .or_insert(FxHashSet::default());
        //        pos_vec.insert(i-offset);

        running_pos += 1;
        running_pos %= w;
    }

    for i in 0..minimizer_seeds.len() {
        if i == minimizer_seeds.len() - 1 {
            minimizer_seeds[i].child_nodes.push(0 as u32);
            //TODO this is incorrect -- why isthis incorrect??
            let dist_on_genome = positions_selected[0] + s.len() as u32 - positions_selected[i];
            minimizer_seeds[i]
                .child_edge_distance
                .push((dist_on_genome as u16, (1, 0)));
        } else {
            minimizer_seeds[i].child_nodes.push((i + 1) as u32);
            let dist_on_genome = positions_selected[i + 1] - positions_selected[i];
            minimizer_seeds[i]
                .child_edge_distance
                .push((dist_on_genome as u16, (1, 0)));
        }
    }

    //TODO
    //    if positions_selected.len() < 10000{
    //        dbg!(&positions_selected);
    //    }
    //    dbg!(num_samp_coord);
    return (minimizer_seeds, positions_selected);
}

pub fn open_sync_seeds(
    string: &DnaString,
    k: usize,
    t: usize,
    s: usize,
    samp_freq: usize,
    dont_use_kmers: &FxHashSet<Kmer16>,
    frequent_kmers: &FxHashMap<Kmer16, usize>,
    is_primary: bool,
) -> (Vec<KmerNode>, Vec<u32>) {
    let mut syncmer_seeds = vec![];
    let mut positions_selected: Vec<u32> = Vec::new();
    //hash all s-mers

    let w = k - s + 1;
    let mut running_pos = 0;
    let mut min_running_pos = usize::MAX;
    let mut window_hashes = vec![0; w];

    for i in 0..string.len() - s + 1 {
        let smer;
        if s == 8{
            smer = string.slice(i, i + s).to_owned();
        }
        else if s == 10{
            smer = string.slice(i, i + s).to_owned();
        }
        else{
            panic!("s must be set to 8 or 10 right now");
        }
        let rc_smer = smer.rc();
        let hash_smer;
        if smer < rc_smer {
            hash_smer = smer;
        } else {
            hash_smer = rc_smer;
        }
        window_hashes[running_pos] = hash(&hash_smer);
        if i < w - 1 {
            continue;
        }

        if min_running_pos == usize::MAX {
            min_running_pos = position_min(&window_hashes).unwrap();
        } else {
            if min_running_pos == running_pos {
                min_running_pos = position_min(&window_hashes).unwrap();
            } else {
                if window_hashes[running_pos] < window_hashes[min_running_pos] {
                    min_running_pos = running_pos;
                }
            }
        }

        if running_pos > min_running_pos {
            if running_pos - min_running_pos == t - 1 {
                let kmer: Kmer16 = string.slice(i - w + 1, i - w + 1 + k).get_kmer(0);
                let canonical;
                let mut node_kmer = kmer.rc();
                if node_kmer < kmer {
                    canonical = false;
                } else {
                    canonical = true;
                    node_kmer = kmer;
                }

                let mut sample_coord = false;
                let mut distance_from_last = 0;
                let mut distance_from_start = 0;
                if !positions_selected.last().is_none() {
                    distance_from_last =
                        (i + 1 - w) as u32 - *positions_selected.last().unwrap() as u32;
                    distance_from_start = i + 1 - w;
                }
                if distance_from_last > 500 {
                    sample_coord = true;
                    //            dbg!(distance_from_last);
                }

                if !dont_use_kmers.contains(&node_kmer) || (sample_coord && is_primary) {
                    positions_selected.push((i + 1 - w) as u32);
                    let mut kmer_node = KmerNode {
                        kmer: node_kmer,
                        id: positions_selected.len() as u32 - 1,
                        order: positions_selected.len() as u32 - 1,
                        order_val: distance_from_start as u32,
                        color: 1,
                        child_nodes: SmallVec::<[u32; 1]>::new(),
                        child_edge_distance: SmallVec::<[(u16, (Color, u8)); 1]>::new(),
                        canonical: canonical, //                    child_nodes: vec![],
                        actual_ref_positions: SmallVec::<[usize; 0]>::new(),
                        repetitive: sample_coord,
                        primary_base: Some(*positions_selected.last().unwrap() as u32),
                    };

                    if positions_selected.len() % samp_freq == 0 || sample_coord {
                        kmer_node
                            .actual_ref_positions
                            .push(*positions_selected.last().unwrap() as usize);
                    }

                    syncmer_seeds.push(kmer_node);
                }
            }
        } else {
            if w - (min_running_pos - running_pos) == t - 1 {
                let mut distance_from_last = 0;
                let mut distance_from_start = 0;
                if !positions_selected.last().is_none() {
                    distance_from_last =
                        (i + 1 - w) as u32 - *positions_selected.last().unwrap() as u32;
                    distance_from_start = i + 1 - w;
                }
                let kmer: Kmer16 = string.slice(i + 1 - w, i + 1 + k - w).get_kmer(0);
                let canonical;
                let mut node_kmer = kmer.rc();
                if node_kmer < kmer {
                    canonical = false;
                } else {
                    canonical = true;
                    node_kmer = kmer;
                }

                let mut sample_coord = false;
                if distance_from_last > 500 {
                    sample_coord = true;
                    //            dbg!(distance_from_last);
                }

                if !dont_use_kmers.contains(&node_kmer) || (sample_coord && is_primary) {
                    positions_selected.push((i + 1 - w) as u32);
                    let mut kmer_node = KmerNode {
                        kmer: node_kmer,
                        id: positions_selected.len() as u32 - 1,
                        order: positions_selected.len() as u32 - 1,
                        order_val: distance_from_start as u32,
                        color: 1,
                        child_nodes: SmallVec::<[u32; 1]>::new(),
                        child_edge_distance: SmallVec::<[(u16, (Color, u8)); 1]>::new(),
                        canonical: canonical, //                    child_nodes: vec![],
                        actual_ref_positions: SmallVec::<[usize; 0]>::new(),
                        repetitive: false,
                        primary_base: Some(*positions_selected.last().unwrap() as u32),
                    };
                    if positions_selected.len() % samp_freq == 0 {
                        kmer_node
                            .actual_ref_positions
                            .push(*positions_selected.last().unwrap() as usize);
                    }

                    syncmer_seeds.push(kmer_node);
                }
            }
        }

        running_pos += 1;
        running_pos %= w;
    }

    for i in 0..syncmer_seeds.len() {
        if i == syncmer_seeds.len() - 1 {
//            syncmer_seeds[i].child_nodes.push(0 as u32);
//            //TODO this is incorrect
//            let dist_on_genome = 1;
//            syncmer_seeds[i]
//                .child_edge_distance
//                .push((dist_on_genome as u16, (1, 0)));
        } else {
            syncmer_seeds[i].child_nodes.push((i + 1) as u32);
            let dist_on_genome = positions_selected[i + 1] - positions_selected[i];
            syncmer_seeds[i]
                .child_edge_distance
                .push((dist_on_genome as u16, (1, 0)));
        }
    }

    return (syncmer_seeds, positions_selected);
}

pub fn read_minimizer_count_file(file: &str) -> FxHashMap<Kmer16, usize> {
    let f = File::open(file).expect("Unable to open file");
    let f = BufReader::new(f);
    let mut return_set = FxHashMap::default();
    for line in f.lines() {
        let l = line.unwrap();
        let splitted: Vec<&str> = l.split_whitespace().collect();
        let kmer = Kmer16::from_ascii(splitted[0].as_bytes());
        return_set.insert(kmer, splitted[1].parse::<usize>().unwrap());
    }
    return return_set;
}
