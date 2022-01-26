use crate::data_structs::{KmerNode,Color};
use debruijn::dna_string::*;
use debruijn::kmer::Kmer10;
use debruijn::kmer::Kmer12;
use debruijn::kmer::Kmer16;
use debruijn::kmer::Kmer8;
use debruijn::Mer;
use debruijn::Vmer;
use fnv::FnvHasher;
use fxhash::hash;
use fxhash::{FxHashSet, FxHashMap};
use smallvec::SmallVec;
use std::hash::{Hash as _Hash, Hasher as _Hasher};

pub fn get_masked_kmers(s: &DnaString, w: usize, k: usize, fraction_mask_f64 : f64) -> FxHashSet<Kmer16>{
    //Get the discarded k-mers here and don't use these k-mers when seeding
    let (seeds1,_p1) = minimizer_seeds(s, w, k, 100, &FxHashSet::default());
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
) -> (Vec<KmerNode>, Vec<u32>) {
    let use_fnv = true;
    let mut minimizer_seeds: Vec<KmerNode> = vec![];
    let mut positions_selected: Vec<u32> = Vec::new();

    //look at windows
    let mut running_pos = 0;
    let mut min_running_pos = usize::MAX;
    let mut window_hashes = vec![0; w];
    if s.len() < k + 1 {
        return (vec![], vec![]);
    }
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

        if !dont_use_kmers.contains(&node_kmer) {
            positions_selected.push((i - offset) as u32);
            let mut kmer_node = KmerNode {
                kmer: node_kmer,
                id: positions_selected.len() as u32 - 1,
                order: positions_selected.len() as u32 - 1,
                color: 1,
                child_nodes: SmallVec::<[u32; 1]>::new(),
                child_edge_distance: SmallVec::<[(u16, (Color, u8)); 1]>::new(),
                //            child_nodes: vec![],
                canonical: canonical,
                actual_ref_positions: SmallVec::<[usize; 0]>::new(),
            };
            if positions_selected.len() % samp_freq == 0 {
                kmer_node
                    .actual_ref_positions
                    .push(*positions_selected.last().unwrap() as usize);
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

    return (minimizer_seeds, positions_selected);
}

pub fn open_sync_seeds(
    string: &DnaString,
    k: usize,
    t: usize,
    samp_freq: usize,
    dont_use_kmers: &FxHashSet<Kmer16>,
) -> (Vec<KmerNode>, Vec<u32>) {
    let s = 10;
    let mut syncmer_seeds = vec![];
    let mut positions_selected: Vec<u32> = Vec::new();
    //hash all s-mers

    let w = k - s + 1;
    let mut running_pos = 0;
    let mut min_running_pos = usize::MAX;
    let mut window_hashes = vec![0; w];

    for i in 0..string.len() - s + 1 {
        let smer: Kmer10 = string.slice(i, i + s).get_kmer(0);
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

                if !dont_use_kmers.contains(&node_kmer) {
                    positions_selected.push((i + 1 - w) as u32);
                    let mut kmer_node = KmerNode {
                        kmer: node_kmer,
                        id: positions_selected.len() as u32 - 1,
                        order: positions_selected.len() as u32 - 1,
                        color: 1,
                        child_nodes: SmallVec::<[u32; 1]>::new(),
                        child_edge_distance: SmallVec::<[(u16, (Color, u8)); 1]>::new(),
                        canonical: canonical, //                    child_nodes: vec![],
                        actual_ref_positions: SmallVec::<[usize; 0]>::new(),
                    };

                    if positions_selected.len() % samp_freq == 0 {
                        kmer_node
                            .actual_ref_positions
                            .push(*positions_selected.last().unwrap() as usize);
                    }

                    syncmer_seeds.push(kmer_node);
                }
            }
        } else {
            if w - (min_running_pos - running_pos) == t - 1 {
                let kmer: Kmer16 = string.slice(i + 1 - w, i + 1 + k - w).get_kmer(0);
                let canonical;
                let mut node_kmer = kmer.rc();
                if node_kmer < kmer {
                    canonical = false;
                } else {
                    canonical = true;
                    node_kmer = kmer;
                }

                if !dont_use_kmers.contains(&node_kmer) {
                    positions_selected.push((i + 1 - w) as u32);
                    let mut kmer_node = KmerNode {
                        kmer: node_kmer,
                        id: positions_selected.len() as u32 - 1,
                        order: positions_selected.len() as u32 - 1,
                        color: 1,
                        child_nodes: SmallVec::<[u32; 1]>::new(),
                        child_edge_distance: SmallVec::<[(u16, (Color, u8)); 1]>::new(),
                        canonical: canonical, //                    child_nodes: vec![],
                        actual_ref_positions: SmallVec::<[usize; 0]>::new(),
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
            syncmer_seeds[i].child_nodes.push(0 as u32);
            //TODO this is incorrect
            let dist_on_genome = 1;
            syncmer_seeds[i]
                .child_edge_distance
                .push((dist_on_genome as u16, (1, 0)));
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
