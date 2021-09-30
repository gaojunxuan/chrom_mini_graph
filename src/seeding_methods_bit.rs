use debruijn::dna_string::*;
use smallvec::SmallVec;
use debruijn::kmer::Kmer12;
use debruijn::kmer::Kmer16;
use debruijn::Vmer;
use fxhash::hash;
use crate::data_structs::KmerNode;

fn position_min<T: Ord>(slice: &[T]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value1.cmp(value0))
        .map(|(idx, _)| idx)
}

pub fn minimizer_seeds(s: &DnaString, w: usize, k: usize) -> (Vec<KmerNode>, Vec<u32>) {
    let mut minimizer_seeds:Vec<KmerNode> = vec![];
    let mut positions_selected: Vec<u32> = Vec::new();

    //look at windows
    let mut running_pos = 0;
    let mut min_running_pos = usize::MAX;
    let mut window_hashes = vec![0; w];
    for i in 0..s.len() - k + 1{
        let kmer: Kmer16 = s.slice(i, i + k).get_kmer(0);
        window_hashes[running_pos] = hash(&kmer);
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
        let kmer_node = KmerNode{
            kmer: kmer,
            id: positions_selected.len() as u32,
            order: positions_selected.len() as u32,
            color: 1,
            child_nodes: SmallVec::<[u32; 1]>::new(),
        };
        minimizer_seeds.push(kmer_node);
        //        let pos_vec = minimizer_seeds
        //            .entry(kmer)
        //            .or_insert(FxHashSet::default());
        //        pos_vec.insert(i-offset);
        positions_selected.push((i - offset) as u32);

        running_pos += 1;
        running_pos %= w;
    }

    for i in 0..minimizer_seeds.len()-1{
        minimizer_seeds[i].child_nodes.push((i+1) as u32);
    }

    return (minimizer_seeds, positions_selected);
}

pub fn open_sync_seeds(string: &DnaString, k: usize, t: usize) -> (Vec<KmerNode>, Vec<u32>) {
    let s = 12;
    let mut syncmer_seeds = vec![];
    let mut positions_selected: Vec<u32> = Vec::new();
    //hash all s-mers

    let w = k - s + 1;
    let mut running_pos = 0;
    let mut min_running_pos = usize::MAX;
    let mut window_hashes = vec![0; w];

    for i in 0..string.len() - s + 1{
        let smer: Kmer12 = string.slice(i, i + s).get_kmer(0);
        window_hashes[running_pos] = hash(&smer);
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
            if running_pos - min_running_pos == t-1 {
                let kmer: Kmer16 = string.slice(i-w+1, i-w+1+k).get_kmer(0);
                let kmer_node = KmerNode{
                    kmer: kmer,
                    id: positions_selected.len() as u32,
                    order: positions_selected.len() as u32,
                    color: 1,
                    child_nodes: SmallVec::<[u32; 1]>::new(),
                    };

                positions_selected.push(i as u32);
                syncmer_seeds.push(kmer_node);
            }
        } else {
            if w - (min_running_pos - running_pos) == t-1 {
                let kmer: Kmer16 = string.slice(i+1-w, i+1+k-w).get_kmer(0);
                let kmer_node = KmerNode{
                    kmer: kmer,
                    id: positions_selected.len() as u32,
                    order: positions_selected.len() as u32,
                    color: 1,
                    child_nodes: SmallVec::<[u32; 1]>::new(),
                    };
                positions_selected.push(i as u32);
                syncmer_seeds.push(kmer_node);
            }
        }

        running_pos += 1;
        running_pos %= w;
    }

    for i in 0..syncmer_seeds.len()-1{
        syncmer_seeds[i].child_nodes.push((i+1) as u32);
    }

    return (syncmer_seeds, positions_selected);
}
