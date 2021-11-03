use crate::data_structs::KmerNode;
use block_aligner::cigar::*;
use block_aligner::scan_block::*;
use block_aligner::scores::*;
use debruijn::dna_string::*;
use debruijn::Kmer;
use debruijn::Mer;
use fxhash::FxHashSet;

pub fn get_first_nonzero_bit(n: usize) -> usize {
    let mut temp = n;
    let mut bit_pos = 0;

    assert!(n != 0);
    loop {
        if temp % 2 == 0 {
            temp /= 2;
            bit_pos += 1
        } else {
            break;
        }
    }

    bit_pos
}

pub fn get_coords(
    anchors: &Vec<(u32, u32)>,
    ref_nodes: &Vec<KmerNode>,
    q_nodes: &Vec<KmerNode>,
    color: u64,
    chroms: &Vec<(DnaString, bool)>,
) -> ((i64, i64), Vec<(i64, usize)>) {
    //    for anchor in anchors.iter() {
    //        println!(
    //            "{},{},{},{}",
    //            ref_nodes[anchor.0 as usize].kmer.to_string(),
    //            q_nodes[anchor.1 as usize].kmer.to_string(),
    //            ref_nodes[anchor.0 as usize].canonical,
    //            q_nodes[anchor.1 as usize].canonical
    //        );
    //    }
    if color == u64::MAX {
        return ((i64::MAX, i64::MAX), vec![]);
    }
    if anchors.len() == 0 {
        return ((i64::MAX, i64::MAX), vec![]);
    }
    let kmer_length = 16;
    let debug = false;
    let mut kmer_hit_positions = vec![];
    let mut temp = color;
    let mut bit = 1;
    let mut bit_pos = 0;

    assert!(color != 0);
    loop {
        if temp % 2 == 0 {
            temp /= 2;
            bit *= 2;
            bit_pos += 1
        } else {
            break;
        }
    }

    let mut path_dist = 0;
    let mut abs_pos_index = 0;
    let first_node = &ref_nodes[anchors[0].0 as usize];
    let last_node = &ref_nodes[anchors[anchors.len() - 1].0 as usize];
    let mut parent_node = first_node;
    let mut chrom = &chroms[chroms.len() - bit_pos - 1].0;
    let strand = chroms[chroms.len() - bit_pos - 1].1;
    let mut num_trav = 0;
    let mut current_anchor = 0;
    let mut visited_nodes = FxHashSet::default();
    loop {
        if parent_node.id == anchors[current_anchor].0 {
            kmer_hit_positions.push((
                path_dist as i64,
                q_nodes[anchors[current_anchor].1 as usize].actual_ref_positions[0],
            ));
            current_anchor += 1;
        }
        if num_trav > 50000 {
            dbg!(&parent_node);
            dbg!(&ref_nodes[parent_node.child_nodes[0] as usize]);
        }
        visited_nodes.insert(parent_node.id);
        num_trav += 1;

        if parent_node.actual_ref_positions.len() > 0 && abs_pos_index == 0 {
            let mut tmp = parent_node.color;
            let mut shift_count = 0;
            let mut offset_count = 0;
            //            dbg!(bit, tmp, bit_pos);
            loop {
                //                dbg!(shift_count, tmp, offset_count);
                if tmp == 0 {
                    panic!();
                }
                if tmp % 2 == 0 {
                    tmp = tmp >> 1;
                    shift_count += 1;
                } else {
                    if shift_count == bit_pos {
                        break;
                    } else {
                        tmp = tmp >> 1;
                        shift_count += 1;
                        offset_count += 1;
                    }
                }
            }

            offset_count = parent_node.actual_ref_positions.len() - offset_count - 1;

            if abs_pos_index == 0 {
                if strand {
                    abs_pos_index = parent_node.actual_ref_positions[offset_count] - path_dist;
                } else {
                    abs_pos_index = parent_node.actual_ref_positions[offset_count] + path_dist;
                }
            }
            //            dbg!(&parent_node.actual_ref_positions, path_dist, offset_count, bit, abs_pos_index);
        }
        if abs_pos_index != 0 && debug {
            let rep_kmer;
            if parent_node.canonical {
                rep_kmer = parent_node.kmer;
            } else {
                rep_kmer = parent_node.kmer.rc();
            }
            if strand {
                println!(
                    "{}, {}, {}",
                    chrom
                        .slice(abs_pos_index + path_dist, abs_pos_index + 16 + path_dist)
                        .to_string(),
                    rep_kmer.to_string(),
                    path_dist
                );
                if chrom
                    .slice(abs_pos_index + path_dist, abs_pos_index + path_dist + 16)
                    .to_string()
                    != rep_kmer.to_string()
                {
                    dbg!(bit, bit_pos, strand);
                    dbg!(abs_pos_index, abs_pos_index + path_dist);
                    dbg!(&parent_node, &last_node, &first_node);
                    dbg!(&ref_nodes[parent_node.child_nodes[0] as usize]);
                    panic!();
                }
            } else {
                dbg!(
                    &parent_node,
                    bit_pos,
                    abs_pos_index,
                    path_dist,
                    abs_pos_index - path_dist
                );
                println!(
                    "{}, {}, {}",
                    chrom
                        .slice(abs_pos_index - path_dist, abs_pos_index - path_dist + 16)
                        .to_string(),
                    rep_kmer.rc().to_string(),
                    path_dist
                );
                if chrom
                    .slice(abs_pos_index - path_dist, abs_pos_index - path_dist + 16)
                    .to_string()
                    != rep_kmer.rc().to_string()
                {
                    panic!();
                }
            }
        }
        if parent_node.id == last_node.id {
            break;
        }
        let mut found = false;
        //        for (a,child_id) in parent_node.child_nodes.iter().enumerate() {
        //            if &ref_nodes[*child_id as usize].color & bit != 0 {
        for dist in parent_node.child_edge_distance.iter() {
            if dist.1 .0 & bit != 0 {
                path_dist += dist.0 as usize;
                //                dbg!(dist, &parent_node);
                found = true;
                parent_node = &ref_nodes[parent_node.child_nodes[dist.1 .1 as usize] as usize];
                break;
            }
            //                }
            //                dbg!(&parent_node);

            //            }
        }
        //        if !found{
        //            if &ref_nodes[parent_node.child_nodes[0] as usize].color & bit != 0 &&
        //            parent_node.child_nodes.len() == 1{
        //                dbg!("missed connection?", &parent_node, bit);
        //                parent_node = &ref_nodes[parent_node.child_nodes[0] as usize];
        //                continue;
        //            }
        //        }
        if !found
            || parent_node.order > last_node.order
            || path_dist > 500000
            || visited_nodes.contains(&parent_node.id)
            || num_trav > 80000
        {
            dbg!(num_trav, found, bit);
            dbg!(&parent_node, &last_node, &first_node);
            dbg!(&ref_nodes[parent_node.child_nodes[0] as usize]);
            panic!();
        }
    }

    //    dbg!(abs_pos_index);
    if abs_pos_index == 0 {
        return ((i64::MAX, i64::MAX), vec![]);
    }
    assert!(abs_pos_index != 0);
    let total_interval;
    if strand {
        total_interval = (abs_pos_index as i64, abs_pos_index + path_dist + kmer_length);
        for pos in kmer_hit_positions.iter_mut() {
            pos.0 += abs_pos_index as i64;
        }
    } else {
//        if abs_pos_index < path_dist + kmer_length {
//            dbg!(bit, bit_pos, strand);
//            dbg!(abs_pos_index, path_dist + kmer_length);
//            dbg!(&parent_node, &last_node, &first_node);
//            dbg!(&ref_nodes[parent_node.child_nodes[0] as usize]);
//
//            panic!();
//        }
        total_interval = (abs_pos_index as i64 - path_dist as i64 - kmer_length as i64, abs_pos_index);
        for pos in kmer_hit_positions.iter_mut() {
            pos.0 = abs_pos_index  as i64 - pos.0;
        }
    }
    return ((total_interval.0 as i64, total_interval.1 as i64), kmer_hit_positions);
}
