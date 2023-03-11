use crate::align;
use crate::data_structs::{Anchors, Color};
use crate::data_structs::{BamInfo, KmerNode};
use block_aligner::cigar::*;
use block_aligner::scan_block::*;
use block_aligner::scores::*;
use debruijn::dna_string::DnaString;
use debruijn::dna_string::*;
use debruijn::Kmer;
use debruijn::Mer;
use fxhash::FxHashSet;
use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::record::{Cigar as hts_Cigar, CigarString, Record};
use rust_htslib::bam::{Format, HeaderView, Writer};
use simple_logger::SimpleLogger;
use std::time::Instant;

#[inline]
pub fn get_nonzero_bits_fast(n: Color) -> Vec<usize> {
    let mut bit_poses = vec![];
    let mut temp = n;
    bit_poses.reserve(128);
    assert!(n != 0);
    loop {
        if temp == 0 {
            break;
        }
        let index = temp.trailing_zeros();
        bit_poses.push(index as usize);
        temp ^= 1 << index;
    }

    bit_poses
}
#[inline]
pub fn get_nonzero_bits(n: Color) -> Vec<usize> {
    let mut temp = n;
    let mut bit_poses = vec![];
    let mut bit_pos = 0;

    assert!(n != 0);
    let mut count = 0;
    loop {
        if count == 500 {
            break;
        }
        if temp % 2 != 0 {
            bit_poses.push(bit_pos);
        }
        temp /= 2;
        bit_pos += 1;
        count += 1;
    }

    bit_poses
}

#[inline]
pub fn get_first_nonzero_bit(n: Color) -> usize {
    return n.trailing_zeros() as usize;
}

pub fn get_coords(
    anchors: &Vec<(u32, u32)>,
    ref_nodes: &Vec<KmerNode>,
    q_nodes: &Vec<KmerNode>,
    color: Color,
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
    if color == Color::MAX {
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
    let chrom = &chroms[chroms.len() - bit_pos - 1].0;
    let strand = chroms[chroms.len() - bit_pos - 1].1;
    let mut num_trav = 0;
    let mut current_anchor = 0;
    let mut visited_nodes = FxHashSet::default();
    loop {
        if current_anchor < anchors.len() {
            if parent_node.id == anchors[current_anchor].0 {
                kmer_hit_positions.push((
                    path_dist as i64,
                    q_nodes[anchors[current_anchor].1 as usize].actual_ref_positions[0],
                ));
                current_anchor += 1;
            }
        }
        if num_trav > 50000 {
            //            dbg!(&parent_node);
            //            dbg!(&ref_nodes[parent_node.child_nodes[0] as usize]);
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
        //Debugging
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
                        .slice(abs_pos_index + path_dist, abs_pos_index + 456 + path_dist)
                        .to_string(),
                    rep_kmer.to_string(),
                    path_dist + abs_pos_index
                );
                dbg!(&parent_node.child_edge_distance);
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
        if parent_node.order >= last_node.order && abs_pos_index != 0 {
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
            //|| path_dist > 10000000
            || visited_nodes.contains(&parent_node.id)
        //|| num_trav > 80000
        {
            log::warn!("Path localization failed for read; TODO fix this.");
            //            dbg!(num_trav, found, bit, bit_pos);
            //            dbg!(&parent_node, &last_node, &first_node);
            //            dbg!(&ref_nodes[parent_node.child_nodes[0] as usize]);
            //            dbg!(&parent_node.order, last_node.order);
            //            dbg!(visited_nodes.contains(&parent_node.id));
            //            dbg!(path_dist);
            //            dbg!(num_trav);
            break;
        }
    }

    //    dbg!(abs_pos_index);
    if abs_pos_index == 0 {
        return ((i64::MAX, i64::MAX), vec![]);
    }
    assert!(abs_pos_index != 0);
    let total_interval;
    if strand {
        total_interval = (
            abs_pos_index as i64,
            abs_pos_index + path_dist + kmer_length,
        );
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
        total_interval = (
            abs_pos_index as i64 - path_dist as i64 - kmer_length as i64,
            abs_pos_index,
        );
        for pos in kmer_hit_positions.iter_mut() {
            pos.0 = abs_pos_index as i64 - pos.0;
        }
    }
    return (
        (total_interval.0 as i64, total_interval.1 as i64),
        kmer_hit_positions,
    );
}

pub fn write_bam_header(
    chroms: &Vec<(DnaString, bool)>,
    ref_names: &Vec<String>,
    bam_name: String,
) -> (HeaderView, Writer) {
    let mut header = Header::new();
    for (i, (reference, _strand)) in chroms.iter().enumerate() {
        let mut new_rec = HeaderRecord::new(b"SQ");
        new_rec.push_tag(b"SN", &ref_names[i]);
        new_rec.push_tag(b"LN", &reference.len());
        header.push_record(&new_rec);
    }
    let writer = Writer::from_path(bam_name, &header, Format::Bam).unwrap();
    return (HeaderView::from_header(&header), writer);
}

pub fn get_bam_record(bam_info: BamInfo, headerview: &HeaderView) -> Record {
    let mut rec = Record::new();
    let cigar_vec = &bam_info.cigar;
    let mut hts_cigar_vec = vec![];
    for (i, oplen) in cigar_vec.iter().enumerate() {
        let op = oplen.op;
        let len = oplen.len;
        let hts_op;

        if op == Operation::M {
            hts_op = hts_Cigar::Match(len as u32);
        } else if op == Operation::D {
            if i == 0 || i == cigar_vec.len() - 1 {
                //                hts_op = hts_Cigar::HardClip(len as u32);
                hts_op = hts_Cigar::Del(len as u32);
            } else {
                hts_op = hts_Cigar::Del(len as u32);
            }
        } else if op == Operation::I {
            if i == 0 || i == cigar_vec.len() - 1 {
                //                hts_op = hts_Cigar::HardClip(len as u32);
                hts_op = hts_Cigar::Ins(len as u32);
            } else {
                hts_op = hts_Cigar::Ins(len as u32);
            }
        } else {
            hts_op = hts_Cigar::RefSkip(len as u32);
        }

        hts_cigar_vec.push(hts_op);
    }
    let mut scaled_quals = vec![0; bam_info.quals.len()];
    for (i, val) in bam_info.quals.iter().enumerate() {
        scaled_quals[i] = val - 33;
    }
    let hts_cigar_view = CigarString(hts_cigar_vec);
    let hts_cigar = Some(hts_cigar_view);
    rec.set(
        bam_info.qname.as_bytes(),
        hts_cigar.as_ref(),
        &bam_info.sequence.as_bytes(),
        &scaled_quals,
    );

    if bam_info.strand == false {
        rec.set_reverse();
    }

    let tid = headerview.tid(bam_info.ref_name.as_bytes()).unwrap();
    rec.set_tid(tid as i32);
    rec.set_pos(bam_info.map_pos);
    rec.set_mapq(bam_info.mapq);
    return rec;
}

pub fn align_from_chain(
    anchors: &Anchors,
    chroms: &Vec<(DnaString, bool)>,
    color: Color,
    ref_graph: &Vec<KmerNode>,
    read_seeds: &Vec<KmerNode>,
    read: &DnaString,
    read_strand: bool,
    quals: &[u8],
    chrom_names: &Vec<String>,
    read_id: &String,
    //    headerview: &HeaderView,
    //    writer: &mut Writer,
) -> Option<BamInfo> {
    log::trace!(
        "Aligning to genome corresponding to colour {} (or {})",
        color.trailing_zeros(),
        chrom_names[(chroms.len() as u32 - color.trailing_zeros() - 1) as usize]
    );
    if anchors.len() < 3 {
        log::trace!("Less than 3 anchors, bad align");
        log::trace!("Alignment score: NA");
        return None;
    }

    let ref_chrom = &chroms[chroms.len() - align::get_first_nonzero_bit(color) - 1].0;
    let strand_chrom = &chroms[chroms.len() - align::get_first_nonzero_bit(color) - 1].1;
    let now = Instant::now();
    let (_ref_coords, kmer_hit_coords) =
        align::get_coords(&anchors, &ref_graph, &read_seeds, color, &chroms);
    if kmer_hit_coords.len() < 3 {
        log::trace!("Less than 3 kmer hits, bad align");
        log::trace!("Alignment score: NA");
        return None;
    }
    //            dbg!(kmer_hit_coords[0]);
    //            dbg!(kmer_hit_coords[1]);
    //            dbg!(kmer_hit_coords.last().unwrap());
    for (i, (c1, c2)) in kmer_hit_coords.iter().enumerate() {
        let adj_c1;
        if *c1 >= ref_chrom.len() as i64 {
            //We don't want to do this right now because I believe minimizer distance
            //across the ends of a reference genome isn't done properly...
            adj_c1 = c1 - ref_chrom.len() as i64;
        } else if *c1 < 0 {
            adj_c1 = ref_chrom.len() as i64 + *c1;
        } else {
            adj_c1 = *c1;
        }
        let node = &read_seeds[anchors[i].1 as usize];
        let r_node = &ref_graph[anchors[i].0 as usize];
        let d1 = read.slice(*c2, *c2 + 16).to_string();
        let d2 = node.kmer.to_string();
        let d3 = node.kmer.rc().to_string();
        let d4 = ref_chrom
            .slice(adj_c1 as usize, adj_c1 as usize + 16)
            .to_string();
        let d5 = ref_chrom
            .slice(adj_c1 as usize, adj_c1 as usize + 16)
            .rc()
            .to_string();
        let d6 = read.slice(*c2, *c2 + 16).rc().to_string();
        if d1 != d4 && d1 != d5 {
            println!(
                "NOT FOUND {} ! read and ref {} {} {} {}, len {}, c1 {}",
                i,
                d1,
                d6,
                d4,
                d5,
                kmer_hit_coords.len(),
                c1
            );
        }
        if d1 != d2 && d1 != d3 {
            println!("NOT FOUND {} ! read and node {} {} {}", i, d1, d2, d3);
        }
    }

    //GET THE REFERENCE STRING
    let ref_map_string;
    let a;
    let b;
    if *strand_chrom {
        a = kmer_hit_coords[0].0;
        if kmer_hit_coords.last().unwrap().0 < a {
            b = ref_chrom.len() as i64;
            println!("Circular mapping b. Cut off TODO. Not implemented yet?");
        } else {
            b = kmer_hit_coords.last().unwrap().0 + 16;
        }
        if b as usize > ref_chrom.len() {
            println!("End of chromosome mapping issue. Continue");
            return None;
        }
    } else {
        b = kmer_hit_coords[0].0 + 16;
        if kmer_hit_coords.last().unwrap().0 < 0 {
            a = 0;
            println!("Circular mapping a. Cut off TODO. Not implemented yet?");
        } else {
            a = kmer_hit_coords.last().unwrap().0;
        }
        if b as usize > ref_chrom.len() {
            println!("End of chromosome mapping issue. Continue");
            return None;
        }
    }

    let read_map_string;
    let read_map_string_slice;
    let qual_map_string: Vec<u8>;
    let c;
    let d;
    if read_strand {
        c = kmer_hit_coords[0].1;
        d = kmer_hit_coords.last().unwrap().1 + 16;
    } else {
        c = kmer_hit_coords.last().unwrap().1;
        d = kmer_hit_coords[0].1 + 16;
        //Need to reverse the quals TODO
    }
    let mut num_mm = 0;
    let mut break_mm = 3;
    let mut i_left = 1;
    let mut i_right = 1;
    if read_strand && *strand_chrom {
        while a as usize >= i_left && c >= i_left {
            if ref_chrom.get(a as usize - i_left) == read.get(c - i_left) {
            } else {
                num_mm += 1;
                if num_mm > break_mm {
                    break;
                }
            }
            i_left += 1;
        }
        while b as usize + i_right <= ref_chrom.len() && d + i_right <= read.len() {
            if ref_chrom.get(b as usize + i_right - 1) == read.get(d + i_right - 1) {
            } else {
                num_mm += 1;
                if num_mm > break_mm{
                    break;
                }
            }
            i_right += 1;
        }
    } else if read_strand && !*strand_chrom {
        while b as usize + i_left <= ref_chrom.len() && c >= i_left {
            if ref_chrom.get(b as usize + i_left - 1) == reverse_comp_04(read.get(c - i_left)) {
            } else {
                num_mm += 1;
                if num_mm > break_mm{
                    break;
                }
            }
            i_left += 1;
        }
        while a as usize >= i_right && d + i_right <= read.len() {
            if ref_chrom.get(a as usize - i_right) == reverse_comp_04(read.get(d - 1 + i_right)) {
            } else {
                num_mm += 1;
                if num_mm >  break_mm{
                    break;
                }
            }
                i_right += 1;
        }
    } else if !read_strand && *strand_chrom {
        while a as usize >= i_left && d + i_left <= read.len() {
            if ref_chrom.get(a as usize - i_left) == reverse_comp_04(read.get(d - 1 + i_left)) {
            } else {
                num_mm += 1;
                if num_mm > break_mm{
                    break;
                }
            }
                i_left += 1;
        }
        while b as usize + i_right <= ref_chrom.len() && c >= i_right {
            if ref_chrom.get(b as usize + i_right - 1) == reverse_comp_04(read.get(c - i_right)) {
            } else {
                num_mm += 1;
                if num_mm > break_mm {
                    break;
                }
            }
                i_right += 1;
        }
    } else {
        while b as usize + i_left <= ref_chrom.len() && d + i_left <= read.len() {
            if ref_chrom.get(b as usize + i_left - 1) == read.get(d + i_left - 1) {
            } else {
                num_mm += 1;
                if num_mm > break_mm{
                    break;
                }
            }
                i_left += 1;
        }
        while a as usize >= i_right && c >= i_right {
            if ref_chrom.get(a as usize - i_right) == read.get(c - i_right) {
            } else {
                num_mm += 1;
                if num_mm >  break_mm{
                    break;
                }
            }
                i_right += 1;
        }
    }
    i_left -= 1;
    i_right -= 1;
    log::trace!("chrom, read, {}, {}", *strand_chrom, read_strand);
    log::trace!(
        "i_left {}, i_right {}, a {}, b {}, c {}, d{}",
        i_left,
        i_right,
        a,
        b,
        c,
        d
    );
    if *strand_chrom {
        ref_map_string = ref_chrom
            .slice(a as usize - i_left, b as usize + i_right)
            .to_string();
    } else {
        ref_map_string = ref_chrom
            .slice(a as usize - i_right, b as usize + i_left)
            .rc()
            .to_string();
    }
    if read_strand {
        read_map_string = read.slice(c - i_left, d + i_right).to_string();
        read_map_string_slice = read.slice(c - i_left, d + i_right);
        qual_map_string = quals[c - i_left..d + i_right].to_vec();
    } else {
        read_map_string = read.slice(c - i_right, d + i_left).rc().to_string();
        read_map_string_slice = read.slice(c - i_right, d + i_left).rc();
        qual_map_string = quals[c - i_right..d + i_left]
            .iter()
            .rev()
            .map(|x| *x)
            .collect::<Vec<u8>>()
            .to_owned();
    }

    let start_pos_chrom;
    if read_strand == *strand_chrom {
        start_pos_chrom = a;
    } else {
        start_pos_chrom = a;
    }
    let block_align = true;

    //ALIGNMENT
    // if block_align {
        let now = Instant::now();
        let block_size = 512;
        //        let block_size = 64;
        let r_cpy = ref_map_string.clone();
        let q_cpy = read_map_string.clone();

        //        dbg!(&ref_map_string, &read_map_string.clone());
        let r = PaddedBytes::from_string::<NucMatrix>(ref_map_string, block_size);
        let q = PaddedBytes::from_string::<NucMatrix>(read_map_string.clone(), block_size);
        let gaps = Gaps {
            open: -2,
            extend: -1,
        };

        let nuc_mat: NucMatrix = NucMatrix::new_simple(1, -2);
        let a = Block::<_, true, false>::align(&q, &r, &nuc_mat, gaps, block_size..=block_size, 50);

        let res = a.res();
        let cigar = a.trace().cigar(res.query_idx, res.reference_idx);
        //        let fmt_string = cigar.format(&q_cpy.as_bytes(), &r_cpy.as_bytes());
        //        println!("{}\n{}", fmt_string.0, fmt_string.1);
        //        println!("Aligning time: {}", now.elapsed().as_secs_f32());
        log::trace!("Alignment score: {}", res.score);
        let ref_chrom_name = &chrom_names[chroms.len() - align::get_first_nonzero_bit(color) - 1];
        let seq;
        let cigar_vec: Vec<OpLen>;
        if *strand_chrom {
            seq = read_map_string;
            cigar_vec = cigar.to_vec();
        } else {
            seq = read_map_string_slice.rc().to_string();
            cigar_vec = cigar.to_vec().into_iter().rev().collect();
        }
        let mapq_score = res.score as f64 / (q.len() + (read.len() - q.len()) * 2) as f64;
        let mapq;
        if mapq_score < 0. {
            mapq = 0;
        } else {
            let mapq_temp = mapq_score.ln() * 40. + 65.;
            if mapq_temp < 0. {
                mapq = 0;
            } else {
                mapq = mapq_temp as i32;
            }
        }
        let mapq = u8::min(mapq as u8, 60);

        let bam_info = BamInfo {
            cigar: cigar_vec,
            sequence: seq,
            quals: qual_map_string,
            qname: read_id.clone(),
            strand: read_strand,
            ref_name: ref_chrom_name.clone(),
            map_pos: start_pos_chrom,
            mapq: mapq,
        };

        log::trace!("Read align time: {}", now.elapsed().as_secs_f32());
        return Some(bam_info);
    //        let bam_rec = align::get_bam_record(
    //            bam_info,
    //            &headerview,
    //        );
    //
    //        writer.write(&bam_rec).unwrap();
    // } 
    // else {
    //     //Testing out WFA alignment. Seems to be worse than blockaligner.
    //     use libwfa::{affine_wavefront::*, bindings::*, mm_allocator::*, penalties::*};

    //     let alloc = MMAllocator::new(BUFFER_SIZE_1G as u64);
    //     let min_wavefront_length = 10;
    //     let max_distance_threshold = 50;

    //     let pattern = ref_map_string.to_string();
    //     let text = read_map_string.to_string();

    //     let mut penalties = AffinePenalties {
    //         match_: 0,
    //         mismatch: 4,
    //         gap_opening: 6,
    //         gap_extension: 2,
    //     };

    //     let pat_len = pattern.as_bytes().len();
    //     let text_len = text.as_bytes().len();

    //     let mut wavefronts = AffineWavefronts::new_reduced(
    //         pat_len,
    //         text_len,
    //         &mut penalties,
    //         min_wavefront_length,
    //         max_distance_threshold,
    //         &alloc,
    //     );

    //     wavefronts
    //         .align(pattern.as_bytes(), text.as_bytes())
    //         .unwrap();

    //     let score = wavefronts.edit_cigar_score(&mut penalties);

    //     println!("score: {}", score);

    //     // The cigar can also be extracted as a byte vector
    //     let cigar = wavefronts.cigar_bytes_raw();
    //     let cg_str = std::str::from_utf8(&cigar).unwrap();

    //     let mut op_vec = vec![];
    //     let mut prev_char = 'A';
    //     let mut running_len = 0;
    //     for c in cg_str.chars() {
    //         if c != prev_char && running_len != 0 {
    //             let op;
    //             if prev_char == 'X' {
    //                 op = Operation::M;
    //             } else if prev_char == 'M' {
    //                 op = Operation::M;
    //             } else if prev_char == 'I' {
    //                 op = Operation::I;
    //             } else if prev_char == 'D' {
    //                 op = Operation::D;
    //             } else {
    //                 panic!("{}\n{}", prev_char, cg_str);
    //                 op = Operation::Sentinel;
    //             }

    //             op_vec.push(OpLen {
    //                 op: op,
    //                 len: running_len,
    //             });
    //             running_len = 0;
    //             prev_char = c;
    //         }
    //         running_len += 1;
    //         prev_char = c;
    //     }
    //     let op;
    //     if prev_char == 'X' {
    //         op = Operation::M;
    //     } else if prev_char == 'M' {
    //         op = Operation::M;
    //     } else if prev_char == 'I' {
    //         op = Operation::I;
    //     } else if prev_char == 'D' {
    //         op = Operation::D;
    //     } else {
    //         panic!("{}\n{}", prev_char, cg_str);
    //     }
    //     op_vec.push(OpLen {
    //         op: op,
    //         len: running_len,
    //     });

    //     let ref_chrom_name = &chrom_names[chroms.len() - align::get_first_nonzero_bit(color) - 1];
    //     println!("cigar_str: {}", cg_str);
    //     let bam_info = BamInfo {
    //         cigar: op_vec,
    //         sequence: read_map_string,
    //         quals: qual_map_string,
    //         qname: read_id.clone(),
    //         strand: read_strand,
    //         ref_name: ref_chrom_name.clone(),
    //         map_pos: start_pos_chrom,
    //         mapq: 60,
    //     };

    //     println!("Read align time: {}", now.elapsed().as_secs_f32());
    //     return Some(bam_info);
    // }
}

fn reverse_comp_04(a: u8) -> u8 {
    if a == 0 {
        return 3;
    } else if a == 1 {
        return 2;
    } else if a == 2 {
        return 1;
    } else if a == 3 {
        return 0;
    } else {
        panic!("Only 0-4 can be input in this function");
    }
}
