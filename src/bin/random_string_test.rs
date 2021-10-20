use bincode;
use std::io::BufWriter;
use fxhash::FxHashSet;
use chrom_mini_graph::chain;
use chrom_mini_graph::graph_utils;
use chrom_mini_graph::seeding_methods_bit;
use chrom_mini_graph::simulation_utils_bit;
use std::fs::File;
use std::io::prelude::*;
use std::time::Instant;

fn main() {
    let num_iters = 1;
    //don't use chain_heuristic; it is the minimap2 heuristic and works poorly
    //on long contigs.
    let chain_heuristic = false;
    let now = Instant::now();
    let print_file = true;
    let theta = 0.001;
    let w = 16;
    let k = 16;
    let t = 3;
    let h = 10;
    //use syncmers if not using minimizers
    let use_minimizers = true;

    let mut seeds1;
    let n = 6 * usize::pow(10, 7);
    let s = simulation_utils_bit::gen_rand_string(n);

    println!("Generating random strings {}", now.elapsed().as_secs_f32());
    //    println!("{}",s.to_string());
    //    println!("{:?}",printable_p);

    //    let seeds3;
    //    let positions_selected_p2;

    let now = Instant::now();

    if use_minimizers {
        let (s1, _p1) = seeding_methods_bit::minimizer_seeds(&s, w, k);
        seeds1 = s1;
    } else {
        let (s1, _p1) = seeding_methods_bit::open_sync_seeds(&s, k, t);
        seeds1 = s1;
    }
    println!("Generating kmers {}", now.elapsed().as_secs_f32());
    println!("Starting graph size {}", seeds1.len());

    let mut old_graph_len = seeds1.len();

    for _i in 0..num_iters {
        let (s_p, mut _mutated_pos) = simulation_utils_bit::gen_mutated_string(&s, theta);
        let seeds2;
        if use_minimizers {
            let (s2, _p2) = seeding_methods_bit::minimizer_seeds(&s_p, w, k);
            seeds2 = s2;
        } else {
            let (s2, _p2) = seeding_methods_bit::open_sync_seeds(&s_p, k, t);
            seeds2 = s2;
        }

        let now = Instant::now();
        let ref_hash_map = chain::get_kmer_dict(&seeds1);
        let q_hash_map = chain::get_kmer_dict(&seeds2);
        let best_anchors = chain::chain_seeds(&mut seeds1, &seeds2, &ref_hash_map, &q_hash_map,h, chain_heuristic, false, &FxHashSet::default());
        println!("Chaining {}", now.elapsed().as_secs_f32());

        let now = Instant::now();
        chain::add_align_to_graph(&mut seeds1, seeds2, best_anchors);
        println!("Generating new graph {}", now.elapsed().as_secs_f32());
        let now = Instant::now();

        println!(
            "New graph now has {} nodes. Difference is {}",
            seeds1.len(),
            seeds1.len() - old_graph_len
        );
        old_graph_len = seeds1.len();

        chain::top_sort(&mut seeds1);
        println!("Top sort {}", now.elapsed().as_secs_f32());
    }
    if print_file {
        let concat_graph = graph_utils::concat_graph(&seeds1[0], &seeds1);
        let mut file = File::create("graph_concat.txt").unwrap();

        for (n1, n2, _weight) in concat_graph.0.iter() {
            let towrite = format!("{},{}\n", n1, n2);
            write!(&mut file, "{}", towrite).unwrap();
        }

        let mut file = File::create("graph.txt").unwrap();

        for node in seeds1.iter() {
            for child in node.child_nodes.iter() {
                let towrite = format!(
                    "{}-{},{}-{}\n",
                    node.id, node.color, seeds1[*child as usize].id, seeds1[*child as usize].color
                );
                write!(&mut file, "{}", towrite).unwrap();
            }
        }

        let encoded_graph: Vec<u8> = bincode::serialize(&seeds1).unwrap();
        let now = Instant::now();
        let mut file_bin = BufWriter::new(File::create("test.bin").unwrap());
        bincode::serialize_into(&mut file_bin, &encoded_graph).unwrap();
        println!(
            "Serializing and writing time {}.",
            now.elapsed().as_secs_f32()
        );
    }
}
