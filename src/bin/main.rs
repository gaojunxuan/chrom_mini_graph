use debruijn::dna_string::*;
use mini_alph::chain;
use mini_alph::graph_utils;
use mini_alph::seeding_methods_bit;
use mini_alph::simulation_utils_bit;
use serde_json::{Result, Value};
use std::env;
use std::fs::File;
use std::io::prelude::*;
use std::time::Instant;

fn main() {
    let args: Vec<String> = env::args().collect();

    //Input real genomes to create pangenom
    let using_genomes;
    if args.len() > 1 {
        using_genomes = true;
    } else {
        using_genomes = false;
    }

    let num_iters = 1;
    let chain_heuristic = false;
    let now = Instant::now();
    let print_file = true;
    let theta = 0.001;
    let w = 16;
    let t = 3;
    let k = 16;
    let h = 10;
    let use_minimizers = true;

    if !using_genomes {
        let mut seeds1;
        let n = 6 * usize::pow(10, 6);
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
            let best_anchors = chain::chain_seeds(&mut seeds1, &seeds2, h, chain_heuristic);
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

            for (n1, n2, weight) in concat_graph.iter() {
                let towrite = format!("{},{}\n", n1, n2);
                write!(&mut file, "{}", towrite).unwrap();
            }

            let mut file = File::create("graph.txt").unwrap();

            for node in seeds1.iter() {
                for child in node.child_nodes.iter() {
                    let towrite = format!(
                        "{}-{},{}-{}\n",
                        node.id,
                        node.color,
                        seeds1[*child as usize].id,
                        seeds1[*child as usize].color
                    );
                    write!(&mut file, "{}", towrite).unwrap();
                }
            }
        }
    } else {
        use bio::io::fasta;

        let mut chroms = vec![];

        for i in 1..args.len() {
            let reader = fasta::Reader::from_file(&args[i]);
            for record in reader.unwrap().records() {
                let rec = record.unwrap();
                println!("{}", rec.id());

                let chrom = DnaString::from_acgt_bytes(rec.seq());
                chroms.push(chrom);
            }
        }

        let mut seeds1;
        let (s1, _p1) = seeding_methods_bit::minimizer_seeds(&chroms[0], w, k);
        dbg!(s1.len());
        seeds1 = s1;

        for i in 1..chroms.len() {
            let mut old_graph_len = seeds1.len();
            let seeds2;
            let now = Instant::now();
            let (s2, _p2) = seeding_methods_bit::minimizer_seeds(&chroms[i], w, k);
            seeds2 = s2;
            println!("Generating kmers {}", now.elapsed().as_secs_f32());
            dbg!(seeds2.len());
            let now = Instant::now();
            let best_anchors = chain::chain_seeds(&mut seeds1, &seeds2, h, chain_heuristic);
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
            let mut file = File::create("simplified_mini_graph.csv").unwrap();

            for (n1, n2, weight) in concat_graph.iter() {
                let towrite = format!("{},{},{}\n", n1, n2, weight);
                write!(&mut file, "{}", towrite).unwrap();
            }
        }

        let j = serde_json::to_string(&seeds1);

        // Print, write to a file, or send to an HTTP server.
        let mut file_json = File::create("serialized_mini_graph.json").unwrap();
        write!(&mut file_json, "{}", j.unwrap()).unwrap();

        let mut file = File::create("full_mini_graph.csv").unwrap();
        for node in seeds1.iter() {
            for child in node.child_nodes.iter() {
                let towrite = format!(
                    "{}-{},{}-{}\n",
                    node.order, node.id, seeds1[*child as usize].order, seeds1[*child as usize].id,
                );
                write!(&mut file, "{}", towrite).unwrap();
            }
        }
    }
}
