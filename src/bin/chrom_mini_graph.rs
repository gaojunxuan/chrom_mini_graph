use bincode;
use bio::io::{fasta,fastq};
use chrom_mini_graph::chain;
use chrom_mini_graph::data_structs::KmerNode;
use chrom_mini_graph::graph_utils;
use chrom_mini_graph::seeding_methods_bit;
use clap::{App, AppSettings, Arg, SubCommand};
use debruijn::dna_string::*;
use debruijn::kmer::Kmer16;
use debruijn::Kmer;
use fxhash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::time::Instant;

fn main() {
    let matches = App::new("chrom_mini_graph")
        .setting(AppSettings::ArgRequiredElseHelp)
        .version("0.1")
        .about("Chromatic minimizer pangenome graph program.")
        .subcommand(
            SubCommand::with_name("generate")
                .about("Generate graph.")
                .version("0.1")
                .arg(
                    Arg::with_name("references")
                        .index(1)
                        .help("Input reference fasta files. Any sequence within every reference is assumed to be homologous with one another.")
                        .takes_value(true)
                        .required(true)
                        .multiple(true),
                ).
                arg(
                    Arg::with_name("output")
                        .short("o")
                        .help("Name of output chromatic reference graph. (Default: serialized_mini_graph)")
                        .takes_value(true),
                ).
                arg(
                    Arg::with_name("syncmer")
                        .short("s")
                        .help("Use syncmers. (Default: use minimizers)")
                )
        )
        .subcommand(
            SubCommand::with_name("map")
                .about("Map sequences onto graph.")
                .version("0.1")
                .arg(
                    Arg::with_name("reference_graph")
                        .required(true)
                        .index(1)
                        .help("Reference graph output from the generate subcommand."),
                ).
                arg(
                    Arg::with_name("reads")
                        .required(true)
                        .index(2)
                        .help("Reads to map to graph"),
                ).
                arg(
                    Arg::with_name("syncmer")
                        .short("s")
                        .help("Use syncmers. (Default: use minimizers)")
                )
        )
        .get_matches();

    let generate;
    let matches_subc;
    if let Some(matches) = matches.subcommand_matches("generate") {
        generate = true;
        matches_subc = matches;
    } else {
        matches_subc = matches.subcommand_matches("map").unwrap();
        generate = false;
    }

    //Input real genomes to create pangenom
    //
    //don't use chain_heuristic; it is the minimap2 heuristic and works poorly
    //on long contigs.
    let chain_heuristic = false;
    let w = 16;
    let k = 16;
    let mask_repet_on_generate = false;
    let s = 10;
    let t = (k - s + 2) / 2 as usize;
    let h = 10;
    //use syncmers if not using minimizers

    let use_minimizers;
    if matches_subc.is_present("syncmer"){
        use_minimizers = false;
    }
    else{
        use_minimizers = true;
    }

    if generate {
        let ref_genomes: Vec<&str> = matches_subc.values_of("references").unwrap().collect();
        let mut chroms = vec![];

        for i in 0..ref_genomes.len() {
            let reader = fasta::Reader::from_file(&ref_genomes[i]);
            for record in reader.unwrap().records() {
                let rec = record.unwrap();
                let chrom = DnaString::from_acgt_bytes(rec.seq());
                chroms.push(chrom);
                println!("{}-th reference is {}.", i, ref_genomes[i]);
            }
        }

        let mut seeds1;
        let seed_p1;
        if use_minimizers {
            seed_p1 = seeding_methods_bit::minimizer_seeds(&chroms[0], w, k);
        } else {
            seed_p1 = seeding_methods_bit::open_sync_seeds(&chroms[0], k, t);
        }
        seeds1 = seed_p1.0;
        let p1 = seed_p1.1;

        println!(
            "Starting reference is {} and has {} nodes.",
            ref_genomes[0],
            seeds1.len()
        );

        for i in 1..chroms.len() {
            println!("-----------------Iteration {}-------------------", i);
            let old_graph_len = seeds1.len();
            let seeds2;
            let now = Instant::now();
            let s2;
            if use_minimizers {
                s2 = seeding_methods_bit::minimizer_seeds(&chroms[i], w, k);
            } else {
                s2 = seeding_methods_bit::open_sync_seeds(&chroms[i], k, t);
            }
            seeds2 = s2.0;
            println!(
                "Generating sketch (minimizers) time: {}",
                now.elapsed().as_secs_f32()
            );
            let now = Instant::now();
            let ref_hash_map = chain::get_kmer_dict_mut(&mut seeds1);
            let q_hash_map = chain::get_kmer_dict(&seeds2);

            ///Test
            let mut kmer_count_dict = FxHashMap::default();
            for node in seeds1.iter() {
                let count_num = kmer_count_dict.entry(node.kmer).or_insert(0);
                *count_num += 1;
            }

            let mut hash_vec: Vec<(&Kmer16, &usize)> = kmer_count_dict.iter().collect();
            hash_vec.sort_by(|a, b| b.1.cmp(a.1));
            let mut dont_use_kmers = FxHashSet::default();

            for i in 0..hash_vec.len() / 10000 as usize {
                dont_use_kmers.insert(hash_vec[i].0);
            }
            ///
            if !mask_repet_on_generate {
                dont_use_kmers = FxHashSet::default();
            }

            let best_anchors = chain::chain_seeds(
                &mut seeds1,
                &seeds2,
                &ref_hash_map,
                &q_hash_map,
                h,
                chain_heuristic,
                false,
                &dont_use_kmers,
            );
            println!("Chaining time: {}", now.elapsed().as_secs_f32());

            let now = Instant::now();
            chain::add_align_to_graph(&mut seeds1, seeds2, best_anchors);
            println!(
                "Generating graph from alignment time: {}",
                now.elapsed().as_secs_f32()
            );
            let now = Instant::now();

            println!(
                "New graph now has {} nodes. Difference is {}.",
                seeds1.len(),
                seeds1.len() - old_graph_len
            );

            chain::top_sort(&mut seeds1);
            println!("Top sort time: {}.", now.elapsed().as_secs_f32());
        }

        let concat_graph = graph_utils::concat_graph(&seeds1[0], &seeds1);
        let mut file = File::create("simplified_mini_graph.csv").unwrap();

        for (n1, n2, weight) in concat_graph.0.iter() {
            let towrite = format!("{},{},{}\n", n1, n2, weight);
            write!(&mut file, "{}", towrite).unwrap();
        }

        let mut file = File::create("simplified_metadata.csv").unwrap();

        for (node_id, vertices) in concat_graph.1.iter() {
            if *node_id as usize > seeds1.len() {
                dbg!(node_id, vertices);
            }
            let node_order = seeds1[*node_id as usize].order;
            let node_color = format!("{:#08b}", seeds1[*node_id as usize].color);
            let mut kmer_list = vec![];
            for vertex in vertices {
                let n = &seeds1[*vertex as usize];
                kmer_list.push(n.kmer.to_string());
            }
            let towrite = format!("{},{},{:?}\n", node_order, node_color, kmer_list);
            write!(&mut file, "{}", towrite).unwrap();
        }

        let mut file_mini_pos = File::create("ref_mini_pos.txt").unwrap();
        for (i, position) in p1.iter().enumerate() {
            write!(&mut file_mini_pos, "{},{}\n", i, position).unwrap();
        }

        let j = serde_json::to_string(&seeds1);

        // Print, write to a file, or send to an HTTP server.
        let serial_name = matches_subc
            .value_of("output")
            .unwrap_or("serialized_mini_graph");
        let serial_json_name = format!("{}.json", serial_name);
        let serial_bin_name = format!("{}.bin", serial_name);

        let mut file_json = File::create(serial_json_name).unwrap();
        write!(&mut file_json, "{}", j.unwrap()).unwrap();

        let now = Instant::now();
        let mut file_bin = BufWriter::new(File::create(serial_bin_name).unwrap());
        bincode::serialize_into(&mut file_bin, &seeds1).unwrap();
        println!(
            "Serializing and writing time {}.",
            now.elapsed().as_secs_f32()
        );

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
    } else {
        let ref_graph_file = matches_subc.value_of("reference_graph").unwrap();
        let ref_graph_f = File::open(ref_graph_file).unwrap();
        let ref_graph_reader = BufReader::new(ref_graph_f);

        let mut ref_graph: Vec<KmerNode> = bincode::deserialize_from(ref_graph_reader).unwrap();

        let order_to_id = chain::top_sort(&mut ref_graph);
        let mut kmer_count_dict = FxHashMap::default();

        for node in ref_graph.iter() {
            let count_num = kmer_count_dict.entry(node.kmer).or_insert(0);
            *count_num += 1;
        }

        let mut hash_vec: Vec<(&Kmer16, &usize)> = kmer_count_dict.iter().collect();
        hash_vec.sort_by(|a, b| b.1.cmp(a.1));
        let mut dont_use_kmers = FxHashSet::default();

        for i in 0..hash_vec.len() / 10000 as usize {
            dont_use_kmers.insert(hash_vec[i].0);
        }
        //        dont_use_kmers.insert(hash_vec[0].0);

        let reads_file = matches_subc.value_of("reads").unwrap();

        let reader = fastq::Reader::from_file(reads_file);
        let mut reads = vec![];
        for record in reader.unwrap().records() {
            let rec = record.unwrap();
            let read = DnaString::from_acgt_bytes(rec.seq());
            let id = rec.id().to_string();
            reads.push((read, id));
        }

        let ref_hash_map = chain::get_kmer_dict(&ref_graph);
        let mut anchor_file = BufWriter::new(File::create("read_anchor_hits.txt").unwrap());

        for (read, id) in reads {
            let read_seeds;
            //            let now = Instant::now();
            let s2;
            if use_minimizers {
                s2 = seeding_methods_bit::minimizer_seeds(&read, w, k);
            } else {
                s2 = seeding_methods_bit::open_sync_seeds(&read, k, t);
            }

            //            println!(
            //                "Generating sketch (minimizers) time: {}",
            //                now.elapsed().as_secs_f32()
            //            );

            read_seeds = s2.0;
            let q_hash_map = chain::get_kmer_dict(&read_seeds);
            let now = Instant::now();
            let best_anchors = chain::chain_seeds(
                &mut ref_graph,
                &read_seeds,
                &ref_hash_map,
                &q_hash_map,
                h,
                chain_heuristic,
                true,
                &dont_use_kmers,
                //&FxHashSet::default(),
            );

            write!(&mut anchor_file, "{}:", id).unwrap();
            for anchor in best_anchors.iter() {
                if anchor == best_anchors.last().unwrap() {
                    write!(&mut anchor_file, "{}\n", anchor.0).unwrap();
                } else {
                    write!(&mut anchor_file, "{},", anchor.0).unwrap();
                }
            }
            println!("Read: {}", id);
            println!("Chaining time: {}", now.elapsed().as_secs_f32());
            let now = Instant::now();
            chain::get_best_path_from_chain(best_anchors, &ref_graph, &order_to_id);
            println!("Path collection time: {}", now.elapsed().as_secs_f32());
        }
    }
}
