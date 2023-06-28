use bincode;
use cmg_shared::data_structs::Bubble;
use clap::ArgAction;
use clap::ArgMatches;
use cmg_shared::data_structs::Cmg;
use simple_logger::SimpleLogger;
use log::LevelFilter;
use bio::io::{fasta, fastq};
use block_aligner::scan_block::*;
use block_aligner::scores::*;
use chrom_mini_graph::align;
use chrom_mini_graph::chain;
use chrom_mini_graph::constants;
use chrom_mini_graph::coord_chain;
use cmg_shared::data_structs::KmerNode;
//use chrom_mini_graph::deconvolution;
use chrom_mini_graph::graph_utils;
use chrom_mini_graph::seeding_methods_bit;
use clap::{Arg, Command};
use debruijn::dna_string::*;
use debruijn::kmer::Kmer16;
use debruijn::Kmer;
use debruijn::Mer;
use fxhash::{FxHashMap, FxHashSet};
use rayon::prelude::*;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::mem;
use std::sync::{Arc, Mutex};
use std::time::Instant;
use std::string::String;

fn main() {
    let matches = Command::new("meta-cmg")
        .arg_required_else_help(true)
        .version("0.1")
        .about("Chromatic minimizer pangenome graph program.")
        .subcommand(
            Command::new("generate")
                .about("Generate graph.")
                .version("0.1")
                .arg(
                    Arg::new("references")
                        .index(1)
                        .help("Input reference fasta files. Any sequence within every reference is assumed to be homologous with one another.")
                        .required(true)
                        .action(ArgAction::Append),
                ).
                arg(
                    Arg::new("output")
                        .short('o')
                        .default_value("serialized_mini_graph")
                        .help("Name of output chromatic reference graph. Produces a .json and .bin file. (Default: serialized_mini_graph)")
                ).
                arg(
                    Arg::new("syncmer")
                        .short('s')
                        .help("Use syncmers. (Default: use minimizers)")
                        .action(ArgAction::SetTrue)
                ).
                arg(
                    Arg::new("circular")
                        .short('c')
                        .help("Assume that the genomes are circular. (Default: not circular)")
                        .action(ArgAction::SetTrue)
                ).
                arg(
                    Arg::new("mask")
                        .short('m')
                        .default_value("0.0002")
                        .help("Mask fraction of k-mers (Default: top 0.0002 of repetitive k-mers)")
                ).
                arg(
                    Arg::new("samp_freq")
                        .short('f')
                        .default_value("30")
                        .help("Graph sampling frequency for strain detect. (Default: 30)")
                ).
                arg(
                    Arg::new("minimizer_weighting")
                        .short('r')
                        .help("marbl minimizer weighting file for k = 16. (Default: none)")
                ).
                arg(
                    Arg::new("w")
                        .short('w')
                        .default_value("16")
                        .help("w value (testing).")
                        .hide(true)
                ).
                arg(
                    Arg::new("h")
                        .short('H')
                        .default_value("50")
                        .help("h value (testing).")
                        .hide(true)
                )
        )
        .subcommand(
            clap::Command::new("map")
                .about("Map sequences onto graph.")
                .version("0.1")
                .arg(
                    Arg::new("reference_graph")
                        .required(true)
                        .index(1)
                        .help("Reference graph (.bin) output from the generate subcommand. E.g. serialized_mini_graph.bin"),
                ).
                arg(
                    Arg::new("reads")
                        .required(true)
                        .index(2)
                        .help("Reads to map to graph. Only support fastq right now."),
                ).
                arg(Arg::new("threads")
                  .short('t')
                  .default_value("10")
                  .help("Number of threads to use. (default: 10).")
                  .value_name("INT")
                ).
                arg(
                    Arg::new("syncmer")
                        .short('s')
                        .help("Use syncmers. (Default: use minimizers)")
                        .action(ArgAction::SetTrue)
                ).
                arg(
                    Arg::new("chain_heuristic")
                        .short('d')
                        .help("Use linearization heuristic instead of DAG-aware heuristic. (Default: use chain heuristic)")
                        .action(ArgAction::SetTrue)
                ).
                arg(
                    Arg::new("dont_output_stuff")
                        .short('u')
                        .help("Output auxillary info (Default: on)")
                        .action(ArgAction::SetTrue)
                ).
                arg(
                    Arg::new("align")
                        .short('a')
                        .help("Output alignment in BAM format (Default: no alignment)")
                        .action(ArgAction::SetTrue)
                ).
                arg(
                    Arg::new("bam_name")
                        .short('b')
                        .default_value("output.bam")
                        .help("Name of output bam file. (Default: output.bam)")
                ).
                arg(
                    Arg::new("short_reads")
                        .short('x')
                        .help("Short read parameters. (Default: Long read params.)")
                        .action(ArgAction::SetTrue)
                        .hide(true)
                ).
                arg(
                    Arg::new("penalty")
                        .short('p')
                        .help("Short read parameters. (Default: Long read params.)")
                        .hide(true)
                ).
                arg(
                    Arg::new("minimizer_weighting")
                        .short('r')
                        .help("marbl minimizer weighting file for k = 16. (Default: none)")
                ).
                arg(
                    Arg::new("h")
                        .short('H')
                        .default_value("50")
                        .help("h value (testing).")
                        .hide(true)
                ).
                arg(
                    Arg::new("w")
                        .short('w')
                        .default_value("16")
                        .help("Value of w for minimizer window size. (Default: 16)")
                ).
                arg(
                    Arg::new("trace")
                        .long("trace")
                        .action(ArgAction::SetTrue)
                        .help("Output debug information for tracing.")
                )
        )
        .get_matches();

    let generate;
    let matches_subc: &ArgMatches;
    if let Some(matches) = matches.subcommand_matches("generate") {
        generate = true;
        matches_subc = matches;
    } else {
        matches_subc = matches.subcommand_matches("map").unwrap();
        generate = false;
    }

    let circular;
    if matches_subc.try_contains_id("circular").is_ok() && matches_subc.get_flag("circular") {
        circular = true;
    } else {
        circular = false;
    }

    let align;
    if matches_subc.try_contains_id("align").is_ok() && matches_subc.get_flag("align") {
        align = true;
    } else {
        align = false;
    }

    let k = 16;
    let s = 8;
    let t = (k - s + 2) / 2 as usize;

    // use syncmers if not using minimizers
    let use_minimizers;
    if matches_subc.get_flag("syncmer") {
        use_minimizers = false;
    } else {
        use_minimizers = true;
    }


    if matches_subc.try_contains_id("trace").is_ok() && matches_subc.get_flag("trace") {
        SimpleLogger::new().with_level(LevelFilter::Trace).env().init().unwrap();
    }
    else{
        SimpleLogger::new().with_level(LevelFilter::Debug).env().init().unwrap();
    }

    let chain_heuristic;
    if matches_subc.try_contains_id("chain_heuristic").is_ok() && matches_subc.get_flag("chain_heuristic") {
        chain_heuristic = false;
    } else {
        chain_heuristic = true;
    }
    
    let h = matches_subc
        .get_one::<String>("h")
        .unwrap()
        .as_str()
        .parse::<usize>()
        .unwrap();

    let w = matches_subc
        .get_one::<String>("w")
        .unwrap()
        .as_str()
        .parse::<usize>()
        .unwrap();

    let minimizer_weight_file = matches_subc.get_one::<String>("minimizer_weighting");
    let frequent_kmers;
    if let Some(file_str) = minimizer_weight_file {
        frequent_kmers = seeding_methods_bit::read_minimizer_count_file(file_str);
    } else {
        frequent_kmers = FxHashMap::default();
    }

    if generate {
        let fraction_mask = matches_subc.get_one::<String>("mask").unwrap().as_str();
        let fraction_mask_f64: f64 = fraction_mask.parse().unwrap();

        let samp_freq = matches_subc
            .get_one::<String>("samp_freq")
            .unwrap()
            .as_str()
            .parse::<usize>()
            .unwrap();

        let ref_genomes: Vec<&str> = matches_subc.get_many::<String>("references").unwrap().map(|x| x.as_str()).collect();
        let mut chroms = vec![];
        let mut good_chroms = vec![];
        let mut chrom_names = vec![];
        let mut good_chrom_names = vec![];

        let overall_start = Instant::now();

        for i in 0..ref_genomes.len() {
            let reader = fasta::Reader::from_file(&ref_genomes[i]);
            for record in reader.unwrap().records() {
                let rec = record.unwrap();
                let rec_desc = rec.desc();
                if let Some(rec_desc_str) = rec_desc {
                    if rec_desc_str.contains("plasmid") {
                        continue;
                    }
                }
                println!(
                    "Iteration: {}, Contig: {}, Reference: {}.",
                    chroms.len(),
                    rec.id(),
                    ref_genomes[i]
                );
                chrom_names.push(rec.id().to_string());
                let chrom = DnaString::from_acgt_bytes(rec.seq());
                chroms.push((chrom, true));
            }
        }

        good_chroms.push((chroms[0].0.clone(), true));
        // good_chroms.push((DnaString::new(), true));
        good_chrom_names.push(chrom_names[0].clone());

        let mut seeds1;
        let seed_p1;
        let dont_use_kmers = seeding_methods_bit::get_masked_kmers(
            &chroms[0].0,
            w,
            k,
            s,
            t,
            fraction_mask_f64,
            use_minimizers,
            &frequent_kmers,
            circular,
        );
        println!("----------------------------------------");
        println!("w = {}, k = {}, s = {}, t = {}", w, k, s, t);
        println!("Sample frequency: {}", samp_freq);
        println!("Fraction of masked kmers: {}", fraction_mask);
        println!("Use minimizers: {}", use_minimizers);
        println!("Use syncmers: {}", !use_minimizers);
        println!("Circular genome: {}", circular);
        println!("h value: {}", h);
        println!("Number of chromosomes: {}", chroms.len());
        println!("Number of masked kmers: {}", dont_use_kmers.len());
        println!("----------------------------------------");
        if use_minimizers {
            seed_p1 = seeding_methods_bit::minimizer_seeds(
                &chroms[0].0,
                w,
                k,
                samp_freq,
                &dont_use_kmers,
                &frequent_kmers,
                true,
                circular
            );
        } else {
            seed_p1 = seeding_methods_bit::open_sync_seeds(
                &chroms[0].0,
                k,
                t,
                s,
                samp_freq,
                &dont_use_kmers,
                &frequent_kmers,
                true,
            );
        }
        seeds1 = seed_p1.0;
        let p1 = seed_p1.1;

        println!(
            "Starting reference is {} and has {} nodes.",
            ref_genomes[0],
            seeds1.len()
        );

        let mut aln_score_array = vec![];
        let mut mean_score = 0.0;
        graph_utils::top_sort_kahns(&mut seeds1);

        for i in 1..chroms.len() {
            let old_graph_len = seeds1.len();
            let genome_string = &chroms[i].0;

            println!("-----------------Iteration {}-------------------", i);
            let mut seeds2;
            let now = Instant::now();
            let s2;
            if use_minimizers {
                s2 = seeding_methods_bit::minimizer_seeds(
                    genome_string,
                    w,
                    k,
                    1,
                    &dont_use_kmers,
                    &frequent_kmers,
                    false,
                    circular
                );
            } else {
                s2 = seeding_methods_bit::open_sync_seeds(
                    genome_string,
                    k,
                    t,
                    s,
                    1,
                    &dont_use_kmers,
                    &frequent_kmers,
                    false,
                );
            }
            seeds2 = s2.0;
            println!(
                "Generating sketch (minimizers) time: {}",
                now.elapsed().as_secs_f32()
            );
            let now = Instant::now();
            let ref_hash_map = chain::get_kmer_dict_mut(&mut seeds1);
            let q_hash_map = chain::get_kmer_dict(&seeds2);

            let mut kmer_count_dict = FxHashMap::default();
            for node in seeds1.iter() {
                let count_num = kmer_count_dict.entry(node.kmer).or_insert(0);
                *count_num += 1;
            }

            let mut hash_vec: Vec<(&Kmer16, &usize)> = kmer_count_dict.iter().collect();
            hash_vec.sort_by(|a, b| b.1.cmp(a.1));
            let qlen = seeds2.len();

            // chain seeds and get the best anchors
            let anc_score_strand_vec = chain::chain_seeds(
                &mut seeds1,
                &mut seeds2,
                &ref_hash_map,
                &q_hash_map,
                h,
                chain_heuristic,
                false,
                &dont_use_kmers,
                circular,
                None,
                None,
                None
            );

            let (best_anchors, aln_score, forward_strand) = anc_score_strand_vec
                .into_iter()
                .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
                .unwrap();
            
            //Need to reverse the read strand so that it is "forward". Bad mutability
            //design pattern here will change TODO
            //Only need this for circula because circular does chaining for both strands,
            //hence mutates the state back to normal. Needs to be reversed if
            //reverse is the best strand. Non-circular already reverses during the chaining.
            let order_val_last = seeds2[qlen - 1].order_val;
            if forward_strand == false && circular {
                for node in seeds2.iter_mut() {
                    node.order = qlen as u32 - node.order - 1;
                    //            for child_id in node.child_nodes.iter_mut(){
                    //                *child_id = q_len as u32 - *child_id - 1;
                    //            }
                }
            }

            if !forward_strand {
                chroms[i].1 = false;
            }

            println!(
                "Chaining time and aln_score and strand: {},{},{}",
                now.elapsed().as_secs_f32(),
                aln_score,
                forward_strand
            );

            println!("Aln score, mean score {},{}", aln_score, mean_score);
            if aln_score < 0.75 * mean_score && circular {
                println!("Bad alignment. Continuing");
                continue;
            }
            // TODO: use actual chromosome data instead of empty string in production
            good_chroms.push((chroms[i].0.clone(), forward_strand));
            // good_chroms.push((DnaString::new(), true));
            good_chrom_names.push(chrom_names[i].clone());

            mean_score = (mean_score * (i - 1) as f64 + aln_score) / (i as f64);
            aln_score_array.push(aln_score);

            let now = Instant::now();
            graph_utils::add_align_to_graph(
                &mut seeds1,
                seeds2,
                best_anchors,
                forward_strand,
                samp_freq,
                circular,
            );
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

            graph_utils::top_sort_kahns(&mut seeds1);
            println!("Top sort time: {}.", now.elapsed().as_secs_f32());
        }

        // compute bubbles
        let mut bubbles: Vec<Bubble> = vec![];
        let mut cmg: Cmg<'_> = Cmg::new(&mut seeds1, k as u8);
        lsd::detector::detect(&mut cmg, &mut bubbles);
        
        let now = Instant::now();
        let top_sort = graph_utils::top_sort(&seeds1);
        for bubble in bubbles.iter_mut() {
            bubble.shortest_path_length = Some(graph_utils::shortest_path_length(&seeds1, &top_sort, &bubble.start, &bubble.end));
            bubble.longest_path_length = Some(graph_utils::longest_path_length(&seeds1, &top_sort, &bubble.start, &bubble.end));
        }
        println!("Computing shortest and longest branches on bubbles time: {}", now.elapsed().as_secs_f32());
        // bubbles_utils::find_superbubbles(&seeds1, &top_sort, &mut bubbles);
        // println!("Found {} bubbles in time {}", bubbles.len(), now.elapsed().as_secs_f32());
        // let mut file = File::create("bubbles.csv").unwrap();
        // write!(
        //     &mut file,
        //     "id,start_node_id,end_node_id,shortest_path,longest_path\n",
        // ).unwrap();
        // for bubble in bubbles.iter() {
        //     let to_write = format!(
        //         "{},{},{},{},{}\n",
        //         bubble.id,
        //         bubble.start,
        //         bubble.end,
        //         bubble.shortest_path_length.unwrap(),
        //         bubble.longest_path_length.unwrap()
        //     );
        //     write!(&mut file, "{}", to_write).unwrap();
        // }
        graph_utils::get_closest_node(&mut seeds1);
        // let mut dist_mat = SparseMatrix::new((seeds1.len(), seeds1.len()));
        // for bubble in bubbles.iter() {
        //     let shortest_path_on_bubble = bubble.shortest_path_length.unwrap();
        //     let longest_path_on_bubble = bubble.longest_path_length.unwrap();
        //     if longest_path_on_bubble as f32 > 1.5 * shortest_path_on_bubble as f32 {
        //         // save estimated reference distance in a sparse matrix
        //         // for i in 0..bubble.kmers.len() {
        //         //     let node1 = &seeds1[bubble.kmers[i] as usize];
        //         //     let closest1 = &closest_node[bubble.kmers[i] as usize];
        //         //     if closest1.is_none() {
        //         //         continue;
        //         //     }
        //         //     let closest1_node = &seeds1[closest1.unwrap().0 as usize];
        //         //     let dist_to_closest1 = closest1.unwrap().1;
    
        //         //     for j in i+1..bubble.kmers.len() {
        //         //         let node2 = &seeds1[bubble.kmers[j] as usize];
        //         //         let closest2 = &closest_node[bubble.kmers[j] as usize];
        //         //         if closest2.is_none() {
        //         //             continue;
        //         //         }
        //         //         let closest2_node = &seeds1[closest2.unwrap().0 as usize];
        //         //         let dist_to_closest2 = closest2.unwrap().1;
        //         //         // find closest distance between two closest ref nodes
        //         //         let mut ref_dist = u16::MAX;
        //         //         for u in 0..closest1_node.actual_ref_positions.len() {
        //         //             for v in 0..closest2_node.actual_ref_positions.len() {
        //         //                 let closest1_offset = closest1_node.actual_ref_positions[u] + (dist_to_closest1 as usize * k);
        //         //                 let closest2_offset = closest2_node.actual_ref_positions[v] + (dist_to_closest2 as usize * k);
        //         //                 let dist = (closest1_offset.abs_diff(closest2_offset)) as u16;
        //         //                 if dist < ref_dist {
        //         //                     ref_dist = dist;
        //         //                 }
        //         //             }
        //         //         }
        //         //         let linearized_dist = node1.order_val.abs_diff(node2.order_val) as u16;
        //         //         if linearized_dist.abs_diff(ref_dist) <= (longest_path_on_bubble - shortest_path_on_bubble) as u16 {
        //         //             continue;
        //         //         }
        //         //         dist_mat.set(node1.id as usize, node2.id as usize, ref_dist as u32);
        //         //         // println!("{}-{}: {}-{}: {}-{}", node1.id, node2.id, node1.order_val, node2.order_val, ref_dist, linearized_dist);
        //         //     }
        //         // }
        //         // save estimated reference distance between the start and end of the bubble
        //         let node1 = &seeds1[bubble.start as usize];
        //         let closest1 = &closest_node[bubble.start as usize];
        //         if closest1.is_none() {
        //             continue;
        //         }
        //         let closest1_node = &seeds1[closest1.unwrap().0 as usize];
        //         let dist_to_closest1 = closest1.unwrap().1;
        //         let node2 = &seeds1[bubble.end as usize];
        //         let closest2 = &closest_node[bubble.end as usize];
        //         if closest2.is_none() {
        //             continue;
        //         }
        //         let closest2_node = &seeds1[closest2.unwrap().0 as usize];
        //         let dist_to_closest2 = closest2.unwrap().1;
        //         // find closest distance between two closest ref nodes
        //         let mut ref_dist = u16::MAX;
        //         for u in 0..closest1_node.actual_ref_positions.len() {
        //             for v in 0..closest2_node.actual_ref_positions.len() {
        //                 let closest1_offset = closest1_node.actual_ref_positions[u] + (dist_to_closest1 as usize * k);
        //                 let closest2_offset = closest2_node.actual_ref_positions[v] + (dist_to_closest2 as usize * k);
        //                 let dist = (closest1_offset.abs_diff(closest2_offset)) as u16;
        //                 if dist < ref_dist {
        //                     ref_dist = dist;
        //                 }
        //             }
        //         }
        //         let linearized_dist = node1.order_val.abs_diff(node2.order_val) as u16;
        //         if linearized_dist.abs_diff(ref_dist) <= (longest_path_on_bubble - shortest_path_on_bubble) as u16 {
        //             continue;
        //         }
        //         dist_mat.set(node1.id as usize, node2.id as usize, ref_dist as u32);
        //     }
        // }
        // println!("{} element set in sparse distance matrix", dist_mat.num_nonzero);
        let now = Instant::now();
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
            .get_one::<String>("output")
            .unwrap()
            .as_str();
        let serial_json_name = format!("{}.json", serial_name);
        let serial_bin_name = format!("{}.bin", serial_name);

        let mut file_json = File::create(serial_json_name).unwrap();
        write!(&mut file_json, "{}", j.unwrap()).unwrap();

        
        let mut file_bin = BufWriter::new(File::create(serial_bin_name).unwrap());
        // bincode::serialize_into(
        //     &mut file_bin,
        //     &(&seeds1, &good_chroms, &good_chrom_names, &dont_use_kmers, &bubbles),
        // )
        // .unwrap();
        bincode::serialize_into(
            &mut file_bin,
            &(&seeds1, &bubbles),
        )
        .unwrap();
        println!(
            "Serializing and writing time {}.",
            now.elapsed().as_secs_f32()
        );

        // let mut file = File::create("full_mini_graph.csv").unwrap();
        // for node in seeds1.iter() {
        //     //            println!("{:?},{}", node.kmer, node.kmer.to_string());
        //     for child in node.child_nodes.iter() {
        //         let towrite = format!(
        //             "{}-{},{}-{}\n",
        //             node.order, node.id, seeds1[*child as usize].order, seeds1[*child as usize].id,
        //         );
        //         write!(&mut file, "{}", towrite).unwrap();
        //     }
        // }
        // // save topological ordering to csv
        // println!("Saving topological ordering to csv");
        // let mut file = File::create("topo_order.csv").unwrap();
        // for node in seeds1.iter() {
        //     let towrite = format!("{},{}\n", node.id, node.order_val);
        //     write!(&mut file, "{}", towrite).unwrap();
        // }
        println!("Total time {}", overall_start.elapsed().as_secs_f32());
    } else {
        let num_t_str = matches_subc.get_one::<String>("threads").unwrap();
        let num_t = match num_t_str.parse::<usize>() {
            Ok(num_t) => num_t,
            Err(_) => panic!("Number of threads must be positive integer"),
        };
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_t)
            .build_global()
            .unwrap();

        let dont_output_stuff = matches_subc.get_flag("dont_output_stuff");
        let short_reads = matches_subc.get_flag("short_reads");
        let ref_graph_file = matches_subc.get_one::<String>("reference_graph").unwrap().as_str();
        let bam_name = matches_subc.get_one::<String>("bam_name").unwrap().as_str();

        let ref_graph_f = File::open(ref_graph_file).unwrap();
        let ref_graph_reader = BufReader::new(ref_graph_f);

        let (mut ref_graph, chroms, chrom_names, dont_use_kmers, bubbles): (
            Vec<KmerNode>,
            Vec<(DnaString, bool)>,
            Vec<String>,
            FxHashSet<Kmer16>,
            Vec<Bubble>,
        ) = bincode::deserialize_from(ref_graph_reader).unwrap();

        let closest_bubble_source = graph_utils::get_closest_bubble_source(&ref_graph, &bubbles);

        let order_to_id = graph_utils::top_sort_kahns(&mut ref_graph);
        let reads_file = matches_subc.get_one::<String>("reads").unwrap().as_str();
        let reader = fastq::Reader::from_file(reads_file);

        let ref_hash_map = chain::get_kmer_dict(&ref_graph);
        let anchor_file = Arc::new(Mutex::new(BufWriter::new(File::create("read_anchor_hits.txt").unwrap())));
        let mut best_genomes_file = File::create("best_genome_reads.txt").unwrap();
        let (headerview, mut writer) =
            align::write_bam_header(&chroms, &chrom_names, bam_name.to_string());

        //        let graph_simp_time = Instant::now();
        //        let simplified_graph = graph_utils::concat_graph(&ref_graph[0], &ref_graph);
        //        println!("Graph simp time {}", graph_simp_time.elapsed().as_secs_f32());

        let start_align = Instant::now();
        let mut num_mapped_reads = Mutex::new(0);
        let mut record_container = vec![];
        let mut bam_info_container: Mutex<Vec<_>> = Mutex::new(vec![]);
        let mut best_hit_for_read: Mutex<FxHashMap<_, _>> = Mutex::new(FxHashMap::default());

        let mut records = reader.unwrap().records().peekable();
        //        let batch = num_t * 200;
        let batch = usize::MAX;
        while let Some(Ok(record)) = records.next() {
            if record_container.len() < batch {
                record_container.push(record);
            }
            if record_container.len() == batch || records.peek().is_none() {
                (0..record_container.len())
                    .collect::<Vec<usize>>()
                    .into_par_iter()
                    .for_each(|i| {
                        let now = Instant::now();
                        let total_time = Instant::now();
                        let rec = &record_container[i];
                        let read = DnaString::from_acgt_bytes(rec.seq());
                        let quals = rec.qual().to_vec();
                        let read_id = rec.id().to_string();

                        let mut best_colors_both_strands = vec![];
                        let mut best_anchors_both_strands = vec![];
                        let mut chain_numbers = vec![];
                        let mut strand_anchor_vec = vec![];

                        log::trace!("Preprocess time: {}", now.elapsed().as_secs_f32());
                        log::trace!("---------------Read: {}---------------", read_id);
                        let now = Instant::now();
                        let mut read_seeds;
                        let s2;
                        if use_minimizers {
                            s2 = seeding_methods_bit::minimizer_seeds(
                                &read,
                                w,
                                k,
                                1,
                                &FxHashSet::default(),
                                &frequent_kmers,
                                false,
                                circular,
                            );
                        } else {
                            s2 = seeding_methods_bit::open_sync_seeds(
                                &read,
                                k,
                                t,
                                s,
                                1,
                                &FxHashSet::default(),
                                &frequent_kmers,
                                false,
                            );
                        }
                        log::trace!("Seeding time: {}", now.elapsed().as_secs_f32());

                        read_seeds = s2.0;
                        let qlen = read_seeds.len();
                        let q_hash_map = chain::get_kmer_dict(&read_seeds);
                        let now = Instant::now();
                        let anc_score_strand_vec = chain::chain_seeds(
                            &ref_graph,
                            &mut read_seeds,
                            &ref_hash_map,
                            &q_hash_map,
                            h,
                            chain_heuristic,
                            true,
                            &FxHashSet::default(),
                            circular,
                            Some(&closest_bubble_source),
                            Some(&bubbles),
                            Some(&ref_graph)
                        );

                        log::trace!("Chaining time: {}", now.elapsed().as_secs_f32());

                        for (k, (best_anchors, _aln_score, read_strand)) in
                            anc_score_strand_vec.iter().enumerate()
                        {
                            if *read_strand == false {
                                for node in read_seeds.iter_mut() {
                                    node.order = qlen as u32 - node.order - 1;
                                }
                            }
                            let now = Instant::now();

                            let (best_colors, best_list_anchors) = chain::get_best_path_from_chain_rewrite(
                                best_anchors,
                                &ref_graph,
                                &order_to_id,
                                &read_seeds,
                                read.len(),
                                false,
                            );

                            log::trace!("Path collection time: {}", now.elapsed().as_secs_f32());

                            //TODO don't want to do clone every anchor list.
                            for i in 0..best_colors.len() {
                                best_colors_both_strands.push(best_colors[i]);
                                best_anchors_both_strands.push(best_list_anchors[i].clone());
                                strand_anchor_vec.push(*read_strand);
                                chain_numbers.push(k);
                            }
                        }

                        if align {
                            let mut do_align = true;
                            let do_base_chain = true;
                            if best_anchors_both_strands.len() == 0 && false{
                                if read.len() > constants::READ_LENGTH_SUPER_CHAIN_CUTOFF
                                    && do_base_chain
                                {
                                    log::trace!("No good alignment found; primary-ref-chaining");
                                    let now = Instant::now();
                                    //Do base chaining if no good alignment is found
                                    let mut base_anc_score_strand_vec =
                                        coord_chain::get_base_chains(
                                            &ref_graph,
                                            &read_seeds,
                                            &ref_hash_map,
                                            &q_hash_map,
                                            h,
                                            &FxHashSet::default(),
                                            read.len(),
                                        );

                                    log::trace!(
                                        "primary-ref chaining time: {}",
                                        now.elapsed().as_secs_f32()
                                    );

                                    for (k, (best_anchors, _aln_score, read_strand)) in
                                        base_anc_score_strand_vec.iter().enumerate()
                                    {
                                        if *read_strand == false {
                                            for node in read_seeds.iter_mut() {
                                                node.order = qlen as u32 - node.order - 1;
                                                //            for child_id in node.child_nodes.iter_mut(){
                                                //                *child_id = q_len as u32 - *child_id - 1;
                                                //            }
                                            }
                                        }
                                        let now = Instant::now();

                                        let (best_colors, best_list_anchors) =
                                            chain::get_best_path_from_chain_rewrite(
                                                best_anchors,
                                                &ref_graph,
                                                &order_to_id,
                                                &read_seeds,
                                                read.len(),
                                                true,
                                            );

                                            log::trace!(
                                                "Path collection time: {}",
                                                now.elapsed().as_secs_f32()
                                            );

                                        //TODO don't want to do clone every anchor list.
                                        for i in 0..best_colors.len() {
                                            best_colors_both_strands.push(best_colors[i]);
                                            best_anchors_both_strands
                                                .push(best_list_anchors[i].clone());
                                            chain_numbers.push(k);
                                            strand_anchor_vec.push(*read_strand);
                                        }

                                        if best_anchors_both_strands.len() > 0 {
                                            do_align = true;
                                        } else {
                                            do_align = false;
                                        }
                                    }
                                } else {
                                    log::trace!("No good alignment found");
                                    do_align = false;
                                }
                            }
                            
                            if do_align {
                                let now = Instant::now();
                                use std::cmp::min;
                                //Print best anchors
                                let mut best_indices = best_anchors_both_strands
                                    .iter()
                                    .enumerate()
                                    .collect::<Vec<(_, _)>>();
                                best_indices
                                    .sort_by(|(_, a), (_, b)| b.1.partial_cmp(&a.1).unwrap());
                                let best_indices: Vec<usize> =
                                    best_indices.iter().map(|(index, _)| *index).collect();
                                let top_n = 10;
                                //                                writeln!(&mut best_genomes_file, ">{}", &read_id).unwrap();
                                for _i in 0..min(top_n, best_indices.len()) {
                                    let mut locked = best_hit_for_read.lock().unwrap();
                                    let ith_color = best_colors_both_strands[best_indices[_i]];
                                    let ith_score = best_anchors_both_strands[best_indices[_i]].1;
                                    let ith_ref_chroms = align::get_nonzero_bits(ith_color);
                                    for bit in ith_ref_chroms {
                                        //                                        writeln!(
                                        //                                            &mut best_genomes_file,
                                        //                                            "{}\t{}",
                                        //                                            &chrom_names[chroms.len() - bit - 1],
                                        //                                            ith_score
                                        //                                        )
                                        //                                        .unwrap();
                                        let vec = locked.entry(read_id.clone()).or_insert(vec![]);
                                        vec.push((ith_score, &chrom_names[chroms.len() - bit - 1]))
                                    }
                                }

                                let map_indices;
                                let map_all = true;
                                if map_all {
                                    map_indices = top_n;
                                } else {
                                    map_indices = 1;
                                }

                                let now = Instant::now();
                                let mut used_chains = FxHashSet::default();
                                for index in 0..usize::min(best_indices.len(), map_indices) {
                                    let best_index = best_indices[index];
                                    let best_score = best_anchors_both_strands[best_indices[0]].1;
                                    let score = best_anchors_both_strands[best_index].1;
                                    let chain_number = chain_numbers[best_index];
                                    if (best_score < score + 50.
                                        || (best_score - score) < 0.05 * (&best_anchors_both_strands[best_index].0.len() * 16) as f64)
                                        && !used_chains.contains(&chain_number)
                                    {
                                        used_chains.insert(chain_number);
                                        let anchors = &best_anchors_both_strands[best_index].0;
                                        let color = &best_colors_both_strands[best_index];
                                        let read_strand = strand_anchor_vec[best_index];

                                        let bam_info = align::align_from_chain(
                                            anchors,
                                            &chroms,
                                            *color,
                                            &ref_graph,
                                            &read_seeds,
                                            &read,
                                            read_strand,
                                            &quals,
                                            &chrom_names,
                                            &read_id,
                                        );
                                        let mut locked = bam_info_container.lock().unwrap();
                                        let mut l2 = num_mapped_reads.lock().unwrap();
                                        *l2 += 1;
                                        if *l2 %  10000 == 0{
                                            log::info!("{} mapped reads", *l2);
                                        }
                                        locked.push(bam_info);
                                    }
                                }
                                log::trace!("Total align time {}", now.elapsed().as_secs_f32());
                            }
                        }
                        log::trace!(
                            "Total time mapping read {} is {}",
                            read_id,
                            total_time.elapsed().as_secs_f32()
                        );
                        if !dont_output_stuff {
                            for (color, anchors_color_pair) in best_anchors_both_strands.iter().enumerate() {
                                write!(
                                    anchor_file.lock().unwrap(),
                                    "{}:{}:",
                                    read_id,
                                    color
                                )
                                .unwrap();
                                for anchor in &anchors_color_pair.0 {
                                    if anchor == anchors_color_pair.0.last().unwrap() {
                                        writeln!(
                                            anchor_file.lock().unwrap(),
                                            "{}\n",
                                            anchor.0
                                        )
                                        .unwrap();
                                    } else {
                                        write!(
                                            anchor_file.lock().unwrap(),
                                            "{},",
                                            anchor.0
                                        )
                                        .unwrap();
                                    }
                                }
                            }
                        }
                    });
                for bam_info in bam_info_container.into_inner().unwrap() {
                    if !bam_info.is_none() {
                        let bam_rec = align::get_bam_record(bam_info.unwrap(), &headerview);
                        writer.write(&bam_rec).unwrap();
                    }
                }
                record_container = vec![];
                bam_info_container = Mutex::new(vec![]);
            }
        }
        if !dont_output_stuff{
            for (read_id, hits) in best_hit_for_read.into_inner().unwrap() {
                writeln!(&mut best_genomes_file, ">{}", &read_id).unwrap();
                for hit in hits {
                    writeln!(&mut best_genomes_file, "{}\t{}", hit.1, hit.0).unwrap();
                }
            }
        }
        log::trace!(
            "Alignment took {} seconds",
            start_align.elapsed().as_secs_f32()
        );
    }
}
