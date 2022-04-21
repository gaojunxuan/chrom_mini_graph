use bincode;
use bio::io::{fasta, fastq};
use block_aligner::scan_block::*;
use block_aligner::scores::*;
use chrom_mini_graph::align;
use chrom_mini_graph::chain;
use chrom_mini_graph::data_structs::KmerNode;
use chrom_mini_graph::deconvolution;
use chrom_mini_graph::graph_utils;
use chrom_mini_graph::seeding_methods_bit;
use clap::{App, AppSettings, Arg, SubCommand};
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
use std::sync::Mutex;
use std::time::Instant;

fn main() {
    let matches = App::new("meta-cmg")
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
                        .help("Name of output chromatic reference graph. Produces a .json and .bin file. (Default: serialized_mini_graph)")
                        .takes_value(true),
                ).
                arg(
                    Arg::with_name("syncmer")
                        .short("s")
                        .help("Use syncmers. (Default: use minimizers)")
                ).
                arg(
                    Arg::with_name("circular")
                        .short("c")
                        .help("Assume that the genomes are circular. (Default: not circular)")
                ).
                arg(
                    Arg::with_name("mask")
                        .short("m")
                        .help("Mask fraction of k-mers (Default: top 0.05% of repetitive k-mers)")
                        .takes_value(true)
                ).
                arg(
                    Arg::with_name("samp_freq")
                        .short("f")
                        .help("Graph sampling frequency for k-mer positions. (Default: 100)")
                        .takes_value(true)
                        .hidden(true)
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
                        .help("Reference graph (.bin) output from the generate subcommand. E.g. serialized_mini_graph.bin"),
                ).
                arg(
                    Arg::with_name("reads")
                        .required(true)
                        .index(2)
                        .help("Reads to map to graph. Only support fastq right now."),
                ).
                arg(Arg::with_name("threads")
                  .short("t")
                  .help("Number of threads to use. (default: 10).")
                  .value_name("INT")
                  .takes_value(true)
                ).
                arg(
                    Arg::with_name("syncmer")
                        .short("s")
                        .help("Use syncmers. (Default: use minimizers)")
                ).
                arg(
                    Arg::with_name("chain_heuristic")
                        .short("d")
                        .help("Use linearization heuristic instead of DAG-aware heuristic. (Default: use chain heuristic)")
                ).
                arg(
                    Arg::with_name("dont_output_stuff")
                        .short("u")
                        .help("Output auxillary info (Default: on)")
                ).
                arg(
                    Arg::with_name("align")
                        .short("a")
                        .help("Output alignment in BAM format (Default: no alignment)")
                ).
                arg(
                    Arg::with_name("bam_name")
                        .short("b")
                        .help("Name of output bam file. (Default: output.bam)")
                        .takes_value(true),
                ).
                arg(
                    Arg::with_name("short_reads")
                        .short("x")
                        .help("Short read parameters. (Default: Long read params.)")
                        .hidden(true)
                ).
                arg(
                    Arg::with_name("deconvolve")
                        .short("D")
                        .help("Deconvolve stuff (Default: off)")
                        .hidden(true)
                ).
                arg(
                    Arg::with_name("penalty")
                        .short("p")
                        .help("Short read parameters. (Default: Long read params.)")
                        .takes_value(true)
                        .hidden(true)
                ).
                arg(
                    Arg::with_name("samp_freq_decon")
                        .short("f")
                        .help("Graph sampling frequency for strain detect. (Default: 30)")
                        .takes_value(true)
                        .hidden(true)
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

    let circular;
    if matches_subc.is_present("circular") {
        circular = true;
    } else {
        circular = false;
    }

    let align;
    if matches_subc.is_present("align") {
        align = true;
    } else {
        align = false;
    }
    let deconvolve;
    if matches_subc.is_present("deconvolve") {
        deconvolve = true;
    } else {
        deconvolve = false;
    }

    let w = 16;
    let k = 16;
    let mask_repet_on_generate = true;
    let s = 10;
    let t = (k - s + 2) / 2 as usize;
    let h = 50;
    let samp_freq = 30;
    //use syncmers if not using minimizers

    let use_minimizers;
    if matches_subc.is_present("syncmer") {
        use_minimizers = false;
    } else {
        use_minimizers = true;
    }

    let chain_heuristic;
    if matches_subc.is_present("chain_heuristic") {
        chain_heuristic = false;
    } else {
        chain_heuristic = true;
    }

    if generate {
        let fraction_mask = matches_subc.value_of("mask").unwrap_or("0.0002");
        let fraction_mask_f64: f64 = fraction_mask.parse().unwrap();

        let ref_genomes: Vec<&str> = matches_subc.values_of("references").unwrap().collect();
        let mut chroms = vec![];
        let mut good_chroms = vec![];
        let mut chrom_names = vec![];
        let mut good_chrom_names = vec![];

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
        good_chrom_names.push(chrom_names[0].clone());

        let mut seeds1;
        let seed_p1;
        let dont_use_kmers =
            seeding_methods_bit::get_masked_kmers(&chroms[0].0, w, k, t, fraction_mask_f64, use_minimizers);
        if use_minimizers {
            seed_p1 = seeding_methods_bit::minimizer_seeds(
                &chroms[0].0,
                w,
                k,
                samp_freq,
                &dont_use_kmers,
            );
        } else {
            seed_p1 = seeding_methods_bit::open_sync_seeds(
                &chroms[0].0,
                k,
                t,
                samp_freq,
                &dont_use_kmers,
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

        for i in 1..chroms.len() {
            let old_graph_len = seeds1.len();
            let genome_string = &chroms[i].0;

            println!("-----------------Iteration {}-------------------", i);
            let mut seeds2;
            let now = Instant::now();
            let s2;
            if use_minimizers {
                s2 = seeding_methods_bit::minimizer_seeds(genome_string, w, k, 1, &dont_use_kmers);
            } else {
                s2 = seeding_methods_bit::open_sync_seeds(genome_string, k, t, 1, &dont_use_kmers);
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

            //            let (best_anchors, aln_score, forward_strand) = chain::chain_seeds(
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
            good_chroms.push((chroms[i].0.clone(), forward_strand));
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

            graph_utils::top_sort(&mut seeds1);
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
        bincode::serialize_into(
            &mut file_bin,
            &(&seeds1, &good_chroms, &good_chrom_names, &dont_use_kmers),
        )
        .unwrap();
        println!(
            "Serializing and writing time {}.",
            now.elapsed().as_secs_f32()
        );

        let mut file = File::create("full_mini_graph.csv").unwrap();
        for node in seeds1.iter() {
            //            println!("{:?},{}", node.kmer, node.kmer.to_string());
            for child in node.child_nodes.iter() {
                let towrite = format!(
                    "{}-{},{}-{}\n",
                    node.order, node.id, seeds1[*child as usize].order, seeds1[*child as usize].id,
                );
                write!(&mut file, "{}", towrite).unwrap();
            }
        }
    } else {
        let num_t_str = matches_subc.value_of("threads").unwrap_or("10");
        let num_t = match num_t_str.parse::<usize>() {
            Ok(num_t) => num_t,
            Err(_) => panic!("Number of threads must be positive integer"),
        };
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_t)
            .build_global()
            .unwrap();

        let dont_output_stuff = matches_subc.is_present("dont_output_stuff");
        let short_reads = matches_subc.is_present("short_reads");
        let ref_graph_file = matches_subc.value_of("reference_graph").unwrap();
        let bam_name = matches_subc.value_of("bam_name").unwrap_or("output.bam");
        let penalty: f64 = matches_subc
            .value_of("penalty")
            .unwrap_or("0.1")
            .parse()
            .unwrap();

        let samp_freq_decon: usize = matches_subc
            .value_of("samp_freq_decon")
            .unwrap_or("30")
            .parse()
            .unwrap();

        let ref_graph_f = File::open(ref_graph_file).unwrap();
        let ref_graph_reader = BufReader::new(ref_graph_f);

        let (mut ref_graph, chroms, chrom_names, dont_use_kmers): (
            Vec<KmerNode>,
            Vec<(DnaString, bool)>,
            Vec<String>,
            FxHashSet<Kmer16>,
        ) = bincode::deserialize_from(ref_graph_reader).unwrap();

        let closest_kmer_vec = graph_utils::get_closest_node(&ref_graph);
        let order_to_id = graph_utils::top_sort(&mut ref_graph);
        let reads_file = matches_subc.value_of("reads").unwrap();
        let reader = fastq::Reader::from_file(reads_file);

        let ref_hash_map = chain::get_kmer_dict(&ref_graph);
        let mut anchor_file = BufWriter::new(File::create("read_anchor_hits.txt").unwrap());
        let mut best_genomes_file = File::create("best_genome_reads.txt").unwrap();
        let (headerview, mut writer) =
            align::write_bam_header(&chroms, &chrom_names, bam_name.to_string());

        //        let graph_simp_time = Instant::now();
        //        let simplified_graph = graph_utils::concat_graph(&ref_graph[0], &ref_graph);
        //        println!("Graph simp time {}", graph_simp_time.elapsed().as_secs_f32());

        let start_align = Instant::now();
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
                        let rec = &record_container[i];
                        let read = DnaString::from_acgt_bytes(rec.seq());
                        let quals = rec.qual().to_vec();
                        let read_id = rec.id().to_string();

                        let mut best_colors_both_strands = vec![];
                        let mut best_anchors_both_strands = vec![];
                        let mut strand_anchor_vec = vec![];
                        println!("---------------Read: {}---------------", read_id);
                        let mut read_seeds;
                        let s2;
                        if use_minimizers {
                            s2 = seeding_methods_bit::minimizer_seeds(
                                &read,
                                w,
                                k,
                                1,
                                &FxHashSet::default(),
                            );
                        } else {
                            s2 = seeding_methods_bit::open_sync_seeds(
                                &read,
                                k,
                                t,
                                1,
                                &FxHashSet::default(),
                            );
                        }

                        read_seeds = s2.0;
                        let qlen = read_seeds.len();
                        let q_hash_map = chain::get_kmer_dict(&read_seeds);
                        let now = Instant::now();
                        //            let (best_anchors, _aln_score, read_strand) = chain::chain_seeds(
                        //
                        chain::get_super_chains(
                            &ref_graph,
                            &read_seeds,
                            &ref_hash_map,
                            &q_hash_map,
                            h,
                            chain_heuristic,
//                            &dont_use_kmers,
                                            &FxHashSet::default(),
                            &closest_kmer_vec,
                            samp_freq,
                            read.len()
                        );
                        println!("Super chaining time: {}", now.elapsed().as_secs_f32());
                        let now = Instant::now();
                        let mut anc_score_strand_vec = chain::chain_seeds(
                            &ref_graph,
                            &mut read_seeds,
                            &ref_hash_map,
                            &q_hash_map,
                            h,
                            chain_heuristic,
                            true,
//                            &dont_use_kmers,
                                            &FxHashSet::default(),
                            circular,
                        );
                        
                        println!("Chaining time: {}", now.elapsed().as_secs_f32());
//                        panic!();

                        for (best_anchors, _aln_score, read_strand) in anc_score_strand_vec.iter() {
                            //                            if !dont_output_stuff {
                            //                                write!(&mut anchor_file, "{}:", read_id).unwrap();
                            //                                for anchor in best_anchors.iter() {
                            //                                    write!(&mut anchor_file, "{},", anchor.0).unwrap();
                            //                                }
                            //                                write!(&mut anchor_file, "\n").unwrap();
                            //                            }

                            //Need to reverse the read strand so that it is "forward". Bad mutability
                            //design pattern here will change TODO. This is needed in alignment for
                            //retrieving positions and in graph construction too.
                            if *read_strand == false {
                                for node in read_seeds.iter_mut() {
                                    node.order = qlen as u32 - node.order - 1;
                                    //            for child_id in node.child_nodes.iter_mut(){
                                    //                *child_id = q_len as u32 - *child_id - 1;
                                    //            }
                                }
                            }
                            let now = Instant::now();

                            let (best_colors, best_list_anchors) = chain::get_best_path_from_chain2(
                                best_anchors,
                                &ref_graph,
                                &order_to_id,
                                &read_seeds,
                            );

                            println!("Path collection time: {}", now.elapsed().as_secs_f32());

                            //TODO don't want to do clone every anchor list.
                            for i in 0..best_colors.len() {
                                best_colors_both_strands.push(best_colors[i]);
                                best_anchors_both_strands.push(best_list_anchors[i].clone());
                                strand_anchor_vec.push(*read_strand);
                            }
                        }

                        if align {
                            if best_anchors_both_strands.len() == 0 {
                                println!("No good alignment found");
                            } else {
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
                                let top_n = 5;
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
                                let map_all = false;
                                if map_all {
                                    map_indices = top_n;
                                } else {
                                    map_indices = 1;
                                }
                                for index in 0..usize::min(best_indices.len(), map_indices) {
                                    let best_index = best_indices[index];
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
                                    locked.push(bam_info);
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
        for (read_id, hits) in best_hit_for_read.into_inner().unwrap() {
            writeln!(&mut best_genomes_file, ">{}", &read_id).unwrap();
            for hit in hits {
                writeln!(&mut best_genomes_file, "{}\t{}", hit.1, hit.0).unwrap();
            }
        }
        println!(
            "Alignment took {} seconds",
            start_align.elapsed().as_secs_f32()
        );
    }
}
