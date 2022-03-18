# chrom_mini_graph 

chrom_mini_graph is a tool for generating and mapping reads onto a chromatic (coloured) minimizer pangenome graph. 

### Requirements 

1. [rust](https://www.rust-lang.org/tools/install) and associated tools such as cargo are required and assumed to be in PATH.
2. cmake version >= 3.12 has to be in path due to dependency on --parallel command (see [here](https://githubmemory.com/repo/alexcrichton/cmake-rs/issues/131?page=1)).

### Install

```
git clone https://github.com/bluenote-1577/chrom_mini_graph
cd chrom_mini_graph
cargo build --release
./target/release/chrom_mini_graph generate test_refs/*
./target/release/chrom_mini_graph map -a -b test_bam.bam serialized_mini_graph.bin test_reads/hg_01243_pacbio_reads.fastq > output.txt
```

1. `cargo build --release` first builds the **chrom_mini_graph** binary, which is found in the ./target/release/ directory.
2. The `chrom_mini_graph generate` command generates a coloured minimizer pangenome graph. 
3. The `chrom_mini_graph map` command chains onto the output graph and produces an alignment. The `-a` option outputs a BAM file with name specified by the `-b` option.

* 6 reference 1M bp segments of chromosome 20 are provided in the test_ref folder. 
* Simulated PacBio CLR reads for hg01243 are available in the test_reads folder. 

# Using chrom_mini_graph

## generate

`chrom_mini_graph generate ref_1.fasta ref_2.fasta ... -o output_from_generate` to create a coloured minimizer pangenome graph for references ref_1.fasta, ref_2.fasta, etc. The output specified by the `-o` option is used for the mapping step. 

* The window size can be easily modified in the `src/bin/chrom_mini_graph.rs` file. The default value is 16.
* Outputs a \*.bin file to be used for mapping and other auxillary information; see below. 
* Each fasta file can have multiple contigs. Each contig will be treated as its own reference genome.

### Ordering for `generate`

The first reference used (i.e. `ref_1.fasta`) serves as the backbone for the minimizer graph. Make sure that this first reference is the most contiguous contig. 

## map

A proof of concept read-to-graph chainer by chaining minimizers in the read onto the graph without knowledge of colour and then finding the best colours (reference genomes) for the chain.

`chrom_mini_graph map -a output_from_generate.bin your_reads.fastq -b bam_name.bam > output.txt` outputs the bam file `bam_name.bam` and directs stdout to a output.txt log. 

The file *best_genome_reads.txt* is also output. The `best_genome_reads.txt` shows the top 5 (or less) best candidate reference genomes for each read.  The format is 
```
>read_1
chrom_1 score_1
chrom_2 score_2
...
>read_2
chrom_1 score_1
chrom_2 score_2
...
```

More than 5 best candidates may be output due to secondary alignments and less than 5 may be output if the read is deemed unmappable to certain references. 

Caveats:

* Only the best candidate genome is aligned to. 
* No supplementary/secondary alignments are output in the bam file; MapQ is defaulted to 60. 

<!--- This is the JSON serialization for the entire graph. See [this document](https://docs.google.com/document/d/1oRHjPgP-Bh9UkySCduWIl5yCpfiLVEoSnRdzdx4a7-Y/edit?usp=sharing) for how to deserialize the graph. **IMPORTANT:** For colouring, the most significant bit corresponds to the first genome in the command, and the least significant bit corresponds to the last genome in the command. --->
