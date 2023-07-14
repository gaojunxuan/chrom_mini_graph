#!/bin/bash
#$ -P simpsonlab
#$ -cwd
#$ -l h_vmem=64G
#$ -l h_rt=144000
#$ -m eas
#$ -M jgao@oicr.on.ca

CMG_FOLDER=/u/jgao/jgao/code/chrom_mini_graph
CMG=/u/jgao/jgao/code/chrom_mini_graph/target/debug/chrom_mini_graph
GENERATED_FILE=/u/jgao/jgao/code/chrom_mini_graph/report/$1/serialized_mini_graph.bin
SIM_FOLDER=/u/jgao/jgao/projects/chrom_minigraph/sim

for j in {0..9}
do
    $CMG map -a --trace $GENERATED_FILE $SIM_FOLDER/$1/reads/ref${j}_len$1.fasta.fastq
done