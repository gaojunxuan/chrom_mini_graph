#!/bin/bash
#$ -P simpsonlab
#$ -N chrom_minigraph_benchmark
#$ -cwd
#$ -l h_vmem=64G
#$ -l h_rt=144000
#$ -m eas
#$ -M jgao@oicr.on.ca

CMG_FOLDER=/u/jgao/jgao/code/chrom_mini_graph
CMG=/u/jgao/jgao/code/chrom_mini_graph/target/debug/chrom_mini_graph
GENERATED_FILE=/u/jgao/jgao/code/chrom_mini_graph/serialized_mini_graph.bin
BAM_FILE=/u/jgao/jgao/code/chrom_mini_graph/output.bam
SORTED_BAM=/u/jgao/jgao/code/chrom_mini_graph/sorted_output.bam
#SIZES=(100000 250000 500000 750000 1000000 2500000 5000000 7500000 10000000)
SIZES=(100000 1000000)
PYTHON3=/u/jgao/miniconda3/bin/python3
SAMTOOLS=/.mounts/labs/simpsonlab/sw/samtools/1.17/bin/samtools
QSUB=/opt/uge-8.6/bin/lx-amd64/qsub
SIM_FOLDER=/u/jgao/jgao/projects/chrom_minigraph/sim

for i in ${SIZES[@]}
do
    $QSUB -N benchmark_size_$i $CMG_FOLDER/benchmark_given_size.sh $i
done