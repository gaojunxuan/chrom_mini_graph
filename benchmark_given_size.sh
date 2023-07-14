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
BAM_FILE=/u/jgao/jgao/code/chrom_mini_graph/report/$1/output.bam
SORTED_BAM=/u/jgao/jgao/code/chrom_mini_graph/report/$1/sorted_output.bam
PYTHON3=/u/jgao/miniconda3/bin/python3
SAMTOOLS=/.mounts/labs/simpsonlab/sw/samtools/1.17/bin/samtools
QSUB=/opt/uge-8.6/bin/lx-amd64/qsub
SIM_FOLDER=/u/jgao/jgao/projects/chrom_minigraph/sim

mkdir -p $CMG_FOLDER/report/
mkdir -p $CMG_FOLDER/report/$1
mkdir -p $CMG_FOLDER/report/$1/logs
mkdir -p $CMG_FOLDER/report/$1/logs/gen
mkdir -p $CMG_FOLDER/report/$1/logs/map

$CMG generate $SIM_FOLDER/$1/refs/*.fasta &> $CMG_FOLDER/report/$1/logs/gen/log_generate_size_${1}.txt
for j in {0..9}
do
    $CMG map -a --trace $GENERATED_FILE $SIM_FOLDER/$1/reads/ref${j}_len$1.fasta.fastq &> $CMG_FOLDER/report/$1/logs/map/log_map_size_${1}_ref_${j}.txt
    $SAMTOOLS sort $BAM_FILE -o $SORTED_BAM
    $SAMTOOLS coverage $SORTED_BAM > $CMG_FOLDER/report/$1/coverage_report_size_$1_ref_${j}.txt
    $PYTHON3 $CMG_FOLDER/analyze_mapping_accuracy.py ref${j}_len$1 $CMG_FOLDER/best_genome_reads.txt > $CMG_FOLDER/report/$1/accuracy_report_size_${1}_ref_${j}.txt
done