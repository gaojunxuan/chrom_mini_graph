#!/bin/bash
#$ -P simpsonlab
#$ -N generate_sim_reads
#$ -cwd
#$ -l h_vmem=32G
#$ -l h_rt=144000
#$ -m eas
#$ -M jgao@oicr.on.ca

BADREAD=/u/jgao/miniconda3/bin/badread
SIZES=(5000000)
for d in ${SIZES[@]}; do
    dir=/u/jgao/jgao/code/rand_seq/$d
    for f in ${dir}/refs/* ; do
        echo $f
        fname="$(basename $f)"
        qsub -N run_badread_${fname} /u/jgao/jgao/code/rand_seq/run_badread.sh $f ${dir}
    done
done

SIZES=(100000 250000 500000 750000 1000000 2500000 5000000 7500000 10000000 25000000)
for i in ${SIZES[@]}; do
    qsub -N benchmark_minigraph_map_$i ./benchmark_minigraph.sh $i
done

SIZES=(100000 250000 500000 750000 1000000 2500000 5000000 7500000 10000000)
for i in ${SIZES[@]}; do
    qsub -N benchmark_cmg_map_noaln_$i run_map.sh $i
done

SIZES=(100000 250000 500000 750000 1000000 2500000 5000000 7500000 10000000 25000000)
for i in ${SIZES[@]}; do
    qsub -N benchmark_cmg_gen_fast_$i run_generate.sh $i
done