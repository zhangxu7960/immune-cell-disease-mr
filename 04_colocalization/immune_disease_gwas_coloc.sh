#!/bin/bash
#SBATCH -J immune_disease_gwas_susie
#SBATCH -p cydell_2
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --array=1-165%4
#SBATCH -o %x_%a.log
#SBATCH -e %x_%a.log
#SBATCH --time=9999:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zx13132357960@163.com

libdir=/public/home/zhangxu02/immune_disease/colocalization/immune_disease_gwas_colocalization

echo "${SLURM_ARRAY_TASK_ID} start time: $(date)"

module load apps/R/4.3.1

Rscript ${libdir}/immune_disease_gwas_coloc.R ${SLURM_ARRAY_TASK_ID}

echo "${SLURM_ARRAY_TASK_ID} end time: $(date)"

date
mkdir -p $libdir/log
mv $(pwd)/${SLURM_JOB_NAME}_${SLURM_ARRAY_TASK_ID}.log $libdir/log/${SLURM_JOB_NAME}_${SLURM_ARRAY_TASK_ID}.log
