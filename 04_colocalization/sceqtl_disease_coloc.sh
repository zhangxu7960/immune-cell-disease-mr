#!/bin/bash
#SBATCH -J sceqtl_disease_susie
#SBATCH -p cydell_2
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --array=1-126%8
#SBATCH -o %x_%a.log
#SBATCH -e %x_%a.log
#SBATCH --time=9999:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zx13132357960@163.com

libdir=/public/home/zhangxu02/immune_disease/colocalization/sceQTL_disease_colocalization

echo "${SLURM_ARRAY_TASK_ID} start time: $(date)"

module load apps/R/4.3.1

Rscript ${libdir}/sceqtl_disease_coloc.R ${SLURM_ARRAY_TASK_ID}

echo "${SLURM_ARRAY_TASK_ID} end time: $(date)"

date
mkdir -p $libdir/log
mv $(pwd)/${SLURM_JOB_NAME}_${SLURM_ARRAY_TASK_ID}.log $libdir/log/${SLURM_JOB_NAME}_${SLURM_ARRAY_TASK_ID}.log
