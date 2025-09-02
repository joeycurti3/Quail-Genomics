#!/bin/bash
#$ -l highp,h_data=10G,h_rt=24:00:00,h_vmem=20G
#$ -pe shared 2
#$ -wd /u/home/1/1joeynik/project-rwayne/CAQU/scripts/step04_analyses/2024/Plink_Fst/2025
#$ -o /u/home/1/1joeynik/project-rwayne/CAQU/logs/step04_c_CAQU_Plink_Fst_unrelated_20250901.out.txt
#$ -e /u/home/1/1joeynik/project-rwayne/CAQU/logs/step04_c_CAQU_Plink_Fst_unrelated_20250901.err.txt
#$ -m bae
#$ -M 1joeynik
#$ -N step04_c_CAQU_Plink_Fst_unrelated

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Version: v3 - rerunning with maf set to 0.012 
# Usage: qsub step04_c_CAQU_Plink_Fst_unrelated_20250901.sh
# Description: Commpute pair-wise Fst values using Hudson estimator via Plink 2.0

## Import Packages

source /u/home/1/1joeynik/project-rwayne/software/miniconda3/etc/profile.d/conda.sh
conda activate my_plink2

set -o pipefail

## Define Variables

WORKDIR=/u/home/1/1joeynik/project-rwayne/CAQU/scripts/step04_analyses/2024/Plink_Fst/2025
INFILE="/u/home/1/1joeynik/project-rwayne/CAQU/preprocessing/VCFs/2024/PassSNPs/GCA_023055505.1_bCalCai1.0.p_PassSNPs.vcf.gz"
POPFILE=${WORKDIR}/popid_unrelated.txt
OUTFILE2=GCA_023055505.1_bCalCai1.0.p_PassSNPs_Plink_Hudson_unrelated_20250901

## Main

# Plink Hudson Estimator

echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Starting Hudson estimator"

plink2 \
        --vcf ${INFILE} \
        --allow-extra-chr \
        --const-fid 0 \
        --within ${POPFILE} \
        --out ${OUTFILE2} \
        --fst CATPHENO method=hudson blocksize=200

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done"

## Clean Up
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done running pairwise F_st in Plink2"
conda deactivate
