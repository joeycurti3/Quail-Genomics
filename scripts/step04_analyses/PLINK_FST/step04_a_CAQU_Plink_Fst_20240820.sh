#!/bin/bash
#$ -l highp,h_data=10G,h_rt=24:00:00,h_vmem=20G
#$ -pe shared 2
#$ -wd <insert workding directory path>
#$ -o <insert log directory path>
#$ -e <insert log directory path>
#$ -m bae
#$ -M 1joeynik
#$ -N step04_a_CAQU_Plink_Fst

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Version: v2 - adding in sample from NE topanga
# Usage: qsub step04_a_CAQU_Plink_Fst_20240820.sh
# Description: Commpute pair-wise Fst values using Weir and Cockerham and Hudson estimators via Plink 2.0

## Import Packages

source <insert miniconda path>
conda activate my_plink2

set -o pipefail

## Define Variables

WORKDIR=<insert working directory path>
INFILE="/u/home/1/1joeynik/project-rwayne/CAQU/preprocessing/VCFs/2024/PassSNPs/GCA_023055505.1_bCalCai1.0.p_PassSNPs.vcf.gz"
POPFILE=<insert path for file 'popid.txt'>
OUTFILE1=GCA_023055505.1_bCalCai1.0.p_PassSNPs_Plink_WC_20240820
OUTFILE2=GCA_023055505.1_bCalCai1.0.p_PassSNPs_Plink_Hudson_20240820

## Main

# Plink WC Estimator

echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Starting WC estimator"

plink2 \
	--vcf ${INFILE} \
	--allow-extra-chr \
	--const-fid 0 \
	--within ${POPFILE} \
	--out ${OUTFILE1} \
        --fst CATPHENO method=wc


exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done"

# Plink Hudson Estimator

echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Starting Hudson estimator"

plink2 \
        --vcf ${INFILE} \
        --allow-extra-chr \
        --const-fid 0 \
        --within ${POPFILE} \
        --out ${OUTFILE2} \
        --fst CATPHENO method=hudson

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done"

## Clean Up
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done running pairwise F_st in Plink2"
conda deactivate
