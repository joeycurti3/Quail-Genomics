#!/bin/bash
#$ -l highp,h_data=20G,h_rt=24:00:00
#$ -wd <insert working directory path>
#$ -o <insert log directory path>
#$ -e <insert log directory path>
#$ -m bae
#$ -M 1joeynik
#$ -N step04_a_CAQU_BCFTools_ROH

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Adapted From: Chris Kyriazis
# Version: v3 - Removing Hollywood Res sample
# Date: THU AUG 22 2024
# Usage: qsub step04_a_CAQU_BCFTools_ROH_20240822.sh
# Description: run BCFTools ROH output, takes VCF file as input

## Import Packages

source <insert path to miniconda>
conda activate my_bcftools

set -o pipefail

## Define Variables

INDIR=<insert path to VCF with merged scaffolds that contains both invariant and variant sites>
INFILE=GCA_023055505.1_bCalCai1.0.p_MergedScaffolds.vcf.gz
OUTDIR=<insert output directory path>
mkdir -p ${OUTDIR}
OUTFILE=GCA_023055505.1_bCalCai1.0.p_BCFTools_ROH_NoHollywoodRes.out

## Main

# for -G30 flag: from manual "use genotypes (FORMAT/GT fields) ignoring genotype likelihoods (FORMAT/PL), setting PL of unseen genotypes to FLOAT. Safe value to use is 30 to account for GT errors."

bcftools roh -G30 ${INDIR}/${INFILE} --include 'FILTER="PASS"' -S samples.txt -o ${OUTDIR}/${OUTFILE}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done with ROH"

## Clean Up
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done running ROH analysis in BCFTools"
conda deactivate
