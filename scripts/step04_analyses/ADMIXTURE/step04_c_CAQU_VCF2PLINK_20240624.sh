#!/bin/bash
#$ -l highp,h_data=20G,h_rt=24:00:00
#$ -wd <insert working directory path>
#$ -o <insert log directory path>
#$ -e <insert log directory path>
#$ -m bae
#$ -M 1joeynik
#$ -N step04_c_CAQU_VCF2PLINK

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Version: v2 - adding samples from MVZ, NHMLA, and WFVZ
# Usage: qsub step04_c_CAQU_VCF2PLINK_20240624.sh
# Description: Convert VCF to PLINK .PED file
# Date: MON JUN 24 2024

## Import Packages

source <insert path to miniconda>
conda activate my_SNPRelate

set -o pipefail

## Define Variables

infile_allsamples="GCA_023055505.1_bCalCai1.0.p_PassSNPs_maf01_ld2_allsamples_GDS2VCF.vcf"
outfile_allsamples="GCA_023055505.1_bCalCai1.0.p_PassSNPs_maf01_ld2_allsamples_PLINK"

infile_SAMO="GCA_023055505.1_bCalCai1.0.p_PassSNPs_maf01_ld2_SAMO_GDS2VCF.vcf"
outfile_SAMO="GCA_023055505.1_bCalCai1.0.p_PassSNPs_maf01_ld2_SAMO_PLINK"

## Main

# Convert VCF to PLINK file for all unrelated samples including outgroups

plink --vcf ${infile_allsamples} --const-fid 0 --allow-extra-chr --recode 12 --out ${outfile_allsamples}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done with VCF2PLINK for all unrelated samples"

# Convert VCF to PLINK file for all unrelated samples in SAMO

plink --vcf ${infile_SAMO} --const-fid 0 --allow-extra-chr --recode 12 --out ${outfile_SAMO}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done with VCF2PLINK for all unrelated SAMO samples"

## Clean Up
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done converting VCF to PLINK file"
conda deactivate

