#!/bin/bash
#$ -l highp,h_data=10G,h_rt=2:00:00
#$ -wd <insert path to working directory>
#$ -o <insert path to log directory>
#$ -e <insert path to log directory>
#$ -m bae
#$ -M 1joeynik
#$ -N step04_a_CAQU_VCF2Plink

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Adapted from: Chris Kyriazis (https://github.com/ckyriazis/)
# Version: v2 - removing hollywood reservoir sample to see how it affects results
# Usage: qsub step04_a_CAQU_VCF2Plink_NoHollywoodRes_20240820.sh
# Description: Convert filtered VCF file with only biallelic 'PASS' sites into a PLINK .bed file for use with eems

## Import Packages
# Note: You will need to change this to path to modules on your computing cluster!
. /u/local/Modules/default/init/modules.sh
module load plink/1.90b624

set -o pipefail

## Define Variables

INFILE_ALL=<insert path to VCF with all biallelic passing SNPs>
WORKDIR=<insert path to working directory>
KEEP_FILE=${WORKDIR}/samples_v2.txt # I needed to filter out a sample with high missingness, so this corresponds to a list with sample names of all samples to keep

## Main

# From Chris:
# you need to use const-fid 0 otherwise it thinks that family name_sample name is structure of ID and tries to split it (and fails)
# allow extra chromosomes: to get it to get over the fact that chr names are non standard
# for faststructure to work you have to filter on maf 0.05

echo -e "[$(date "+%Y-%m-%d %T")] Starting PLINK VCF to BED conversion on all unrelated samples"

plink --vcf ${INFILE_ALL} --keep ${KEEP_FILE} --make-bed --keep-allele-order --const-fid 0 --allow-extra-chr --out GCA_023055505.1_bCalCai1.0.p_PassSNPs_unrelated_SAMO_NoHollywoodRes

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" 
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done with all unrelated CAQU samples"

## Clean up

echo -e "[$(date "+%Y-%m-%d %T")] Done running PLINK VCF2BED"

