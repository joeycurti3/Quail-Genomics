#! /bin/bash
#$ -l highp,h_rt=120:00:00,h_data=20G,h_vmem=40G
#$ -pe shared 2 
#$ -wd <insert path to working directory>
#$ -o <insert path to log directory>
#$ -e <insert path to log directory>
#$ -m abe
#$ -N step03_d_CAQU_FilterTally
#$ -t <adjust according to how many scaffolds you can run simultaneously>

# Version: V2 - adding in samples from MVZ, NHMLA, and WFVZ
# Usage: qsub step03_d_CAQU_FilterTally_20240612.sh
# Description: Run BCFTools filter tally function and stats  on filtered VCFs to determine number of SNPs and what was filtered out
# Adapted By: Joseph Curti (jcurti3@g.ucla.edu)
# Author: Meixi Lin 
# Date: WED JUN 12 2024

## SETUP WORKSPACE  ##

sleep $((RANDOM % 120))
source <insert path to miniconda>
conda activate my_bcftools

set -o pipefail

## Define Variables

WORKDIR=<insert path to output directory from custom python script filtering step>
OUTDIR=<insert output directory path>
mkdir -p ${OUTDIR}

# Define Functions ##

tally_filters() {
    bcftools query -f '%FILTER\n' ${1} | sort --parallel=6 | uniq -c | sort -nr > ${2}
}

## Main

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; Starting BCFTools filter tally and post-filter stats"

cd ${WORKDIR}

# BCFTools filter tally

echo -e "[$(date "+%Y-%m-%d %T")] Generating site-level filter summary for INT${SGE_TASK_ID} ..."

tally_filters ${WORKDIR}/${REF}_INT${SGE_TASK_ID}_CustomPyFilter.vcf.gz ${OUTDIR}/${REF}_INT${SGE_TASK_ID}_CustomPyFilter_SiteFilterTally.txt

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done with BCFTools filter tally"

# BCFTools stats

echo -e "[$(date "+%Y-%m-%d %T")] Generate BCFTOOLS stats for INT${SGE_TASK_ID} ... "

bcftools stats -f 'PASS' -s- --threads 3 ${WORKDIR}/${REF}_INT${SGE_TASK_ID}_CustomPyFilter.vcf.gz > ${OUTDIR}/${REF}_INT${SGE_TASK_ID}_CustomPyFilter_BCFstatsPass.txt 

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done with BCFTools stats"

## Cleanup

echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done running BCFTools filter tally and stats"

conda deactivate

