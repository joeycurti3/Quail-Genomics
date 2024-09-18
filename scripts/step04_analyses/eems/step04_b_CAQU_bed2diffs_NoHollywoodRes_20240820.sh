#!/bin/bash
#$ -l highp,h_data=15G,h_rt=24:00:00,h_vmem=30G
#$ -pe shared 2
#$ -wd <insert path to working directory>
#$ -o <insert path to log directory>
#$ -e <insert path to log directory>
#$ -m bae
#$ -M 1joeynik
#$ -N step04_b_CAQU_bed2diffs_NoHollywoodRes

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Version: v2 - removing hollywood res sample to see how it affects outputs
# Usage: qsub step04_b_CAQU_bed2diffs_NoHollywoodRes_20240820.sh
# Date: TUE AUG 20 2024
# Description: Convert VCF in .bed format to a pairwise genetic distance matrix using bed2diffs
# References: 
# https://github.com/dipetkov/eems/blob/master/Documentation/EEMS-doc.pdf

## Import Packages
# Note: I had an apptainer built with eems installed in it including dependencies. You will need to change this to your install of eems.

. /u/local/Modules/default/init/modules.sh
module load apptainer/1.2.2

set -o pipefail

## Define Variables

INFILE=GCA_023055505.1_bCalCai1.0.p_PassSNPs_unrelated_SAMO_NoHollywoodRes_chr1

## Main

echo -e "[$(date "+%Y-%m-%d %T")] Starting bed2diffs"

apptainer exec $H2_CONTAINER_LOC/h2-eems.sif bed2diffs_v1 --bfile ${INFILE} --nthreads 2

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done"

## Clean Up
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done running pairwise F_st in VCFTools"
