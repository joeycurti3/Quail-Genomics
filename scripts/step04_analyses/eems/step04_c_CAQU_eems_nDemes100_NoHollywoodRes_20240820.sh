#!/bin/bash
#$ -l highp,h_data=10G,h_rt=24:00:00,h_vmem=20G
#$ -pe shared 2
#$ -wd <insert path to working directory>
#$ -o <insert path to log directory>
#$ -e <insert path to log directory>
#$ -m bae
#$ -M 1joeynik
#$ -N step04_c_CAQU_eems_nDemes100_NoHollywoodRes

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Version: v2 - removing Hollywood Res sample to see how that affects output 
# Usage: qsub step04_c_CAQU_eems_nDemes100_NoHollywoodRes_20240820.sh
# Date: TUE AUG 20 2024
# Description: Estimate migration surface using pairwise genetic distance matrix from bed2diffs using eems package
# References: 
# https://github.com/dipetkov/eems/blob/master/Documentation/EEMS-doc.pdf

## Import Packages
# Note: you will need to change to directory for eems package 

. /u/local/Modules/default/init/modules.sh
module load apptainer/1.2.2

set -o pipefail

## Define Variables

## Main

echo -e "[$(date "+%Y-%m-%d %T")] Starting eems, first MCMC run"

apptainer exec $H2_CONTAINER_LOC/h2-eems.sif runeems_snps --params params.chain1.ndemes100.NHR --seed 100

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done"

echo -e "[$(date "+%Y-%m-%d %T")] Second MCMC run"

apptainer exec $H2_CONTAINER_LOC/h2-eems.sif runeems_snps --params params.chain2.ndemes100.NHR --seed 200

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done"

echo -e "[$(date "+%Y-%m-%d %T")] Third MCMC run"

apptainer exec $H2_CONTAINER_LOC/h2-eems.sif runeems_snps --params params.chain3.ndemes100.NHR --seed 300

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done"

## Clean Up
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done running eems with nDemes = 100"
