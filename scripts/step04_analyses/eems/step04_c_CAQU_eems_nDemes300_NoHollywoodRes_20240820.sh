#!/bin/bash
#$ -l highp,h_data=10G,h_rt=24:00:00,h_vmem=20G
#$ -pe shared 2
#$ -wd <insert working directory path>
#$ -o <insert log directory path>
#$ -e <insert log directory path>
#$ -m bae
#$ -M 1joeynik
#$ -N step04_c_CAQU_eems_nDemes300_NoHollywoodRes

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Version: v2 - Removing Hollywood Res samples to see how it affects output
# Usage: qsub step04_c_CAQU_eems_nDemes300_NoHollywoodRes_20240820.sh
# Date: TUE AUG 20 2024
# Description: Estimate migration surface using pairwise genetic distance matrix from bed2diffs using eems package
# References: 
# https://github.com/dipetkov/eems/blob/master/Documentation/EEMS-doc.pdf

## Import Packages
# Note: change to path to eems package

. /u/local/Modules/default/init/modules.sh
module load apptainer/1.2.2

set -o pipefail

## Define Variables

## Main

echo -e "[$(date "+%Y-%m-%d %T")] Starting eems, first MCMC run"

apptainer exec $H2_CONTAINER_LOC/h2-eems.sif runeems_snps --params params.chain1.ndemes300.NHR --seed 100

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done"

echo -e "[$(date "+%Y-%m-%d %T")] Second MCMC run"

apptainer exec $H2_CONTAINER_LOC/h2-eems.sif runeems_snps --params params.chain2.ndemes300.NHR --seed 200

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi
echo -e "[$(date "+%Y-%m-%d %T")] Done"

echo -e "[$(date "+%Y-%m-%d %T")] Third MCMC run"

apptainer exec $H2_CONTAINER_LOC/h2-eems.sif runeems_snps --params params.chain3.ndemes300.NHR --seed 300

exitVal=${?}#if [ ${exitVal} -ne 0 ]; then
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done"

## Clean Up
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done running eems with nDemes = 300"
