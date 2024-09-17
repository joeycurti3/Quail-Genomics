#!/bin/bash
#$ -l h_data=4G,h_rt=6:00:00
#$ -wd <insert working directory>
#$ -o <insert log directory>
#$ -e <insert log directory>
#$ -m abe
#$ -N step01_a_CAQU_checkMD5sums

# Version: v1
# Usage: qsub step01_a_CAQU_checkMD5sums_20220927.sh
# Description: Check if file transfer of California quail genomes from QB3 was successfull using MD5
# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Date: WED SEP 27 2022

## Setup workspace 

source <insert path to miniconda>
conda activate CRLF
 
set -xeo pipefail

# working directories
WORKDIR= <insert directory with sequencing data>

## Main 

# Check current MD5 against what was given to my by sequencing facility
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}: CHECKING MD5 SUMS FOR DOWNLOADED DATA..."
cd ${WORKDIR}
md5sum --check md5sum.txt

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"

## cleanup 

# Add success to prog log
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done with MD5sums"
