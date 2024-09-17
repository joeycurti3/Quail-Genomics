#!/bin/bash
#$ -l highp,h_data=20G,h_rt=20:00:00
#$ -wd <insert working directory>
#$ -o <insert log directory>
#$ -e <insert log directory>
#$ -m abe

# Version: v1 
# Usage: qsub step01_c_CRLF_RefGenDict_GCA_023055505.1_bCalCai1.0.p_20220927.sh
# Description: making sequence dictionary for ref genome GCA_023055505.1_bCalCai1.0.p
# Author: Joey Curti (jcurti3@g.ucla.edu)
# Date: WED Sep 28 2022
# References: 
#https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk

## Setup workspace 

source <insert miniconda directory>
conda activate CRLF

set -xeo pipefail

## Define variables 

HOMEDIR=/u/home/1/1joeynik/project-rwayne/CAQU
WORKDIR=${HOMEDIR}/ref_genome/GCA_023055505.1_bCalCai1.0.p
REFERENCE_SEQ='GCA_023055505.1_bCalCai1.0.p_genomic.fasta' # should be name of reference seq, extension should be fasta but not fna

## Main 

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; Start making refseq dictionary for reference sequence ${REFERENCE_SEQ}"
cd ${WORKDIR}

# generate index for reference sequence from BWA

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; bwa index"

bwa index ${REFERENCE_SEQ} # creates many files such as .amb .ann .bwt .pac .rbwt .rpac .rsa .sa

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${LOG}
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; Done"

# Samtools faidx

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; Samtools index"

samtools faidx ${REFERENCE_SEQ} # Creates .fai file

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${LOG}
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; Done"

# generate dictionary

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; Picard sequence dictionary"

picard -Xmx4G CreateSequenceDictionary \
REFERENCE=${REFERENCE_SEQ} # creates .dict file

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" 
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; Done"

## Cleanup 

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done creating index and dictionary for reference sequence"
conda deactivate
