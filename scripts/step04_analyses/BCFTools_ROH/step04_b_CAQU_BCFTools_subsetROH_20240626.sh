#!/bin/bash
#$ -l highp,h_data=20G,h_rt=24:00:00
#$ -wd <insert working directory path>
#$ -o <insert log directory path>
#$ -e <insert log directory path>
#$ -m bae
#$ -M 1joeynik
#$ -N step04_b_CAQU_BCFTools_subsetROH

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Adapted From: Chris Kyriazis
# Version: V2 - adding samples from MVZ, NHMLA, and WFVZ
# Date: WED JUN 26 2024
# Usage: qsub step04_b_CAQU_BCFTools_subsetROH_20240626.sh
# Description: subset output of BCFTools ROH analysis from previous script by individual

## Import Packages

set -o pipefail

## Define Variables

FILE=GCA_023055505.1_bCalCai1.0.p_BCFTools_ROH.out
INDIR=<insert output directory from previous script path>
SAMPLES=(
  T11B098 T1B075 T1B081 T1B083 T2B088 T3B066 T3B091 T3B092 T3B103 T4B035
  T4B073 T4B074 T6B067 T6B107 T7B036 T8B041 T8B070 T8B094 T1B050 T2B086
  T5B095 T4B006 T9B090 T4B110 T4B005 T1B082 T1B052 T4B009 T1B054 T2B085
  T05B038 T2B084 T2B029 T4B057 T11B111 T5B093 T2B001 T2B089 T4B096 T4B004
  T6B055 T2B087 T2B060 T3B024 T3B109 T3B032 T7B013 LACM107363 LACM107541
  LACM112287 WFVZ52698 WFVZ53206 MVZCCGP-CaOr35_I-C04 MVZCCGP-CaOr33_I-A04
  MVZCCGP-CaOr30_I-F03 MVZCCGP-CaOr104_I-H10 MVZCCGP-CaOr81_II-B02 MVZCCGP-CaOr78_II-H01
  MVZCCGP-CaOr71_II-F01 MVZCCGP-CaOr96_I-A10 MVZCCGP-CaOr45_I-C05 MVZCCGP-CaOr37_II-A01
)

## Main

echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Starting subset of BCFTools ROH output"

cd ${INDIR}

tail -n +6 ${FILE} > tmp.txt #this removes the header that is 6 lines long
 
for i in "${SAMPLES[@]}"; do
  awk '$2 == "'$i'" { print }' tmp.txt > ${FILE}_$i
done

wait

gzip *_${FILE}
rm tmp.txt

## Clean Up
echo "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Done with subsetting ROH file"
