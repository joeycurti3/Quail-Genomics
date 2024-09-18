#!/bin/bash
#$ -l highp,h_data=20G,h_rt=12:00:00
#$ -wd <insert path to working directory>
#$ -o <insert log directory path>
#$ -e <insert log directory path>
#$ -m bae
#$ -M 1joeynik
#$ -N step04_d_CAQU_ADMIXTURE_SAMO

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Adapted from: Meixi Lin (meixilin@ucla.edu)
# Version: v2 - including samples from MVZ, NHMLA, and WFVZ
# Usage: qsub step04_d_CAQU_ADMIXTURE_SAMO_20240624.sh
# Description: With input PLINK BED file that has been LD Pruned (r^2 = .2),MAF=.05,and highly related individuals removedl;  run ADMIXTURE with K=1:10 and with 10 iterations
# Reference: http://dalexander.github.io/admixture/admixture-manual.pdf
# Date: MON JUN 24 2024

## Import Packages

source <insert path to miniconda>
conda activate my_admixture

set -o pipefail

## Define Variables

infile1="GCA_023055505.1_bCalCai1.0.p_PassSNPs_maf01_ld2_SAMO_PLINK.ped"
outfile1="GCA_023055505.1_bCalCai1.0.p_PassSNPs_maf01_ld2_SAMO_PLINK"
outfile2="GCA_023055505.1_bCalCai1.0.p_PassSNPs_maf01_ld2_SAMO_ADMIXTURE_20240619"

## Main

for K in {1..10};do
        for i in {1..10};do
            admixture --cv -s time -j8 ${infile1} ${K} | tee log_K${K}.iter${i}_SAMO.out
            mv ${outfile1}.${K}.Q ${outfile2}.K${K}.iter${i}_SAMO.Q
            mv ${outfile1}.${K}.P ${outfile2}.K${K}.iter${i}_SAMO.P
            # get the CV error and loglikelihood during each run
            CVERROR=$(awk '/^CV/ {print $4}' log_K${K}.iter${i}_SAMO.out)
            LL=$(awk '/^Loglikelihood/ {print $2}' log_K${K}.iter${i}_SAMO.out)
            echo -e "${K},${i},${CVERROR},${LL}" >> Admixture_CV_LLsummary_SAMO_20240624.csv
        done
    done

# move outputs to their own directory, seperate out the Q matrices

mkdir -p <insert output directory path>
mkdir -p <insert subdirectory for Q files path>

mv *ADMIXTURE_2* *log_K* <insert output directory path>
cd <insert output directory path>
mv *.Q <insert Q file directory path>

## Cleanup
echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}.${SGE_TASK_ID}; Finished running admixture for K=1:10 and 10 iterations"
conda deactivate
