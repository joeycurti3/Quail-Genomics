#!/bin/bash
#$ -l highp,h_rt=48:00:00,h_data=6G,h_vmem=48G
#$ -pe shared 8
#$ -wd <insert working directory path>
#$ -o <insert log directory path>
#$ -e <insert log directory path>
#$ -m bae
#$ -M 1joeynik
#$ -N step04_a_CAQU_SlidingWindowHet

# Adapted From: Chris Kyriazis and Jacqueline Robinson
# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Version: V2 - Adding in samples from MVZ, NHMLA, and WFVZ
# Usage: qsub -t 1-<number of scaffolds> run_step04_a_CAQU_SlidingWindowHet_20240623.sh
# Description: Runs custom sliding window heterozygosity script 
# Note: scaffolds that we excluded from this analysis because they didn't contain many SNPs have their values set to 'NULL' below, otherwise the script errors out.

## Import Packages

source <insert path to miniconda>
conda activate my_SlidingWindowHet

## Define variables

INDIR=<insert path to VCF that has been filtered but still contains both variant and invariant sites>
OUTDIR=<insert output directory path>
VCF=${INDIR}/GCA_023055505.1_bCalCai1.0.p_INT${SGE_TASK_ID}_CustomPyFilter.vcf.gz
SCRIPT=${OUTDIR}/step04_a_CAQU_SlidingWindowHet_20240623.py

if [ $SGE_TASK_ID == 1 ]
then
        SCAFFOLD=JALIRH010000001.1
elif [ $SGE_TASK_ID == 2 ]
then
	SCAFFOLD=JALIRH010000002.1
elif [ $SGE_TASK_ID == 3 ]
then
        SCAFFOLD=JALIRH010000003.1
elif [ $SGE_TASK_ID == 4 ]
then
        SCAFFOLD=JALIRH010000004.1
elif [ $SGE_TASK_ID == 5 ]
then
        SCAFFOLD=JALIRH010000005.1
elif [ $SGE_TASK_ID == 6 ]
then
        SCAFFOLD=JALIRH010000006.1
elif [ $SGE_TASK_ID == 7 ]
then
        SCAFFOLD=JALIRH010000007.1
elif [ $SGE_TASK_ID == 8 ]
then
        SCAFFOLD=JALIRH010000008.1
elif [ $SGE_TASK_ID == 9 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 10 ]
then
        SCAFFOLD=JALIRH010000010.1
elif [ $SGE_TASK_ID == 11 ]
then
        SCAFFOLD=JALIRH010000011.1
elif [ $SGE_TASK_ID == 12 ]
then
        SCAFFOLD=JALIRH010000012.1
elif [ $SGE_TASK_ID == 13 ]
then
        SCAFFOLD=JALIRH010000013.1
elif [ $SGE_TASK_ID == 14 ]
then
        SCAFFOLD=JALIRH010000014.1
elif [ $SGE_TASK_ID == 15 ]
then
        SCAFFOLD=JALIRH010000015.1
elif [ $SGE_TASK_ID == 16 ]
then
        SCAFFOLD=JALIRH010000016.1
elif [ $SGE_TASK_ID == 17 ]
then
        SCAFFOLD=JALIRH010000017.1
elif [ $SGE_TASK_ID == 18 ]
then
        SCAFFOLD=JALIRH010000018.1
elif [ $SGE_TASK_ID == 19 ]
then
        SCAFFOLD=JALIRH010000019.1
elif [ $SGE_TASK_ID == 20 ]
then
        SCAFFOLD=JALIRH010000020.1
elif [ $SGE_TASK_ID == 21 ]
then
        SCAFFOLD=JALIRH010000021.1
elif [ $SGE_TASK_ID == 22 ]
then
        SCAFFOLD=JALIRH010000022.1
elif [ $SGE_TASK_ID == 23 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 24 ]
then
        SCAFFOLD=JALIRH010000024.1
elif [ $SGE_TASK_ID == 25 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 26 ]
then
        SCAFFOLD=JALIRH010000026.1
elif [ $SGE_TASK_ID == 27 ]
then
        SCAFFOLD=JALIRH010000027.1
elif [ $SGE_TASK_ID == 28 ]
then
        SCAFFOLD=JALIRH010000028.1
elif [ $SGE_TASK_ID == 29 ]
then
        SCAFFOLD=JALIRH010000029.1
elif [ $SGE_TASK_ID == 30 ]
then
        SCAFFOLD=JALIRH010000030.1
elif [ $SGE_TASK_ID == 31 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 32 ]
then
        SCAFFOLD=JALIRH010000032.1
elif [ $SGE_TASK_ID == 33 ]
then
        SCAFFOLD=JALIRH010000033.1
elif [ $SGE_TASK_ID == 34 ]
then
        SCAFFOLD=JALIRH010000034.1
elif [ $SGE_TASK_ID == 35 ]
then
        SCAFFOLD=JALIRH010000035.1
elif [ $SGE_TASK_ID == 36 ]
then
        SCAFFOLD=JALIRH010000036.1
elif [ $SGE_TASK_ID == 37 ]
then
        SCAFFOLD=JALIRH010000037.1
elif [ $SGE_TASK_ID == 38 ]
then
        SCAFFOLD=JALIRH010000038.1
elif [ $SGE_TASK_ID == 39 ]
then
        SCAFFOLD=JALIRH010000039.1
elif [ $SGE_TASK_ID == 40 ]
then
        SCAFFOLD=JALIRH010000040.1
elif [ $SGE_TASK_ID == 41 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 42 ]
then
        SCAFFOLD=JALIRH010000042.1
elif [ $SGE_TASK_ID == 43 ]
then
        SCAFFOLD=JALIRH010000043.1
elif [ $SGE_TASK_ID == 44 ]
then
        SCAFFOLD=JALIRH010000044.1
elif [ $SGE_TASK_ID == 45 ]
then
        SCAFFOLD=JALIRH010000045.1
elif [ $SGE_TASK_ID == 46 ]
then
        SCAFFOLD=JALIRH010000046.1
elif [ $SGE_TASK_ID == 47 ]
then
        SCAFFOLD=JALIRH010000047.1
elif [ $SGE_TASK_ID == 48 ]
then
        SCAFFOLD=JALIRH010000048.1
elif [ $SGE_TASK_ID == 49 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 50 ]
then
        SCAFFOLD=JALIRH010000050.1
elif [ $SGE_TASK_ID == 51 ]
then
        SCAFFOLD=JALIRH010000051.1
elif [ $SGE_TASK_ID == 52 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 53 ]
then
        SCAFFOLD=JALIRH010000053.1
elif [ $SGE_TASK_ID == 54 ]
then
        SCAFFOLD=JALIRH010000054.1
elif [ $SGE_TASK_ID == 55 ]
then
        SCAFFOLD=JALIRH010000055.1
elif [ $SGE_TASK_ID == 56 ]
then
        SCAFFOLD=JALIRH010000056.1
elif [ $SGE_TASK_ID == 57 ]
then
        SCAFFOLD=JALIRH010000057.1
elif [ $SGE_TASK_ID == 58 ]
then
        SCAFFOLD=JALIRH010000058.1
elif [ $SGE_TASK_ID == 59 ]
then
        SCAFFOLD=JALIRH010000059.1
elif [ $SGE_TASK_ID == 60 ]
then
        SCAFFOLD=JALIRH010000060.1
elif [ $SGE_TASK_ID == 61 ]
then
        SCAFFOLD=JALIRH010000061.1
elif [ $SGE_TASK_ID == 62 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 63 ]
then
        SCAFFOLD=JALIRH010000063.1
elif [ $SGE_TASK_ID == 64 ]
then
        SCAFFOLD=JALIRH010000064.1
elif [ $SGE_TASK_ID == 65 ]
then
        SCAFFOLD=JALIRH010000065.1
elif [ $SGE_TASK_ID == 66 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 67 ]
then
        SCAFFOLD=JALIRH010000067.1
elif [ $SGE_TASK_ID == 68 ]
then
        SCAFFOLD=JALIRH010000068.1
elif [ $SGE_TASK_ID == 69 ]
then
        SCAFFOLD=JALIRH010000069.1
elif [ $SGE_TASK_ID == 70 ]
then
        SCAFFOLD=JALIRH010000070.1
elif [ $SGE_TASK_ID == 71 ]
then
        SCAFFOLD=JALIRH010000071.1
elif [ $SGE_TASK_ID == 72 ]
then
        SCAFFOLD=JALIRH010000072.1
elif [ $SGE_TASK_ID == 73 ]
then
        SCAFFOLD=JALIRH010000073.1
elif [ $SGE_TASK_ID == 74 ]
then
        SCAFFOLD=JALIRH010000074.1
elif [ $SGE_TASK_ID == 75 ]
then
        SCAFFOLD=JALIRH010000075.1
elif [ $SGE_TASK_ID == 76 ]
then
        SCAFFOLD=JALIRH010000076.1
elif [ $SGE_TASK_ID == 77 ]
then
        SCAFFOLD=JALIRH010000077.1
elif [ $SGE_TASK_ID == 78 ]
then
        SCAFFOLD=JALIRH010000078.1
elif [ $SGE_TASK_ID == 79 ]
then
        SCAFFOLD=JALIRH010000079.1
elif [ $SGE_TASK_ID == 80 ]
then
        SCAFFOLD=JALIRH010000080.1
elif [ $SGE_TASK_ID == 81 ]
then
        SCAFFOLD=JALIRH010000081.1
elif [ $SGE_TASK_ID == 82 ]
then
        SCAFFOLD=JALIRH010000082.1
elif [ $SGE_TASK_ID == 83 ]
then
        SCAFFOLD=JALIRH010000083.1
elif [ $SGE_TASK_ID == 84 ]
then
        SCAFFOLD=JALIRH010000084.1
elif [ $SGE_TASK_ID == 85 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 86 ]
then
        SCAFFOLD=JALIRH010000086.1
elif [ $SGE_TASK_ID == 87 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 88 ]
then
        SCAFFOLD=JALIRH010000088.1
elif [ $SGE_TASK_ID == 89 ]
then
        SCAFFOLD=JALIRH010000089.1
elif [ $SGE_TASK_ID == 90 ]
then
        SCAFFOLD=JALIRH010000090.1
elif [ $SGE_TASK_ID == 91 ]
then
        SCAFFOLD=JALIRH010000091.1
elif [ $SGE_TASK_ID == 92 ]
then
        SCAFFOLD=JALIRH010000092.1
elif [ $SGE_TASK_ID == 93 ]
then
        SCAFFOLD=JALIRH010000093.1
elif [ $SGE_TASK_ID == 94 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 95 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 96 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 97 ]
then
        SCAFFOLD=JALIRH010000097.1
elif [ $SGE_TASK_ID == 98 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 99 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 100 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 101 ]
then
        SCAFFOLD=JALIRH010000101.1
elif [ $SGE_TASK_ID == 102 ]
then
        SCAFFOLD=JALIRH010000102.1
elif [ $SGE_TASK_ID == 103 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 104 ]
then
        SCAFFOLD=JALIRH010000104.1
elif [ $SGE_TASK_ID == 105 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 106 ]
then
        SCAFFOLD=JALIRH010000106.1
elif [ $SGE_TASK_ID == 107 ]
then
        SCAFFOLD=NULL
elif [ $SGE_TASK_ID == 108 ]
then
        SCAFFOLD=JALIRH010000108.1
elif [ $SGE_TASK_ID == 109 ]
then
        SCAFFOLD=JALIRH010000109.1
elif [ $SGE_TASK_ID == 110 ]
then
        SCAFFOLD=JALIRH010000110.1
elif [ $SGE_TASK_ID == 111 ]
then
        SCAFFOLD=JALIRH010000111.1
elif [ $SGE_TASK_ID == 112 ]
then
        SCAFFOLD=JALIRH010000112.1
elif [ $SGE_TASK_ID == 113 ]
then
        SCAFFOLD=JALIRH010000113.1
elif [ $SGE_TASK_ID == 114 ]
then
        SCAFFOLD=JALIRH010000114.1
elif [ $SGE_TASK_ID == 115 ]
then
        SCAFFOLD=JALIRH010000115.1
elif [ $SGE_TASK_ID == 116 ]
then
        SCAFFOLD=JALIRH010000116.1
elif [ $SGE_TASK_ID == 117 ]
then
        SCAFFOLD=JALIRH010000117.1
elif [ $SGE_TASK_ID == 118 ]
then
        SCAFFOLD=JALIRH010000118.1
elif [ $SGE_TASK_ID == 119 ]
then
        SCAFFOLD=JALIRH010000119.1
fi

## Main

echo -e "[$(date "+%Y-%m-%d %T")] Starting heterozygosity estimation for INT${SGE_TASK_ID}"

cd ${INDIR}
python ${SCRIPT} ${VCF} 1000000 1000000 ${SCAFFOLD}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" 
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] Done with heterozygosity estimation for INT${SGE_TASK_ID}"

## Cleanup 

echo -e "[$(date "+%Y-%m-%d %T")] Job ID: ${JOB_ID}; Finished running Sliding Window Heterozygosity script"

