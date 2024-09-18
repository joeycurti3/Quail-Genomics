#!/bin/bash
#$ -l highp,h_data=20G,h_rt=24:00:00
#$ -wd <insert working directory path>
#$ -o <insert log directory path>
#$ -e <insert log directory path>
#$ -m bae
#$ -M 1joeynik
#$ -N step04_c_CAQU_IndGeneticDist

## Import Packages

source <insert path to miniconda>
conda activate my_geneticdistance

set -o pipefail

## Main

Rscript step04_c_CAQU_IndGeneticDist_20240820.R

## Clean Up
echo -e "[$(date "+%Y-%m-%d %T")] Done"

conda deactivate

