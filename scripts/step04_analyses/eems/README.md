### <ins>The following files/scripts are in this directory and should be executed in this order:</ins> 

1. step04_a_CAQU_VCF2Plink_NoHollywoodRes_20240820.sh
2. step04_b_CAQU_bed2diffs_NoHollywoodRes_20240820.sh
3. step04_c_CAQU_eems_nDemes100_NoHollywoodRes_20240820.sh
4. step04_c_CAQU_eems_nDemes200_NoHollywoodRes_20240820.sh
5. step04_c_CAQU_eems_nDemes300_NoHollywoodRes_20240820.sh
6. step04_d_CAQU_eemsplots2_20240809.R
7. params.chain1.ndemes100.NHR - EEMS parameters file for nDemes = 100 and for the first MCMC chain
8. params.chain2.ndemes100.NHR - EEMS parameters file for nDemes = 100 and for the second MCMC chain
9. params.chain3.ndemes100.NHR - EEMS parameters file for nDemes = 100 and for the third MCMC chain
10. params.chain1.ndemes200.NHR - EEMS parameters file for nDemes = 200 and for the first MCMC chain
11. params.chain2.ndemes200.NHR - EEMS parameters file for nDemes = 200 and for the second MCMC chain
12. params.chain3.ndemes200.NHR - EEMS parameters file for nDemes = 200 and for the third MCMC chain
13. params.chain1.ndemes300.NHR - EEMS parameters file for nDemes = 300 and for the first MCMC chain
14. params.chain2.ndemes300.NHR - EEMS parameters file for nDemes = 300 and for the second MCMC chain
15. params.chain3.ndemes300.NHR - EEMS parameters file for nDemes = 300 and for the third MCMC chain
16. samples_v2.txt - Sample names for input into Step 1
17. GCA_023055505.1_bCalCai1.0.p_PassSNPs_unrelated_SAMO_chr1.coord - Coordinates for polygon around area where EEMS analysis will run, needs to have same prefix as outputs of Step 1

Note: eems is run for 100, 200, and 300 nDemes, or the number of partitions on the landscape. Further, for each value of nDemes, I also run 3 replicate MCMC chains.

### <ins>Important note to run these scripts:</ins> 
Currently, bed2diffs cannot handle organisms with atypical chromosomal configurations (i.e., scaffold level assemblies).

This is addressed here:
https://github.com/dipetkov/eems/issues/15

To get around this, I set all scaffolds names in the .bim file to "1" and ran eems assuming just this one chromosome, as directed by the package author.

To do this, I ran:
cp datapath.bed newdatapath.bed
cp datapath.fam newdatapath.fam
awk '{print 1,$2,$3,$4,$5,$6}' datapath.bim  > newdatapath.bim

### <ins>Metadata</ins> 
#### *samples_v2.txt*
| Column Name | Data Type | Description | Comments |
| :---: | :---: | :---: | :---: | 
| Column 1 | numeric | Family ID (FID) field required by PLINK, set to zero because no information about family structure | NA |
| Column 2 | character | Sample name | NA |

#### *GCA_023055505.1_bCalCai1.0.p_PassSNPs_unrelated_SAMO_chr1.coord*
| Column Name | Data Type | Description | Comments |
| :---: | :---: | :---: | :---: | 
| Column 1 | numeric | Longitude of point in perimeter polygon | EPSG:4269: NAD83 |
| Column 2 | numeric | Latitude of point in perimeter polygon | EPSG:4269: NAD83 |

#### *params.chain[insert chain number].ndemes[insert demes number].NHR*
| Parameter | Data Type | Description | Example Value |
|-----------|-----------|-------------|---------------|
| datapath | String | Absolute path to the input EEMS dataset files, specified without their file extensions. | /u/home/1/1joeynik/project-rwayne/CAQU/scripts/step04_analyses/2024/eems/GCA_023055505.1_bCalCai1.0.p_PassSNPs_unrelated_SAMO_NoHollywoodRes_chr1 |
| mcmcpath | String | Absolute path to the directory where EEMS writes MCMC output files for the specified analysis and chain. | /u/home/1/1joeynik/project-rwayne/CAQU/scripts/step04_analyses/2024/eems/outputs_NoHollywoodRes/nIndiv56-nSites25199223-nDemes100-chain1 |
| nIndiv | Integer | Number of diploid individuals included in the analysis. | 55 |
| nSites | Integer | Number of genetic sites (SNPs) included in the analysis. | 25199223 |
| nDemes | Integer | Number of demes used to discretize the study area in the EEMS model. | 100 |
| diploid | Boolean | Specifies whether the input genetic data represent diploid individuals. | true |
| numMCMCIter | Integer | Number of MCMC iterations performed after burn-in. | 5000000 |
| numBurnIter | Integer | Number of initial MCMC iterations discarded as burn-in. | 2500000 |
| numThinIter | Integer | Interval at which MCMC samples are retained to reduce autocorrelation. | 9999 |
| qEffctProposalS2 | Numeric | Proposal variance for MCMC updates to effective diversity parameters (q). | 0.007 |
| mEffctProposalS2 | Numeric | Proposal variance for MCMC updates to effective migration parameters (m). | 2.7 |
| mrateMuProposalS2 | Numeric | Proposal variance for MCMC updates to the mean effective migration rate parameter. | 0.03 |
| mSeedsProposalS2 | Numeric | Proposal variance for MCMC updates to the spatial locations of migration-rate Voronoi seeds. | 0.03 |
