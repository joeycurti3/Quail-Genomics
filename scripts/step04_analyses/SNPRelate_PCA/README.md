### <ins>The following files/scripts are in this directory and should be executed in this order:</ins> 

1. run_step04_a_CAQU_SNPRelatePCA_ld2_maf05_unrelated_20240802.sh - Wrapper script needed to submit Step 2
2. step04_a_CAQU_SNPRelatePCA_ld2_maf05_unrelated_20240802.R - Prune input GDS file for linkage disequilibrium and minor allele frequency, then perform PCA on all samples
3. run_step04_a_CAQU_SNPRelatePCA_ld2_maf05_unrelated_SAMO_20240822.sh - Wrapper script needed to submit Step 4
4. step04_a_CAQU_SNPRelatePCA_ld2_maf05_unrelated_SAMO_20240822.R - Prune input GDS file for linkage disequilibrium and minor allele frequency, then perform PCA on subset of samples from Santa Monica Mountains
5. SNPRelate_PopID_unrelated_20240617.csv - Text file of sample names and the population name where the sample came from, input to Step 2
6. SNPRelate_SideID_unrelated_SAMO_20240822.csv - Text file of sample names and the side of the road the sample came from, input to Step 4
7. SNPRelate_PopID_unrelated_SAMO_v2_20240822.csv - Text file of sample names and the population name where the sample came from, input to Step 4

### <ins>Metadata</ins> 
#### *SNPRelate_PopID_unrelated_20240617.csv*
| Column Name | Data Type | Description | Comments |
| :---: | :---: | :---: | :---: | 
| Column 1 | character | Sample name | NA |
| Column 2 | character | Sampling location name | NA |

#### *SNPRelate_SideID_unrelated_SAMO_20240822.csv*
| Column Name | Data Type | Description | Comments |
| :---: | :---: | :---: | :---: | 
| Column 1 | character | Sample name | NA |
| Column 2 | numeric | Number (1,2) indicating what side (North/South or East/West) the sample came from | North/East = 1, South/West = 2 |

#### *SNPRelate_PopID_unrelated_SAMO_v2_20240822.csv*
| Column Name | Data Type | Description | Comments |
| :---: | :---: | :---: | :---: | 
| Column 1 | character | Sample name | NA |
| Column 2 | character | Sampling location name | NA |
