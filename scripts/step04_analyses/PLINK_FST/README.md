### <ins>The following files/scripts are in this directory and should be executed in this order:</ins> 

1. step04_a_CAQU_Plink_Fst_20240820.sh - Script to run Plink on all SAMO samples
2. step04_a_CAQU_Plink_Fst_unrelated_20250901.sh - script to run Plink on all unrelated samples
3. step04_b_CAQU_Plink_Fst_Plot_20250901.R - plotting script for both all and unrelated SAMO samples
4. 20240820_CAQU_Locations.csv - Sampling location spatial information
5. popid.txt - input for Plink based on all SAMO samples
6. popid_unrelated.txt - input for Plink based on unrelated SAMO samples

### <ins>Metadata</ins> 
#### *20240820_CAQU_Locations.csv*
| Column Name | Data Type | Description | Comments |
| :---: | :---: | :---: | :---: | 
| sample_name | character | name of the sample used in analysis | NA |
| longitude | numeric | Longitude of the sampling location | EPSG:4269: NAD83 |
| latitude | numeric | Latitude of the sampling location | EPSG:4269: NAD83 |

#### *popid.txt*
| Column Name | Data Type | Description | Comments |
| :---: | :---: | :---: | :---: | 
| Column 1 | numeric | Family ID (FID) field required by PLINK, set to zero because no information about family structure | NA |
| Column 2 | character | Sample name | NA |
| Column 3 | character | Sampling location name | NA |

#### *popid_unrelated.txt*
| Column Name | Data Type | Description | Comments |
| :---: | :---: | :---: | :---: | 
| Column 1 | numeric | Family ID (FID) field required by PLINK, set to zero because no information about family structure | NA |
| Column 2 | character | Sample name | NA |
| Column 3 | character | Sampling location name | NA |
