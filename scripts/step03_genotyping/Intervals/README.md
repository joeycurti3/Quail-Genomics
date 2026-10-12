These BED files define scaffold-specific genomic intervals used for variant calling. Intervals were generated from the California Quail (CAQU) reference genome assembly (GCA_023055505.1_bCalCai1.0.p) using the script `step03_a_CAQU_GetIntervals_20240918.R`. The script uses the NCBI assembly report (`GCA_023055505.1_bCalCai1.0.p_assembly_report.txt`) as input, available [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/055/505/GCA_023055505.1_bCalCai1.0.p/GCA_023055505.1_bCalCai1.0.p_assembly_report.txt).

Each BED file contains a single genomic interval spanning the full length of an individual scaffold. Intervals follow the standard BED format (`name start end`), using zero-based, half-open coordinates. Consequently, start coordinates are always 0, and end coordinates correspond to the total scaffold length in base pairs.

Scaffolds assigned to CAQU sex chromosomes were excluded from these intervals, specifically:

- JALIRH010000009.1
- JALIRH010000023.1
- JALIRH010000025.1
- JALIRH010000031.1
- JALIRH010000041.1
- JALIRH010000049.1
- JALIRH010000052.1
- JALIRH010000066.1
- JALIRH010000062.1
- JALIRH010000085.1
- JALIRH010000094.1
- JALIRH010000103.1

### <ins>Metadata</ins>
#### *scaffold_[insert_number].bed*
| Column Name | Data Type | Description | Comments |
| :---: | :---: | :---: | :---: |
| name | character | GenBank accession number identifying the scaffold in the CAQU reference genome | Corresponds to `GenBank-Accn` in the NCBI assembly report |
| start | integer | Starting genomic coordinate of the interval | Always 0; BED uses zero-based coordinates |
| end | integer | Ending genomic coordinate of the interval | Corresponds to the scaffold length in base pairs; BED end coordinates are exclusive |

**Example:** `scaffold_68.bed`

```text
JALIRH010000068.1 0 2524086
```

#### *GCA_023055505.1_bCalCai1.0.p_assembly_report.txt*

| Column Name | Data Type | Description | Comments |
| :---: | :---: | :---: | :---: |
| Sequence-Name | character | Name of the chromosome, scaffold, or sequence in the genome assembly | Assigned by the assembly submitter |
| Sequence-Role | character | Role of the sequence within the genome assembly | e.g., assembled-molecule, unlocalized-scaffold, unplaced-scaffold |
| Assigned-Molecule | character | Chromosome or molecule to which the sequence is assigned | NA if the sequence is not assigned to a chromosome or molecule |
| Assigned-Molecule-Location/Type | character | Location or type of the assigned molecule | e.g., Chromosome, Mitochondrion, or na |
| GenBank-Accn | character | GenBank accession number identifying the sequence | Includes accession version |
| Relationship | character | Relationship between the GenBank and RefSeq sequences | e.g., = indicates identical sequences; NA if no relationship is reported |
| RefSeq-Accn | character | RefSeq accession number corresponding to the GenBank sequence | NA if no RefSeq accession is available |
| Assembly-Unit | character | Assembly unit to which the sequence belongs | e.g., Primary Assembly, non-nuclear |
| Sequence-Length | integer | Length of the sequence in base pairs | Units: bp |
| UCSC-style-name | character | UCSC Genome Browser-style sequence identifier | e.g., chr1, chr2, chrM; NA if unavailable |
