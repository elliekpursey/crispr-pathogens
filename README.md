# CRISPR-Cas and MGE detection in bacterial pathogens

[Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline to detect CRISPR-Cas systems using [CRISPRCasTyper](https://github.com/Russel88/CRISPRCasTyper), and search against AMR, plasmid and [ICE](https://db-mml.sjtu.edu.cn/ICEberg/) databases using [abricate](https://github.com/tseemann/abricate).    

Dependencies
====== 
`Snakemake` version ≥ 5.1 and [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). 

Install
====== 
Clone this github repository.

```shell
git clone https://github.com/elliekpursey/crispr-pathogens.git
```
Input
====== 
Genomic fasta files (.fna), which must be placed in resources > genomes. The *Pseudomonas aeruginosa* PAO1 genome is provided as a test file.

Output
====== 
The followed compiled output files are produced:

File | Description 
--- | --- 
*amr.csv* | AMR genes for all input files in [abricate](https://github.com/tseemann/abricate) output format
*crispr_arrays.csv* | CRISPR arrays (crisprs_near_cas.tab from [CRISPRCasTyper](https://github.com/Russel88/CRISPRCasTyper)) for all input files 
*crisprcas.csv* | CRISPR-Cas systems (CRISPR_Cas.tab from [CRISPRCasTyper](https://github.com/Russel88/CRISPRCasTyper)) for all input files 
*genome_lengths.csv* | Lengths of all genomes 
*ICEs.csv* | ICEs for all input files in [abricate](https://github.com/tseemann/abricate) output format
*plasmids.csv* | plasmid replicons for all input files in [abricate](https://github.com/tseemann/abricate) output format 
*intI1.csv* | intI1 copies for all input files in [abricate](https://github.com/tseemann/abricate) output format

Further analysis
====== 

Citation
======
