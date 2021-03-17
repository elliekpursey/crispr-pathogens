# CRISPR-Cas and MGE detection in bacterial pathogens

[Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline to detect CRISPR-Cas systems using [CRISPRCasTyper](https://github.com/Russel88/CRISPRCasTyper), and search against AMR, plasmid and [ICE](https://db-mml.sjtu.edu.cn/ICEberg/) databases using [abricate](https://github.com/tseemann/abricate).    

Dependencies
====== 
`Snakemake` version ≥ 5.1 and [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). This pipeline has been tested in Ubuntu.

Install
====== 
Navigate to your chosen working directory and clone this github repository.

```shell
git clone https://github.com/elliekpursey/crispr-pathogens.git
```
Input
====== 
Genomic fasta files (.fna), which must be placed in resources > genomes > {species-id}. The *Pseudomonas aeruginosa* PAO1 genome is provided as a test file within the folder 287 (its taxonomy id). You must create and name directories within the genomes directory for each species/group of genomes you want to use.

Usage
======

### Local use
Run the following command from within the crispr-pathogens directory:

```shell
snakemake -s workflow/Snakefile --use-conda
```

### Cluster use
This repository includes a slurm profile in workflow > profiles > slurm.default, created using [this cookiecutter](https://github.com/Snakemake-Profiles/slurm). You can also make your own profile, and profiles are available for [other job schedulers](https://github.com/Snakemake-Profiles).

The pipeline can then be run with the command:

```shell
snakemake -s workflow/Snakefile --profile profiles/slurm.default -j 3 --use-conda
``` 

where `-j` specifies the number of jobs to submit simultaneously and `--profile` gives the location of your chosen cluster profile.

Output
====== 
The followed compiled output files are produced:

File | Description 
--- | --- 
*amr.csv* | AMR genes for all input files in [abricate](https://github.com/tseemann/abricate#output) output format
*crispr_arrays.csv* | CRISPR arrays (crisprs_near_cas.tab from [CRISPRCasTyper](https://github.com/Russel88/CRISPRCasTyper#output-)) for all input files 
*crisprcas.csv* | CRISPR-Cas systems (CRISPR_Cas.tab from [CRISPRCasTyper](https://github.com/Russel88/CRISPRCasTyper#output-)) for all input files 
*genome_lengths.csv* | Lengths of all genomes 
*ICEs.csv* | ICEs for all input files in [abricate](https://github.com/tseemann/abricate#output) output format
*plasmids.csv* | plasmid replicons for all input files in [abricate](https://github.com/tseemann/abricate#output) output format 
*intI1.csv* | intI1 copies for all input files in [abricate](https://github.com/tseemann/abricate#output) output format

Further analysis
====== 

Citation
======
