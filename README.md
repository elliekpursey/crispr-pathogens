# CRISPR-Cas and MGE detection in bacterial pathogens

[Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline to detect CRISPR-Cas systems using [CRISPRCasTyper](https://github.com/Russel88/CRISPRCasTyper), and search against AMR, plasmid and [ICE](https://db-mml.sjtu.edu.cn/ICEberg/) databases using [ABRicate](https://github.com/tseemann/abricate).    

Dependencies
====== 
`Snakemake` version ≥ 5.1 and [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). This pipeline has been tested in CentOS.

Install
====== 
Navigate to your chosen working directory and clone this github repository.

```shell
git clone https://github.com/elliekpursey/crispr-pathogens.git
```
Input
====== 
Genomic fasta files (.fna), which must be placed in resources > genomes > {species-id}. The *Pseudomonas aeruginosa* PAO1 and PA14 genomes are provided as test files within the folder 287 (the taxonomic id for this species). You must create your own directories within the genomes directory for each {species-id} you want to use.

Usage
======

### Local use
Run the following command from within the crispr-pathogens directory:

```shell
snakemake -s workflow/Snakefile --use-conda --cores 4
```

You may need to adjust the number of cores you use based on your system (--cores argument)

### Cluster use
This repository includes a slurm profile in workflow > profiles > slurm.default, created using [this cookiecutter](https://github.com/Snakemake-Profiles/slurm). You can also make your own profile, and profiles are available for [other job schedulers](https://github.com/Snakemake-Profiles).

The pipeline can then be run with the command:

```shell
snakemake -s workflow/Snakefile --profile profiles/slurm.default -j 4 --use-conda
``` 

where `-j` specifies the number of jobs to submit simultaneously and `--profile` gives the location of your chosen cluster profile.

Output
====== 
The followed compiled output files are produced:

File | Description 
--- | --- 
*amr.csv* | AMR genes for all input files in [ABRicate](https://github.com/tseemann/abricate#output) output format
*crispr_arrays.csv* | CRISPR arrays (crisprs_near_cas.tab from [CRISPRCasTyper](https://github.com/Russel88/CRISPRCasTyper#output-)) for all input files 
*crisprcas.csv* | CRISPR-Cas systems (CRISPR_Cas.tab from [CRISPRCasTyper](https://github.com/Russel88/CRISPRCasTyper#output-)) for all input files 
*genome_lengths.csv* | Lengths of all genomes 
*ICEs.csv* | ICEs for all input files in [abricate](https://github.com/tseemann/abricate#output) output format
*plasmids.csv* | plasmid replicons for all input files in [abricate](https://github.com/tseemann/abricate#output) output format 
*intI1.csv* | intI1 copies for all input files in [abricate](https://github.com/tseemann/abricate#output) output format

Additional analysis
====== 

To detect spacer targets and find which of these targets carry ABR genes, the script *find_spacer_targets.py* was used.

Initial trees were made using mashtree with the script *phylogeny.sh*. The script *root_trees.R* was used to produce the trees in folder *final_trees*.

The following R scripts were used to perform all statistical analyses and generate all figures, supplementary figures and supplementary tables:

Script | Function 
--- | --- 
*descriptive_plots.R* | Descriptive plots, generates figs S1-3
*spacers_descriptive_plots.R* | Descriptive plots, generates fig S5
*CRISPR_binom_models.R* | Binomial GLM of CRISPR-Cas vs. MGE presence, generates figs 1 and S4 and tables S1-4
*bayesian_phylo_models.R* | Bayesian phylogenetically-controlled GLM of ABR counts vs. CRISPR-Cas type and spacer counts, generates fig 2 and tables S5-6
*spacer_target_analysis.R* | Bayesian GLMs of ABR counts vs. counts of spacers targeting MGEs with and without ABR, generates fig 3 and tables S7-8

In addition, the rebuttal analyses in response to reviewer comments are in the folder *rebuttal*. This includes a snakemake pipeline that was used for detecting Type I-C acrs in *Pseudomonas aeruginosa*, as well as R scripts to analyse its output. These produce figures S6 and tables S7 and S8 in the final manuscript.

Citation
======

If you use any of the code in this repository, please cite the paper. 

The work produced by this additional analyis, and using this pipeline, is free to read as a preprint:
https://www.biorxiv.org/content/10.1101/2021.04.12.439454v1

The final version is published in Proceedings of the Royal Society B:
https://royalsocietypublishing.org/doi/10.1098/rstb.2020.0464


