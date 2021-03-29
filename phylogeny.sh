#!/bin/bash

#SBATCH --partition=hmq
#SBATCH --job-name="s_mash"
#SBATCH --time=250:00:00
#SBATCH --output standard_sketch_size_mash-%j.out
#SBATCH --mail-user=ep458@exeter.ac.uk   
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=24
##SBATCH --nodes=1

module load Conda/Python3/3.7.2

conda create --name phylo_env python=3.7.0 
source activate phylo_env
conda install -c conda-forge -c bioconda mash=2.0.0
conda install -c bioconda mashtree 
conda install -c bioconda ncbi-genome-download 

cd resources

# SP

# find outliers for removal later
mash dist genomes/1314/*.fna | tee tree/1314_dist.tab > /dev/null
mkdir tree/reference_genome_SP

# download outgroup (Streptococcus agalactiae)
ncbi-genome-download --formats fasta --species-taxids 1311 bacteria --flat-output --output-folder tree/reference_genome_SP --verbose --refseq-categories representative 

gunzip tree/reference_genome_SP/*.gz

mkdir tree/SP_sketches

cd genomes/1314

for genome in *.fna
    do mash sketch $genome -o ../../tree/SP_sketches/${genome%.fna}
    done

cd ../../tree/reference_genome_SP

for genome in *.fna
    do mash sketch $genome -o ../SP_sketches/${genome%.fna}
    done

cd ../

mkdir SP_tmp

mashtree --tempdir SP_tmp --numcpus 12 SP_sketches/*.fna > 1314.dnd

# SE

# find outliers for removal later
mash dist 1282/*.fna | tee 1282_dist.tab > /dev/null
mkdir reference_genome_SE

# download outgroup (SA)
ncbi-genome-download --formats fasta --species-taxids 1280 bacteria --flat-output --output-folder reference_genomes --verbose --refseq-categories reference 

gunzip reference_genome_SE/*.gz

mkdir SE_sketches

cd 1282
for genome in *.fna
    do mash sketch $genome -o ../SE_sketches/${genome%.fna}
    done

cd ../

cd reference_genome_SE

for genome in *.fna
    done

cd ../

mashtree --numcpus 12 SE_sketches/*.msh > 1282.dnd

# PA

# find outliers for removal later
mash dist 287/*.fna | tee 287_dist.tab > /dev/null
mkdir reference_genome_PA

# download outgroup 
ncbi-genome-download --formats fasta --species-taxids 480 bacteria --flat-output --output-folder reference_genome_PA --verbose --refseq-categories representative 

gunzip reference_genome_PA/*.gz

mkdir PA_sketches 

cd 287
for genome in *.fna
    do mash sketch $genome -o ../PA_sketches/${genome%.fna}
    done

cd ../

cd reference_genome_PA

for genome in *.fna
    do mash sketch $genome -o ../PA_sketches/${genome%.fna}
    done

cd ../

mashtree --numcpus 12 PA_sketches/*.msh > 287.dnd

# AB

mkdir reference_genome_AB

# find outliers for removal later
mash dist 470/*.fna | tee 470_dist.tab > /dev/null

mkdir AB_sketches 

# download outgroup 
ncbi-genome-download --formats fasta --species-taxids 480 bacteria --flat-output --output-folder reference_genome_AB --verbose --refseq-categories representative
ncbi-genome-download --formats fasta --species-taxids 287 bacteria --flat-output --output-folder reference_genome_AB --verbose --refseq-categories reference
ncbi-genome-download --formats fasta --species-taxids 352 bacteria --flat-output --output-folder reference_genome_AB --verbose --refseq-categories representative

gunzip reference_genome_AB/*.gz

cd 470
for genome in *.fna
    do mash sketch $genome -o ../AB_sketches/${genome%.fna}
    done

cd ../

cd reference_genome_AB

for genome in *.fna
    do mash sketch $genome -o ../AB_sketches/${genome%.fna}
    done

cd ../

mashtree --numcpus 12 AB_sketches/*.msh > 470.dnd

# EFm

mkdir tree/reference_genome_Efm

# find outliers for removal later
mash dist genomes/1352/*.fna | tee tree/1352_dist.tab > /dev/null

mkdir tree/Efm_sketches 

# download outgroup (E. faecalis)
ncbi-genome-download --formats fasta --species-taxids 1351 bacteria --flat-output --output-folder tree/reference_genome_Efm --verbose --refseq-categories representative

gunzip tree/reference_genome_Efm/*.gz

cd genomes/1352

for genome in *.fna
    do mash sketch $genome -o ../../tree/Efm_sketches/${genome%.fna}
    done

cd ../../tree/reference_genome_Efm

for genome in *.fna
    do mash sketch $genome -o ../Efm_sketches/${genome%.fna}
    done

cd ../

mkdir Efm_tmp

mashtree --tempdir Efm_tmp --numcpus 12 Efm_sketches/*.msh > 1352.dnd

cd ../

# EFs
mkdir tree/reference_genome_Efs

# find outliers for removal later
mash dist genomes/1351/*.fna | tee tree/1351_dist.tab > /dev/null

mkdir tree/Efs_sketches 

# download outgroup (E. faecium)
ncbi-genome-download --formats fasta --species-taxids 1352 bacteria --flat-output --output-folder tree/reference_genome_Efs --verbose --refseq-categories representative

gunzip tree/reference_genome_Efs/*.gz

cd genomes/1351

for genome in *.fna
    do mash sketch $genome -o ../../tree/Efs_sketches/${genome%.fna}
    done

cd ../../tree/reference_genome_Efs

for genome in *.fna
    do mash sketch $genome -o ../Efs_sketches/${genome%.fna}
    done

cd ../

mkdir Efs_tmp

mashtree --tempdir Efs_tmp --numcpus 12 Efs_sketches/*.msh > 1351.dnd
