#!/bin/bash

#SBATCH --partition=hmq  
#SBATCH --time=100:00:00   
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --job-name="ncbi"
#SBATCH --mail-user=ep458@exeter.ac.uk   
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load Conda/Python3/3.7.2
module load BioPerl/1.6.924

# create conda environment, activate it and install ncbi-genome-download and dependencies
conda create --prefix /nobackup/beegfs/home/ISAD/ep458/.conda/envs/ncbi_env python=3.7.0 
source activate ncbi_env
conda install -c bioconda ncbi-genome-download=0.3.0
conda install pandas
conda install -c bioconda blast 
conda install -c conda-forge biopython 

# make output directories
cd proc_b/resources

# mkdir blast_dbs
cd blast_dbs
mkdir virus 

# download all ncbi viral genomes
ncbi-genome-download --formats fasta viral --parallel 24 --flat-output --output-folder virus --verbose

# download plasmid database
wget -O All_plasmids.fna "https://dryad-assetstore-merritt-west.s3.us-west-2.amazonaws.com/ark%3A/13030/m5741nwt%7C1%7Cproducer/All_plasmids.fna?response-content-type=text%2Fplain&X-Amz-Security-Token=IQoJb3JpZ2luX2VjEN%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJGMEQCIHa5d1n7TPp%2FWK5Hwj8XQuFBkQ%2FEh0IDmvokyf3hnpuRAiA%2F368TNVF4LuWeQsZfwhrNfs9cs0mNk9RhpdxYeoAoiCq9AwjY%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F8BEAAaDDQ1MTgyNjkxNDE1NyIMIbnVf%2B8IabJIP%2BW3KpEDMmeDQfI359lbw8NhaT8l8982Qzb60S%2Br37LSg1aWyDlgp6dLuqiwSRea05iQBCTPrEs54pH19MSVgE1LxR%2F22xJdLMu1DjvVtz4ezi4r9l%2Fr8bSYjsGXjFOjn1Ukfx39H3%2BeSCWiYhXja45JomK7vudtVv9LqU67R7YFUH5h%2FSW2WNNm2TbbnvCynHNQKAb9oE1bNrwnE7oouI3vVn4QadDlwAHMHxNCESushfL0O7QYBwxj1dvP9ahcvtHFPNYAnLkcVz4xw2yPzmFBXHMUTFwIm%2F8jcL5XS6XHe9jNYfBgAGq%2FbP9A8dc74FbVHQmRaGbZyU3qOH2y%2BGH4lb8%2BzdaBC6EGMYgFJ2yrW9y9tXowAM2cpNf9T9WpXaZYSuaYNQ%2BVg7oh2W%2B0EgxHgFLXnWcpKLigeSNBVofmgM0EpmT9P6wLub5Cd%2F%2Byqy4%2BLIMF8nQ30T0Payj9u3QV8T1ZR837OGrZUNvTLoPXGwmsD5dBdMKGxjPw8GZ5CrPv5SgEk7Dg5H%2Fyi4jdDo31EU5pzTEwgLaagQY67AGJDfD%2FbgcExH4znpHKwNtFs%2FdZywIQUcYXdGoi6t0%2BpfiPPpDWg8Gg6OMZTu9zO4C1LvGBGxUGpzVgCA6ZTjAgUDMl9jAZ7Es8zj4aj6Pil%2FtISfp1NVFgVDwKURid3PMJ9hpa4uomNFLzSnfsaRjcgTWq0BpzQZZsJr4%2FGNHqMMsMmQq%2BGWQUWqBH6KxMdD9tJTaZZf1nkUuMMN0BJZJRR3Yh4WjySBPPj2kxXLw3H4aBOvQKhAkc4iLR37lNcZmXF1kn6Z3Wz0urNCeAdkrCDQ4oLk2jJ0yi%2FGBLBeWnho%2B15Yc87Rl7f9fuUw%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20210212T163306Z&X-Amz-SignedHeaders=host&X-Amz-Expires=14400&X-Amz-Credential=ASIAWSMX3SNW6G2SS4UI%2F20210212%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Signature=302645dd134462a0cc9567fbdb613635a937797e24a82e273d3eba4413474b43" 

# make viral blast database
cd virus

gunzip *.gz

cat *.fna | makeblastdb -dbtype nucl -title virus_db -out ../virus_db -parse_seqids

cd ../

# make plasmid blast database
makeblastdb -in All_plasmids.fna -dbtype nucl -title plasmid_db -out plasmid_db -parse_seqids

# make ICE blast database
makeblastdb -in ../ICE_seq_all.fa -dbtype nucl -title ICE_db -out ICE_db -parse_seqids

cd ../../

# make outfile for python
touch results/crispr/all_spacers.fna

cd ../

# run python script to blast selected genomes
python find_spacer_targets.py >> spacers.out

