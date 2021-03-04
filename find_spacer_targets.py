#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import subprocess
from Bio.Blast.Applications import NcbiblastnCommandline

# import dataframe of CRISPR arrays near Cas for spacer searching
crisprs = pd.read_csv('results/crispr/crispr_arrays.csv', names=["Contig","CRISPR","Start","End","Consensus_repeat", "N_repeats", "Repeat_len", "Spacer_len_avg", "Repeat_identity",
                 "Spacer_identity", "Spacer_len_sem", "Trusted", "Prediction", "Subtype", "Subtype_probability", "id"]) 

crisprs['path'] = crisprs['id'] + 'spacers/' + crisprs['CRISPR'] + '.fa'

crispr_paths = list(crisprs['path'])

# write all arrays into one file 
spacer_file = 'results/crispr/all_spacers.fna'

with open(spacer_file, "wb") as outfile:
    for path in crispr_paths:
        with open(path, "rb") as infile:
            outfile.write(infile.read())

# run single blast search for all arrays for each db
ICE_db = "ICE_db"
plasmid_db = "plasmid_db"
virus_db = "virus_db"
amr_db = "amr_db"

list_dbs = [ICE_db, plasmid_db, virus_db, amr_db]

# blast against all dbs
for database in list_dbs:
    print('running blast against ' + database)
    blastn_handle = NcbiblastnCommandline(cmd='blastn',
                                          query=spacer_file,
                                          db='resources/blast_dbs/' + database,
                                          evalue=0.001,
                                          out="results/crispr/" + database + "_spacers.csv",
                                          outfmt=6)
    stdout, stderr = blastn_handle()


# read all spacer dfs back in for merging
ICE = pd.read_csv("results/crispr/ICE_db_spacers.csv", sep='\t', names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]) 
plasmid = pd.read_csv("results/crispr/plasmid_db_spacers.csv", sep='\t', names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]) 
virus = pd.read_csv('results/crispr/virus_db_spacers.csv', sep='\t', names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]) 
amr = pd.read_csv('results/crispr/amr_db_spacers.csv', sep='\t', names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]) 

# add columns specifying database type 
ICE['target'] = 'ICE'
plasmid['target'] = 'plasmid'
virus['target'] = 'virus'
amr['target'] = 'amr'

# merge
frames = [ICE, plasmid, virus, amr]

final = pd.concat(frames)

# write to csv
final.to_csv("results/crispr/spacer_targets.csv")

# blast MGEs against AMR database to find MGEs that carry AMR
list_dbs_2 = [ICE_db, plasmid_db, virus_db]

# blast against all dbs
for database in list_dbs_2:
    print('running amr blast against ' + database)
    blastn_handle = NcbiblastnCommandline(cmd='blastn',
                                          query='resources/blast_dbs/amr_db',
                                          db='resources/blast_dbs/' + database,
                                          evalue=0.001,
                                          out="results/crispr/" + database + "_with_AMR.csv",
                                          outfmt=6)
    stdout, stderr = blastn_handle()

# read all dfs back in for merging
ICE_with_AMR = pd.read_csv("results/crispr/ICE_db_with_AMR.csv", sep='\t', names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]) 
plasmid_with_AMR = pd.read_csv("results/crispr/plasmid_db_with_AMR.csv", sep='\t', names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]) 
virus_with_AMR = pd.read_csv('results/crispr/virus_db_with_AMR.csv', sep='\t', names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]) 

# add columns specifying database type 
ICE_with_AMR['target'] = 'ICE'
plasmid_with_AMR['target'] = 'plasmid'
virus_with_AMR['target'] = 'virus'

# merge
frames_AMR = [ICE_with_AMR, plasmid_with_AMR, virus_with_AMR]

final_MGEs_with_AMR = pd.concat(frames_AMR)

# write to csv
final_MGEs_with_AMR.to_csv("results/crispr/MGEs_with_AMR.csv")



