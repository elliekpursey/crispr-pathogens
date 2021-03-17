#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 11:17:59 2020

@author: elizabethpursey
"""

import pandas as pd
from Bio import SeqIO
import os
import glob

names = []
lengths = []

for fasta in list(snakemake.input):
    with open(fasta, "r") as handle:
        name = fasta
        length = 0
        for record in SeqIO.parse(handle, "fasta"):
            length = length + len(record)
    lengths.append(length)
    names.append(name)

length_df = pd.DataFrame(data=zip(names, lengths), columns=['id','length'])
length_df.to_csv(snakemake.output[0], index=False)

