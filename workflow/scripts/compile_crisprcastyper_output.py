#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 11:17:59 2020

@author: elizabethpursey
"""

import pandas as pd

# create output files to ensure pipeline runs in case no CRISPR-Cas systems detected
open(snakemake.output[0], 'a').close()
open(snakemake.output[1], 'a').close()

# compile CRISPR-Cas systems
for input_file in list(snakemake.input):
    try:
        results = pd.read_csv(input_file + "/CRISPR_Cas.tab", sep = '\t')
        print(results)
        results['id'] = input_file
        results.to_csv(snakemake.output[0], mode='a', header=False)
    except:
        print("no result for " + input_file)

# compile CRISPR-Cas arrays
for input_file in list(snakemake.input):
    try:
        results = pd.read_csv(input_file + "/crisprs_near_cas.tab", sep = '\t')
        results['id'] = input_file
        results.to_csv(snakemake.output[1], mode='a', header=False)
    except:
        print("no result for " + input_file)
