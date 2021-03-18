#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 11:17:59 2020

@author: elizabethpursey
"""

import pandas as pd

f1 = open(snakemake.output[0], "w")
f1.close()        
f2 = open(snakemake.output[1], "w")
f2.close()

for input_file in list(snakemake.input):
    print(input_file + "CRISPR_Cas.tab")
    try:
        results = pd.read_csv(input_file + "CRISPR_Cas.tab", sep = '\t')
        results['id'] = input_file
        results.to_csv(snakemake.output[0], mode='a', header=False)
    except:
        print("no result for " + input_file)

for input_file in list(snakemake.input):
    try:
        results = pd.read_csv(input_file + "crisprs_near_cas.tab", sep = '\t')
        results['id'] = input_file
        results.to_csv(snakemake.output[1], mode='a', header=False)
    except:
        print("no result for " + input_file)
