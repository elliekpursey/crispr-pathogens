#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 11:17:59 2020

@author: elizabethpursey
"""

import pandas as pd

for input_file in list(snakemake.input):
    results['id'] = input_file
    results = pd.read_csv(input_file, sep = '\t')
    results.to_csv(snakemake.output[0], mode='a', header=False)
