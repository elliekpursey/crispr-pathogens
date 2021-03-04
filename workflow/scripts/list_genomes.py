#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 11:17:59 2020

@author: elizabethpursey
"""

import pandas as pd
import os

genomes = list(snakemake.input)

genome_df = pd.DataFrame(genomes)

genome_df.to_csv(snakemake.output[0], index=False)