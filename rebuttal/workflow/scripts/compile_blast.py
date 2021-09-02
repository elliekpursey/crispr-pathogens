import pandas as pd
import numpy

df = pd.DataFrame(columns=['qseqid','seqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qseq'])
df.to_csv(snakemake.output[0], mode='w')

for input_file in list(snakemake.input):
    results = pd.read_csv(input_file, sep = '\t', names=['qseqid','seqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qseq'])
    results['genome'] = input_file
    filtered = results[results['evalue'] < 0.001]
    filtered.to_csv(snakemake.output[0], header=False, mode='a')

