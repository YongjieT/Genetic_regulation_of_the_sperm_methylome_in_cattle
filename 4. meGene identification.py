import sys
import pandas as pd
import numpy as np
import tensorqtl
from tensorqtl import genotypeio, cis, trans


input_phenotype_file=sys.argv[1]
input_genotype_file=sys.argv[2]
input_chr=sys.argv[3]

#input phe
phenotype_file_df = pd.read_table(input_phenotype_file, sep="\t",index_col = 3)
#phe matrix
phenotype_pos_df = phenotype_file_df[["#Chr", "start"]]
phenotype_df = phenotype_file_df.drop(columns=['#Chr','start','end','gid','strand'])
phenotype_pos_df.columns = ['chr', 'pos']

#input GT
genotype_file_df = pd.read_table(input_genotype_file, sep="\t",index_col = 2)
#GT matrix
variant_df = genotype_file_df[["#CHROM", "POS"]]
genotype_df = genotype_file_df.drop(columns=['#CHROM','POS','REF','ALT','QUAL','FILTER','INFO','FORMAT'])
variant_df.columns = ['chrom', 'pos']


# permutations
cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, window=500000)

filename = f"tensorQTL_permut_input_chr_{input_chr}.txt"
cis_df.to_csv(filename, index=True, sep='\t')
