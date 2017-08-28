import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import dataset_ccle
from pdb import set_trace
from numpy import number 
from os.path import exists


df1 = dataset_ccle.load_gexp()
df2 = df1[df1.columns[1:4]]
df3 = df2 

df3 = df3.dropna()

df3 = df3.loc[df2.index[0:500]]

cls_vector = ['MUT' for x in range(1)]
cls_vector += ['WT' for x in range(1)]

# xxx

# set_trace()

df = df3 

# df.drop_duplicates(subset=df.columns[0], inplace=True) #drop duplicate gene_names.
# df.set_index(keys=df.columns[0], inplace=True)
# df.dropna(how='all', inplace=True)                     #drop rows with all NAs
# df2 = df.select_dtypes(include=[number])
# #drop any genes which std ==0

# # set_trace() 
# df_std =  df2.groupby(by=cls_vector, axis=1).std()
# df2 =  df2[~df_std.isin([0]).any(axis=1)]
# df2 = df2 + 0.00001 # we don't like zeros!!!

# run call
# enrichr library are supported by call module. Just provide the name
# you may also provide a gene_sets file in gmt format, just like GSEA do.
gs_res = gp.gsea(data=df3, gene_sets='KEGG_2016', cls=cls_vector, 
                 permutation_type='phenotype',
                 outdir='gsea_reprot', method='signal_to_noise', format='png')


xxx