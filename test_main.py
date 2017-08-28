# reference
# http://pythonhosted.org/gseapy/gseapy_example.html#gsea-example

import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import dataset_ccle
from pdb import set_trace
from numpy import number 
from os.path import exists

def test(): 
    classfile = 'GSEApy/data/P53.cls'
    # 50 2 1
    # #MUT WT
    # MUT MUT MUT MUT MUT MUT MUT MUT MUT MUT MUT MUT MUT MUT MUT MUT MUT MUT MUT MUT \
    # MUT MUT MUT MUT MUT MUT MUT MUT MUT MUT MUT MUT MUT WT WT WT WT WT WT WT WT WT \ 
    # WT WT WT WT WT WT WT WT

    geneexpfile = "GSEApy/data/P53_resampling_data.txt"
    #            NAME   786-0  BT-549  CCRF-CEM  COLO 205    EKVX  HCC-2998  HCT-15  \
    # 0        CTLA2B  111.19   86.22    121.85     75.19  208.62    130.59  124.72
    # 1        SCARA3  460.30  558.34    183.55     37.29  158.00     43.61   80.83
    # 2  LOC100044683   97.25  118.94     81.17    119.51  119.88    107.73  165.57
    # 3          CMBL   33.45   55.10    221.67     50.30   35.12     75.70   84.01
    # 4         CLIC6   35.75   41.26     63.04    219.86   42.53     54.19   86.98

    phenoA,phenoB, class_vector =  gp.parser.gsea_cls_parser(classfile)

    gene_exp = pd.read_table(geneexpfile)

    # gene_exp.head()

    gs_res = gp.gsea(
        data=gene_exp,
        gene_sets='KEGG_2016',
        # cls=class_vector,
        cls=['Control' for i in range(25)] + ['Drug Treatment' for i in range(25)], 
        permutation_type='phenotype',
        outdir='output',
        method='signal_to_noise',
        format='png'
        )

    gsea_results= gs_res.res2d
    # gs_res.res2d.head()

    with plt.style.context('ggplot'):
        gsea_results = gsea_results.reset_index()
        gsea_results.head(5).plot.barh(y='fdr',x='Term',fontsize=10)

        plt.savefig('figure-gsea.pdf')

    # xxx
   

