#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 11:24:29 2019

@author: pbousounis
"""


from datetime import datetime
import os
import pandas as pd
import pybedtools
from pybedtools import BedTool
import re


# today's date
today = datetime.today().strftime('%Y-%m-%d')

# set base directories
basedir = '/Users/pbousounis/Experiments/2019-10-29_hg19mod/2019-11-10_RefSeq-ClinVar_GRCh37_slop_region_recovery'
slomp_dir = os.path.join(basedir, 'output/2019-11-18_slops_merged-by-tx')
data_dir = os.path.join(basedir, 'data')

# get list of merged slop tx files
slomp_files = os.listdir(slomp_dir)

#file = os.path.join(slomp_dir, slomp_files[1])

gene2refseq_id = os.path.join(basedir, 'output/2019-11-18_gene2refseq_tx.tsv')
clinvar_bed_file = os.path.join(data_dir, '2019-11-07_ClinVar-GRCh37_path-likely_path.bed')


def slomp_proc(slompfile, clinvar_file, refseq_slop_file):
    
    # import merged slop bed file
    slomp_df = pd.read_csv(os.path.join(slomp_dir, slompfile, sep='\t'))
    
    
    # import gene symbol - RefSeq tx acc conversion table
    gene2refseq = pd.read_csv(gene2refseq_id, sep='\t')
    
    
    # import clinvar p-lp bed file
    clinvar_bed = pd.read_csv(clinvar_bed_file, sep='\t')
    
    # extract gene names; calculate gene recovery
    clinvar_bed['gene'] = clinvar_bed.name.str.split('_').str[1]
    clinvar_bed['gene_recovery'] = clinvar_bed.groupby('gene')['gene'].transform('count')
    gene_recovery = clinvar_bed[['gene', 'gene_recovery']].drop_duplicates()
    
    
    
    # run bedtools intersect -c 
    slop_bedtool = BedTool(slomp_df)
    cv_bedtool = BedTool(clinvar_bed_file)
    ovl = slop_bedtool.intersect(b=cv_bedtool, c=True).to_dataframe()
    ovl.columns = ['chr', 'start', 'end', 'exon_name', 'tx_id', 'slop_recovery_ex']
    ovl = ovl[ovl.slop_recovery_ex != 0].sort_values(['chr', 'start'])
    
    # add gene names and gene recovery values to intersected table
    ovl = ovl.merge(gene2refseq, how='left', on='tx_id')
    ovl = ovl.merge(gene_recovery, how='left', on='gene')
    
    # build transcript slop recovery table
    tx_slop_recovery = ovl[['gene', 'tx_id', 'gene_recovery', 'slop_recovery_ex']]
    tx_slop_recovery['slop_recovery_tx'] = tx_slop_recovery.groupby('tx_id')['slop_recovery_ex'].transform('sum')
    
    # calculate slopped transcript percent recovery by gene
    tx_slop_recovery['slop_recovery%'] = tx_slop_recovery.loc[:, 'slop_recovery_tx'] / tx_slop_recovery.loc[:, 'gene_recovery'] * 100
    
    # add column with slop file ID
    slop_length = re.findall(r'slop\d+', slomp_df)[0]
    tx_slop_recovery['slop_len'] = slop_length