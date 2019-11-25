# -*- coding: utf-8 -*-
"""
2019-11-15
Pavlos Bousounis
Calculate recovery of ClinVar pathogenic/likely-pathogenic variants in RefSeq exonic regions (GRCh37) by slop length.
"""

from datetime import datetime
import glob
import numpy as np
import os
import pandas as pd
import pathlib
import pybedtools
from pybedtools import BedTool
import re


# today's date
today = datetime.today().strftime('%Y-%m-%d')

# specify and set base directory
basedir = '/Users/pbousounis/Experiments/2019-10-29_hg19mod/2019-11-10_RefSeq-ClinVar_GRCh37_slop_region_recovery'
os.chdir(basedir)

# define function bed_overlap()
def bed_overlap(ref_bed_filepath, test_bed_filepath, dir_out='.'):
    
    """ Given two bed files, ref_bed and test_bed, perform bedtools intersect -c to return the number of 
    test_bed regions that overlap any regions in the ref_bed file. Returns the ref_bed file with a column 
    of overlap counts """
    
    cwd = os.getcwd()
    
    # specify the reference bed file
    ref_bedtool = BedTool(ref_bed_filepath)
    prfx_ref = ref_bed_filepath.split('/')[-1]
    prfx_ref = prfx_ref.split('.')[0]
    
    # specify the new ClinVar bed file
    test_bedtool = BedTool(test_bed_filepath)
    prfx_test = test_bed_filepath.split('/')[-1]
    prfx_test = prfx_test.split('.')[0]

    # specify name/path of output bed file
    bed_out = dir_out + '/{}_IN_{}.bed'.format(prfx_test, prfx_ref)
    
    # run bedtools intersect to get all test_bed regions NOT found in ref_bed (-v option)
    ref_in_test = test_bedtool.intersect(b=ref_bedtool, c=True)
    
    # save the bed overlap file
    ref_in_test.saveas(bed_out)
    
    # confirm file saved    
    if os.path.isfile(bed_out):
        print('Success!\nFile saved to: \n{}.\n'.format(os.path.join(cwd, bed_out)))

    return(ref_in_test)


# the clinvar bed file to intersect
bed_a_file = 'data/2019-11-07-GRCh37_path-likely_path_headless.bed'

# list of slopped refseq bed files to intersect
slop_files = glob.glob('data/2019-11-13_slop/2019-11-13_RefSeqGRCh37_exon_slop_bed/*RefSeqGRCh37*slop*bed')[6:]


bedcols = ['chr', 'start', 'end', 'name']
df = pd.DataFrame(columns=bedcols)
for file in slop_files:

    # import file as dataframe
    tmp = pd.read_csv(file, sep='\t', header=None, names=bedcols)
        
    # extract transcript IDs and gene names to their own columns
    tmp['tx_id'] = tmp['name'].str.extract(r'(\w+-)(N(M|R)_\d+)')[1]
    #tmp['gene'] = tmp['name'].str.split('-').str[0]
    
    # add column with slop file ID
    slop_length = re.findall(r'slop\d+', file)[0]
    tmp['slop_len'] = slop_length
    
    # split by transcript accession (RefSeq)
    tx_table = pd.DataFrame(columns=bedcols) 
    for tx in tmp.tx_id.unique():
        
        # subset by tx_id
        tmp_tx = tmp[tmp.tx_id == tx].copy()
        
        # extract exon numbers (arbitrary; relative to transcript)
        tmp_tx.loc[:, 'exon_num'] = tmp_tx.name.str.split('-', n=1).str[-1]
        
        # convert to bedtool 
        tx_bedtool = BedTool.from_dataframe(tmp_tx).sort()
        
        # bedtools merge -c 7 -o 'collapse'
        tx_merged = tx_bedtool.merge(c=7, o='collapse')
        print('\n{} regions merged.\n'.format(tx))

        
        # convert to dataframe
        tx_df = pd.read_table(tx_merged.fn, names=bedcols)
        
        # add tx_id column
        tx_df.loc[:, 'tx_id'] = tx_df['name'].str.split(':').str[0]

        # append to master table
        tx_table = pd.concat([tx_df, tx_table], sort=False)
        print('Success! Merged {0} regions added to {1} master table.\n'.format(tx, slop_length))
    
    file_out = os.path.join('output', (today + '_RefSeqGRCh37_' + slop_length + '_merged.bed'))
    tx_table.to_csv(file_out, sep='\t', index=None)



## test overlap of slop50 merged (above)
        
tmp_overlap = bed_overlap(bed_a_file, '2019-11-17_RefSeqGRCh37_slop50_transcripts.bed').to_dataframe()
tmp_overlap.columns = ['chrom', 'start', 'end', 'exon_name', 'tx_acc', 'exon_recovery']
tx_ovl = tmp_overlap[tmp_overlap['exon_recovery'] != 0 ]
tx_ovl.loc['tx_recovery'] = tx_ovl.groupby("tx_acc")["exon_recovery"].transform('sum')


# import and process ClinVar bed file
cv_path = '/Users/pbousounis/Experiments/2019-10-29_hg19mod/2019-11-10_RefSeq-ClinVar_GRCh37_slop_region_recovery/data/2019-11-07_ClinVar-GRCh37_path-likely_path.bed'
cv_bed = pd.read_csv(cv_path, sep='\t')

# remove chr MT entries and extract gene names to new column
cv_bed = cv_bed[cv_bed.chr != 'MT']
cv_bed['gene'] = cv_bed['name'].str.split('_').str[1]

# get ClinVar variants by gene
cv_bed['gene_recovery'] = cv_bed.groupby('gene')['gene'].transform('count')





# create/specify output directory
slop_int_out = os.path.join(('../output/' + today + '_RefSeqGRCh37_ClinVar_slop-intersect-merge'))
pathlib.Path(slop_int_out).mkdir(exist_ok=True)

# create empty dict 
for i in slop_files:
    slop_num = re.findall(r'\d+', i)[-1]
    fname = '{}_RefSeqGRCh37_ClinVar_recovery_{}.bed'.format(today, slop_num)
    bed_overlap(bed_a_file, i, slop_int_out)
    
    
# list overlapped slop files
ovl_slop_files = glob.glob(slop_int_out + '/*RefSeqGRCh37_slop*IN*GRCh37_path-likely_path_headless.bed')



