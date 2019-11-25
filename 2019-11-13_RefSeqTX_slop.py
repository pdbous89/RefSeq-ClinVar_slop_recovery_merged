#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 22:10:15 2019

@author: pbousounis
"""

os.chdir('/Users/pbousounis/Experiments/2019-10-29_hg19mod/RefSeq-ClinVar_GRCh37_slop_region_recovery/data')

bed_a_filepath = 

rsbed_file = '2019-11-08_RefSeq-GRCh37_latest_genomicGFF3.bed'
rsbed = pd.read_csv(rsbed_file, sep='\t', low_memory=False, names=['chr', 'start', 'end', 'name'])

# extract genes
rsbed['gene'] = rsbed['name'].str.split('-').str[0]

# extract transcripts 
rsbed['rs_tx'] = rsbed['name'].str.extract(r'(\w+-)(N(M|R)_\d+)')[1]

# create tx_slop file output directory
slop_tx_dir_out = os.path.join(basedir, (today + "_RefSeqGRCh37_exon_TX_slop_bed"))
pathlib.Path(slop_tx_dir_out).mkdir(exist_ok=True)


for tx in rsbed.rs_tx.unique():
    
    tmp_df = rsbed[rsbed.rs_tx == tx]
    tmp_df['chr'] = 'chr' + tmp_df['chr']
    
    
    bedtool_tx = BedTool.from_dataframe(tmp_df)
    
    # create tx output directory
    tx_dir_out = os.path.join(slop_tx_dir_out, tx)
    pathlib.Path(tx_dir_out, ).mkdir(exist_ok=True)
    
    # for each N in the interval 0, 5, 10, .., 155 (non-inclusive) run slop adding N bases to each end of each region
    for i in list(range(5, 155, 5)):
        
        # create the output filename
        file_out = os.path.join(tx_dir_out, '{}_RefSeqGRCh37_slop{}_{}.bed'.format(today, i, tx))
            
        # perform slop: add 'i' bases to each end of each region
        bedtool_tx.slop(g=hg19_genome, b=i).saveas(file_out)
            
        # check file save
        if os.path.isfile(os.path.join(os.getcwd(), file_out)):
            print('\nSuccess! {} saved.\n'.format(file_out))

            
tmp = clinvar_bedtool.intersect(b=bedtool_tx, v=True)
tmp.head()