import sys, glob
import argparse
import numpy as np
import pandas as pd



######################## PATH CONFIG ################################


parser = argparse.ArgumentParser()
parser.add_argument('-vdj', required=True,  help="vdj data frame")
parser.add_argument('-gex', required=True,  help="gex data frame")
parser.add_argument('-output',required=True, help="output filepath")

args = parser.parse_args()

vdj_loc = parser.vdj
gex_loc = parser.gex
out_loc = parser.output 

######################## ############ ################################

sample_relationships = pd.read_table(sample_relationships_loc)
sample_relationships = sample_relationships.set_index('sample_uid_vdj')

# read and prepare gex label df
gex_label_df = pd.read_table(gex_loc)
gex_label_df = gex_label_df.rename(columns={'Unnamed: 0':'barcode', 
                                            'sample_uid':'sample_uid_gex'})

print('dropping non-bursa cells...', file=sys.stderr)
print(f'original dataframe contained {gex_label_df.shape[0]} cells', file=sys.stderr)
gex_label_df = gex_label_df[gex_label_df.sample_uid.str.startswith("TBd")]
print(f'after dropping, dataframe contains {gex_label_df.shape[0]} cells', file=sys.stderr)

gex_label_df['cb'] = gex_label_df.barcode.str.split("-").map(lambda x: x[0])


# read and prepare vdj df
igh_df = pd.read_table(vdj_loc)
igh_df['sample_uid_vdj'] = igh_df['sample_uid'].copy()
igh_df['sample_uid_gex'] = igh_df['sample_uid_vdj'].map(lambda x: sample_relationships.loc[x,'sample_uid_gex'])

########### MERGE ################
igh_df = igh_df.merge(label_df, on=['cb', 
                                    'sample_uid_gex',
                                    'donor',
                                    'tissue'], 
                                how='outer',
                                suffixes=('_vdj','_gex'))

print('Merged dataframe has ', igh_df.shape[0], 'total cells', file=sys.stderr)
# # REMOVE CELLS THAT ARE NEITHER B CELLS BY GEX NOR VDJ
# igh_df = igh_df[igh_df.vdj_sequence.notna() | igh_df.bcell_subtype.notna()]
# print(igh_df.shape[0], 'after removing non-B cells')

igh_df.to_csv(out_loc, sep='\t', index=False)