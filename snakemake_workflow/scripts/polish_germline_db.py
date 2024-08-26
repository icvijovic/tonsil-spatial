import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

from copy import deepcopy

from Levenshtein import distance, hamming
from subprocess import Popen, PIPE

from pacbio_vdj_utils.sequence_utils import fasta_to_dict
import argparse

parser = argparse.ArgumentParser(description='Polishes germline database')
parser.add_argument('input_file', help="airr-formatted file containing lineage info")
parser.add_argument('-germline_db', required=True, help="database containing constant region sequences")
parser.add_argument('-imgt_db', required=True, help="path to existing database for species")
parser.add_argument('-imgt_allele_info', required=True, help="path to igblast internal data")
parser.add_argument('-outdir', default='.', help="output directory (default: working directory)")
args = parser.parse_args()

FILENAME = args.input_file
GERMLINE_IGHV_DB_PATH = args.germline_db
IMGT_DB_PATH = args.imgt_db
IMGT_ALLELE_INFO_PATH = args.imgt_allele_info
outdir = args.outdir

samplename=GERMLINE_IGHV_DB_PATH.split("/")[-1].split('.fasta')[0]
OUTFILE = "{}/{}_polished.fasta".format(outdir,samplename)

GERMLINE_IGHV_DB = fasta_to_dict(GERMLINE_IGHV_DB_PATH)
IMGT_DB = fasta_to_dict(IMGT_DB_PATH)

###############################
with open(IMGT_ALLELE_INFO_PATH) as f:
    header = f.readline().rstrip()
allele_columns = ["_".join(x.lower().split(" ")) for x in header[1:-1].split(", ")]

imgt_allele_df = pd.read_table(IMGT_ALLELE_INFO_PATH, comment='#',
                                                      header=None,
                                                      names=allele_columns, 
                                                      delim_whitespace=True)
imgt_allele_df = imgt_allele_df.set_index("gene/allele_name")


df=pd.read_table(FILENAME, usecols=['vdj_sequence','lineage_id','v_call',
                                    'v_sequence_alignment', 
                                    'v_germline_alignment', 
                                    'v_db_call', 'lineage_id', 'cdr3', 'v_mismatch', 'tissue'])
germline_genes = df.v_db_call.unique()

df['cdr3_length'] = df.cdr3.str.len()

df.v_call = df.v_call.str.split(",").map(lambda x: x[0])
df['imgt_db_mismatch'] = df.apply(lambda x: hamming(x.v_sequence_alignment, x.v_germline_alignment), axis=1)

germline_hits = df[df['imgt_db_mismatch'] == 0].groupby(['v_call','v_db_call'])[['lineage_id']].nunique()

germline_hits = germline_hits.reset_index()
sequences_by_imgt_gene = df.groupby('v_call')['vdj_sequence'].nunique()
lineages_by_db_gene = df.groupby('v_db_call')['lineage_id'].nunique()
unmutated_lineages_by_db_gene = df[df.v_mismatch==0].groupby('v_db_call')['lineage_id'].nunique()
otherwise_called_to = df.groupby('v_call')['v_db_call'].agg(set)
germline_hits['nseqs'] = germline_hits.v_call.map(sequences_by_imgt_gene)
germline_hits = germline_hits.rename(columns={'lineage_id':'nlins'})


print('identical hits:', (germline_hits.v_call == germline_hits.v_db_call).sum())
revised = 0
missing = 0

germline_hits=germline_hits[germline_hits.v_call != germline_hits.v_db_call]
hypermutant = germline_hits.apply(lambda x: x.v_call == x.v_db_call.split("_")[0], axis=1)

print("revised:", hypermutant.sum())
germline_hits = germline_hits[~hypermutant]
germline_hits['other_calls'] = germline_hits.v_call.map(otherwise_called_to)

distance_dict = {}
templated_seqs = {}

for it, hit in germline_hits.iterrows():

        imgt_call = IMGT_DB[hit.v_call]
        db_call = GERMLINE_IGHV_DB[hit.v_db_call]

        try:
            cdr3_position = imgt_allele_df.loc[hit.v_call,'fwr3_stop'].astype(int)
        except KeyError:
            sys.stderr.write(f"Gene {hit.v_call} missing from db\n")
            gene=hit.v_call.split('*')[0]
            related_alleles = [x for x in imgt_allele_df.index if x.startswith(gene.rstrip('D'))]
            if len(related_alleles) == 0:
                sys.stderr.write(f"Found no related alleles. Rearrangement likely not functional. Continuing without clipping at CDR3\n")
                cdr3_position = -1
            else:
                sys.stderr.write(f"Found related alleles: str({related_alleles})\n")
                cdr3_positions = [imgt_allele_df.loc[x,'fwr3_stop'].astype(int) for x in related_alleles]
                cdr3_positions = sorted([(cdr3_positions.count(x),x) for x in set(cdr3_positions)])
                cdr3_position = cdr3_positions[0][1]
        
        d = distance(imgt_call[1:cdr3_position], db_call[:cdr3_position-1])
        distance_dict.update({hit.v_call:d})
        templated_seqs.update({hit.v_call:imgt_call[:cdr3_position]})

germline_hits['d'] = germline_hits.v_call.map(distance_dict)
germline_hits['alt_germline_lineages'] = germline_hits.v_db_call.map(unmutated_lineages_by_db_gene)
germline_hits['r'] = germline_hits.nlins / germline_hits.alt_germline_lineages

germline_hits = germline_hits.sort_values(['d','nlins','r'], ascending=False, ignore_index=True)
new_ds = {}
new_rs = {}
for it, row in germline_hits.iterrows():
    if it == 0:
        new_ds.update({it:row['d']})
        new_rs.update({it:row['r']})
    else:
        new_min = sorted([(distance(templated_seqs[row.v_call],
                                           templated_seqs[germline_hits.loc[it2,'v_call']]),
                            it2)
                                    for it2 in range(0, it)])[0]

        new_ds.update({it:min([new_min[0], row['d']])})
        if new_min[0] < row['d']:
            new_rs.update({it:germline_hits.loc[it,'nlins']/
                              germline_hits.loc[new_min[1],'nlins']})
        else:
            new_rs.update({it:germline_hits.loc[it,'r']})

germline_hits['d_new'] = germline_hits.index.map(new_ds)
germline_hits['r_new'] = germline_hits.index.map(new_rs)

germline_filter = (germline_hits[['d','nlins']].max(axis=1) > 4) \
                              & (germline_hits['d'] > 0) \
                              & (germline_hits['nlins'] > 1) \
                              & (germline_hits['d_new'] * germline_hits['r_new'] > 0.3)
print("{} of {} pass filter".format(germline_filter.sum(), germline_filter.shape[0]))

germline_hits = germline_hits[germline_filter]
germline_hits['nalts'] = germline_hits.other_calls.str.len()
print(germline_hits[['v_call', 'v_db_call', 'd', 'nlins', 'r', 'd_new', 'r_new', 'nalts','nseqs']])

polished_database = deepcopy(GERMLINE_IGHV_DB)
for allele in germline_hits.v_call.values:
    polished_database.update({allele:templated_seqs[allele]})

with open(OUTFILE, 'w') as writer:
    for key in sorted(polished_database.keys()):
        writer.write(">{}\n{}\n".format(key,polished_database[key]))

blast_command = "makeblastdb -in {} -dbtype nucl".format(OUTFILE)
p = Popen(blast_command, stdout=PIPE, stderr=PIPE, shell=True, executable='/bin/bash')
blast_out, blast_err = p.communicate()

if len(blast_err) > 0:
    print("Uh oh, encountered error:")
    print(blast_err.decode('utf-8'))
    sys.exit(1)
else:
    print("Success!")
    print(blast_out.decode('utf-8'))
