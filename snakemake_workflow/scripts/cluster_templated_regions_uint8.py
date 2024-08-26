import sys
import os

import time
import argparse
import numpy as np
import pandas as pd

import scipy.spatial.distance as ssd

from cython_packages.MST_SINGLE_UINT8 import *
from pacbio_vdj_utils.cluster_vdj import *

parser = argparse.ArgumentParser(description='Cluster sequences for large datasets')
parser.add_argument('-airr', help="airr-formatted file containing contig information")
parser.add_argument('-cdr3clusters', help="table with unique vdjs and corresponding cdr3 clusters")
#parser.add_argument('-matrixfiles', required=True, nargs='+',
#                    help='.npy files containing pairwise levenshtein distance matrices')
parser.add_argument('-outdir', default='.', help="output directory (default: working directory)")
parser.add_argument('-matrixdir', default='.', help="input distance matrix directory "
                                                    "(default: same as output directory)")
parser.add_argument('-samplename', default='.', help='samplename (default: filename retained)')
parser.add_argument('-fractional_cutoff', type=float, required=False, default=0.15)

args = parser.parse_args()

CONTIGFILE = args.airr
UNIQUE_VDJ_FILE=args.cdr3clusters
OUTDIR = args.outdir
#MATRIXFILES = args.matrixfiles
MATRIXDIR = args.matrixdir
FRACTIONAL_CUTOFF = args.fractional_cutoff
SAMPLENAME = args.samplename

if SAMPLENAME==".":
    SAMPLENAME=CONTIGFILE.split("/")[-1].split('.')[0]

if MATRIXDIR == ".":
    MATRIXDIR = OUTDIR


####################################################################################################
def get_current_memory_usage():
    ''' Memory usage in GB '''
    with open('/proc/self/status') as memusage_file:
        memusage = memusage_file.read().split('VmRSS:')[1].split('\n')[0][:-3]
    return int(memusage.strip())/1024/1024

####################################################################################################

# read only unique vdjs first

start = time.time()
unique_vdjs = pd.read_table(UNIQUE_VDJ_FILE)

sys.stderr.write("Verifying that all distance matrices are available...\n")
for cluster_id in unique_vdjs.cluster_id.unique():
    for s in ['v', 'j']:
        BINARY_MATRIX_FILENAME = f'{MATRIXDIR}/{SAMPLENAME}_{cluster_id}_templated_{s}.npy'
        if os.path.exists(BINARY_MATRIX_FILENAME):
            pass
        else:
            sys.stderr.write(f"Cannot find the following file {BINARY_MATRIX_FILENAME}. Aborting...\n")
            sys.exit(1)

sys.stderr.write(f"Clustering {unique_vdjs.shape[0]} unique variable sequences...\n")

unique_vdjs['templated_v'] = unique_vdjs.apply(lambda x: x.vdj_sequence[0:int(x.v_templated_len)], axis=1)
unique_vdjs['templated_j'] = unique_vdjs.apply(lambda x: x.vdj_sequence[-int(x.j_templated_len):], axis=1)
unique_vdjs['templated_vj'] = unique_vdjs['templated_v'] + "+" + unique_vdjs['templated_j']

TOTAL_LINEAGES = 0
unique_vdjs['lineage_id'] = -1

for cluster_id in unique_vdjs.cluster_id.value_counts().index:
    subset = unique_vdjs[unique_vdjs['cluster_id'] == cluster_id]

    V_BINARY_MATRIX_FILENAME = f'{MATRIXDIR}/{SAMPLENAME}_{cluster_id}_templated_v.npy'
    J_BINARY_MATRIX_FILENAME = f'{MATRIXDIR}/{SAMPLENAME}_{cluster_id}_templated_j.npy'
   
    templated_seqs = subset.templated_vj.values
    n = len(templated_seqs)
    subset_idx = subset.index

    longest_templated_sequence = subset.templated_vj.str.len().max()

    sys.stderr.write(
        f"Processing VDJ sequence subset: cluster_id={cluster_id}, n={n}\n")

    Ds = np.load(V_BINARY_MATRIX_FILENAME, allow_pickle = False)
    DJs = np.load(J_BINARY_MATRIX_FILENAME, allow_pickle = False)
    
    R = np.uint8(255) * np.ones((n,n), np.uint8) - Ds
    Ds = Ds + DJs
    Ds[R<DJs] = np.uint8(255)

    n=Ds.shape[0]
    diag_indices = np.diag_indices(n)
    Ds[diag_indices] = np.uint8(255)
    closest_d = np.min(Ds, axis=1)

    Ds[diag_indices] = np.uint8(0)

    if n > 1:
        # report on memory usage when dealing with very large lineages
        if n > 8000:
            matrix_size =  sys.getsizeof(Ds)/(1024)**3
            matrix_size = round(matrix_size*100)/100
            total_memory_usage = round(get_current_memory_usage()*100)/100

            print(f"Clustering large subset:\n n = {n}, matrix_size = {matrix_size} GB", file=sys.stderr)
            print(f"    Current memory usage: {total_memory_usage} GB", file=sys.stderr)
            print(f"    Reading {BINARY_MATRIX_FILENAME}", file=sys.stderr)

        #cluster CDR3s
        cutoff = np.uint8(FRACTIONAL_CUTOFF*longest_templated_sequence)
        #Ds = ssd.squareform(Ds)
        upper_triangular_indices = np.triu_indices(n, 1)
        condensed_distance_matrix = Ds[upper_triangular_indices]

        new_lineage_ids = get_cluster_ids(condensed_distance_matrix,
                                          cutoff=cutoff,
                                          method='single')
    else:
        new_lineage_ids = np.zeros(1)
    unique_vdjs.loc[subset_idx, 'lineage_id'] = TOTAL_LINEAGES + new_lineage_ids
    TOTAL_LINEAGES += max(new_lineage_ids) + 1

    sys.stderr.write(f"[{time.time() - start}s]  found {new_lineage_ids.max() + 1} new lineages...\n")

sys.stderr.write(f"Clustering templated sequences took {round(time.time() - start)} seconds\n")


vdj_families = unique_vdjs[['vdj_sequence','lineage_id']].set_index('vdj_sequence').to_dict()
vdj_families = vdj_families['lineage_id']

# read whole dataframe and append lineage ids
df = pd.read_table(CONTIGFILE)
df['vdj_sequence'] = df.apply(lambda x: x.sequence[int(x.v_sequence_start)-1:
                                       int(x.j_sequence_end)], axis = 1)

df['lineage_id'] = df.vdj_sequence.map(vdj_families)
df = df[df.lineage_id.notna()]

df.to_csv(f'{OUTDIR}/{SAMPLENAME}.tsv.gz', sep = '\t', index = False)
