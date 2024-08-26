import sys, os

import time
import argparse
import numpy as np
import pandas as pd

import scipy.spatial.distance as ssd

from cython_packages.MST_SINGLE_UINT8 import *
from pacbio_vdj_utils.cluster_vdj import *

parser = argparse.ArgumentParser(description='Cluster sequences for large datasets')
parser.add_argument('input_file', help="airr-formatted file")
#parser.add_argument('-matrixfiles', required=True, nargs='+', 
#                    help='.npy files containing pairwise hamming distance matrices for all cdr3 groups')
parser.add_argument('-outdir', default='.', help="output directory (default: working directory)")
parser.add_argument('-matrixdir', default='.', help="input distance matrix directory "
													"(default: same as output directory)")
parser.add_argument('-samplename', default='.', help='samplename (default: filename retained)')
parser.add_argument('-max_distance_from_vstart', type=int, required=False, default=0)
parser.add_argument('-max_distance_from_cstart', type=int, required=False, default=0)
parser.add_argument('-fractional_cutoff', type=float, required=False, default=0.15)

args = parser.parse_args()

FILENAME = args.input_file
OUTDIR = args.outdir
#MATRIXFILES = args.matrixfiles
MATRIXDIR = args.matrixdir
V_SSTART_MAX = args.max_distance_from_vstart
C_SSTART_MAX = args.max_distance_from_cstart
FRACTIONAL_CUTOFF = args.fractional_cutoff
SAMPLENAME = args.samplename

if SAMPLENAME==".":
    SAMPLENAME=FILENAME.split("/")[-1].split('.')[0]

if MATRIXDIR == ".":
    MATRIXDIR = OUTDIR


#########################################################################################################
def get_current_memory_usage():
    ''' Memory usage in GB '''
    with open('/proc/self/status') as memusage_file:
        memusage = memusage_file.read().split('VmRSS:')[1].split('\n')[0][:-3]
    return int(memusage.strip())/1024/1024

#########################################################################################################

NROWS = None
# read only relevant columns

df = pd.read_table(FILENAME, usecols=['sequence',
                              'v_family',
                              'v_sequence_alignment',
                              'j_sequence_alignment',
                              'v_sequence_start',
                              'v_germline_start',
                              'v_germline_alignment',
                              'j_germline_alignment',
                              'cdr3_start',
                              'cdr3_end',
                              'cdr3',
                              'j_sequence_end'], nrows=NROWS)

df['cdr3_length'] = df['cdr3'].str.len()

sys.stderr.write(f"Starting with {df.shape[0]} reads...\n")

df = df[df.v_germline_start <= V_SSTART_MAX + 1]
try:
    df = df[df['c_sstart'] <= C_SSTART_MAX + 1]
except KeyError:
    pass
sys.stderr.write(f"Keeping {df.shape[0]} reads that contain entire VDJ region...\n")

#df = df[df['v_family'].str.startswith("IGH")]
#sys.stderr.write("Subsetting to {} reads that map to heavy chain...\n".format(df.shape[0]))

df['vdj_sequence'] = df.apply(lambda x: x.sequence[int(x.v_sequence_start)-1:
                                       int(x.j_sequence_end)], axis = 1)

df['v_templated_len'] = df['cdr3_start'] - df['v_sequence_start']
df['j_templated_len'] = df['j_sequence_end'] - df['cdr3_end']

# drop reads with gap in templated alignment
df['sum_gap_len'] = df.v_sequence_alignment.map(lambda x: x.count("-"))
df['sum_gap_len'] = df['sum_gap_len'] + df.j_sequence_alignment.map(lambda x: x.count("-"))
df = df[df['sum_gap_len'] < 1]

sys.stderr.write(f"Keeping {df.shape[0]} reads with no deletions in templated sequences...\n")

df['sum_insertion_len'] = df.v_germline_alignment.map(lambda x: x.count("-"))
df['sum_insertion_len'] = df['sum_insertion_len'] + df.j_germline_alignment.map(lambda x: x.count("-"))
df = df[df['sum_insertion_len'] < 1]

sys.stderr.write(f"Keeping {df.shape[0]} reads with no insertions in templated sequences...\n")

unique_vdjs = df[['v_family',
                  'vdj_sequence',
                  'v_templated_len',
                  'j_templated_len',
                  'cdr3',
                  'cdr3_length']].drop_duplicates(ignore_index=True)

unique_vdjs['cdr3_group'] = unique_vdjs['v_family'] + "_" + unique_vdjs.cdr3_length.astype(str)

cdr3_group_sizes = unique_vdjs.cdr3_group.value_counts()
unique_vdjs['cdr3_group_size'] = unique_vdjs.cdr3_group.map(cdr3_group_sizes)

sys.stderr.write("Verifying that all distance matrices are available...\n")
for cdr3_group in unique_vdjs.cdr3_group.unique():
    vfam, cdr3_len = cdr3_group.split("_")[0], cdr3_group.split("_")[1]
    binary_matrix_filename = f'{MATRIXDIR}/{SAMPLENAME}_{vfam}_{cdr3_len}_cdr3.npy'
    if os.path.exists(binary_matrix_filename):
        pass
    else:
        sys.stderr.write(f"Cannot find the following file {binary_matrix_filename}. Aborting...\n")
        sys.exit(1)

sys.stderr.write(f"Clustering {unique_vdjs.shape[0]} unique variable sequences...\n")

TOTAL_CLUSTERS = 0
unique_vdjs['cluster_id'] = -1

start = time.time()

for cdr3_group in unique_vdjs.cdr3_group.unique():
    vfam, cdr3_len = cdr3_group.split("_")[0], int(cdr3_group.split("_")[1])

    subset = unique_vdjs[unique_vdjs['cdr3_group'] == cdr3_group]

    cdr3_seqs = subset.cdr3.values
    subset_idx = subset.index

    sys.stderr.write(f"Processing VDJ sequence subset: v_family={vfam}, cdr3_length={cdr3_len}, n={subset.shape[0]}\n")

    binary_matrix_filename = f'{MATRIXDIR}/{SAMPLENAME}_{vfam}_{cdr3_len}_cdr3.npy'

    Ds = np.load(binary_matrix_filename,allow_pickle = False)
    n=Ds.shape[0]
    diag_indices = np.diag_indices(n)
    Ds[diag_indices] = np.uint8(255)
    closest_d = np.min(Ds, axis=1)

    Ds[diag_indices] = np.uint8(0)

    if n > 1:
        # report on memory usage when dealing with very large lineages
        if n > 10000:
            matrix_size =  sys.getsizeof(Ds)/(1024)**3
            matrix_size = round(matrix_size*100)/100
            total_memory_usage = round(get_current_memory_usage()*100)/100
            
            print(f"Clustering large subset:\n n = {n}, matrix_size = {matrix_size} GB", file=sys.stderr)
            print(f"    Current memory usage: {total_memory_usage} GB", file=sys.stderr)
            print(f"    Reading {binary_matrix_filename}", file=sys.stderr)

        #cluster CDR3s
        cutoff = np.uint8(FRACTIONAL_CUTOFF*cdr3_len)
        #Ds = ssd.squareform(Ds)
        upper_triangular_indices = np.triu_indices(n, 1)
        condensed_distance_matrix = Ds[upper_triangular_indices]

        new_cluster_ids = get_cluster_ids(condensed_distance_matrix,
                                          cutoff=cutoff,
		                                  method='single')
    else:
        new_cluster_ids = np.zeros(1)

    unique_vdjs.loc[subset_idx, 'cluster_id'] = TOTAL_CLUSTERS + new_cluster_ids
    TOTAL_CLUSTERS += max(new_cluster_ids) + 1

    sys.stderr.write(f"[{time.time() - start}s]  found {new_cluster_ids.max() + 1} new cdr3 clusters, continuing...\n")

sys.stderr.write(f"Clustering CDR3s took {round(time.time() - start)} seconds\n")

unique_vdjs['templated_v'] = unique_vdjs.apply(lambda x: x.vdj_sequence[0:int(x.v_templated_len)], axis=1)
unique_vdjs['templated_j'] = unique_vdjs.apply(lambda x: x.vdj_sequence[-int(x.j_templated_len):], axis=1)

unique_vdjs.to_csv(f"{OUTDIR}/{SAMPLENAME}_unique_vdjs_cdr3_clusters.tsv.gz", sep= '\t')

LARGE_GROUP_CUTOFF=8*10**3

sys.stderr.write("Preparing to cluster templated sequences within cdr3 clusters...\n")

cluster_sizes = unique_vdjs.cluster_id.value_counts()
unique_vdjs['cluster_size'] = unique_vdjs.cluster_id.map(cluster_sizes)

large_groups = unique_vdjs['cluster_id'][
                    unique_vdjs['cluster_size'] > LARGE_GROUP_CUTOFF].unique()
small_groups = unique_vdjs['cluster_id'][
                    unique_vdjs['cluster_size'] <= LARGE_GROUP_CUTOFF].unique()

start = time.time()

sys.stderr.write(
    f"Calculating distance matrices for {len(unique_vdjs.cluster_id.unique())} unique clusters...\n")

for cluster_id in large_groups:
    subset = unique_vdjs[unique_vdjs['cluster_id'] == cluster_id]
    
    for seq_element in ['templated_v', 'templated_j']:
        BINARY_MATRIX_FILENAME = f'{OUTDIR}/{SAMPLENAME}_{cluster_id}_{seq_element}.npy'

        templated_seqs = subset[seq_element].values
        n = len(templated_seqs)
        subset_idx = subset.index

        sys.stderr.write(
            f"Processing VDJ sequence subset: cluster_id={cluster_id}, n={n}\n")


        sys.stderr.write("\t writing sequences to disk for cluster execution...\n" )
        FASTA_FILENAME = f'{OUTDIR}/{SAMPLENAME}_{cluster_id}_{seq_element}.fasta'

        with open(FASTA_FILENAME, 'w') as writer:
            items = [f">{it}\n{templated}" for it, templated in enumerate(templated_seqs)]
            writer.write("\n".join(items))

        sys.stderr.write(f"\t sequences written to {FASTA_FILENAME}!\n" )

for cluster_id in small_groups:
    subset = unique_vdjs[unique_vdjs['cluster_id'] == cluster_id]
    
    for seq_element in ['templated_v', 'templated_j']:
        BINARY_MATRIX_FILENAME = f'{OUTDIR}/{SAMPLENAME}_{cluster_id}_{seq_element}.npy'

        templated_seqs = subset[seq_element].values
        n = len(templated_seqs)
        subset_idx = subset.index

        sys.stderr.write(
            f"Processing VDJ sequence subset: cluster_id={cluster_id}, n={n}\n")

        sys.stderr.write("\t computing distance matrix locally...\n")

        D = np.zeros((n,n),np.uint8)

        for i in range(n):
            for j in range(i):
                d = distance(templated_seqs[i], templated_seqs[j])
                d = np.uint8(min(d, 255))
                D[i,j] = d
                D[j,i] = d

        np.save(BINARY_MATRIX_FILENAME, D, allow_pickle=False)

    sys.stderr.write(f"\t\t this took {time.time() - start} seconds\n")


sys.stderr.write(f'Done! Execution took {(time.time() - start)} seconds\n')
