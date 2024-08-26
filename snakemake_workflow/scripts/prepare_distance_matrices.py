import sys

import time
import argparse
import numpy as np
import pandas as pd

from Levenshtein import hamming

parser = argparse.ArgumentParser(description='Cluster sequences for large datasets')
parser.add_argument('input_file', help="airr-formatted file")
parser.add_argument('-outdir', default='.', help="output fasta and matrix directory "
                                     		"(default: current directory)")
parser.add_argument('-samplename', default='.', help='samplename (default: filename retained)')
parser.add_argument('-max_distance_from_vstart', type=int, required=False, default=0)
parser.add_argument('-max_distance_from_cstart', type=int, required=False, default=0)

args = parser.parse_args()

FILENAME = args.input_file
OUTPUTDIR = args.outdir
V_SSTART_MAX = args.max_distance_from_vstart
C_SSTART_MAX = args.max_distance_from_cstart
SAMPLENAME = args.samplename

if SAMPLENAME==".":
    SAMPLENAME=FILENAME.split("/")[-1].split('.')[0]


###################################################################################

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


LARGE_GROUP_CUTOFF=10**4
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
large_groups = unique_vdjs['cdr3_group'][
					unique_vdjs['cdr3_group_size'] > LARGE_GROUP_CUTOFF].unique()
small_groups = unique_vdjs['cdr3_group'][
					unique_vdjs['cdr3_group_size'] <= LARGE_GROUP_CUTOFF].unique()


start = time.time()

sys.stderr.write(
	f"Calculating distance matrices for {unique_vdjs.shape[0]} unique variable sequences...\n")

for cdr3_group in large_groups:
    vfam, cdr3_len = cdr3_group.split("_")[0], cdr3_group.split("_")[1]

    BINARY_MATRIX_FILENAME = f'{OUTPUTDIR}/{SAMPLENAME}_{vfam}_{cdr3_len}_cdr3.npy'

    subset = unique_vdjs[unique_vdjs['cdr3_group'] == cdr3_group]
    sys.stderr.write(
    	f"Processing VDJ sequence subset: v_family={vfam}, cdr3_length={cdr3_len}, n={len(subset)}\n")
    subset_lengths = subset.cdr3_length.unique()
    print(subset_lengths)
    cdr3_seqs = subset.cdr3.values
    n = len(cdr3_seqs)
    subset_idx = subset.index

    sys.stderr.write("\t writing sequences to disk for cluster execution...\n" )
    FASTA_FILENAME = f'{OUTPUTDIR}/{SAMPLENAME}_{vfam}_{cdr3_len}_cdr3.fasta'

    with open(FASTA_FILENAME, 'w') as writer:
        items = [f">{int(it)}\n{cdr3}" for it, cdr3 in enumerate(cdr3_seqs)]
        writer.write("\n".join(items))

    sys.stderr.write(f"\t sequences written to {FASTA_FILENAME}\n" )

for cdr3_group in small_groups:
    new_start = time.time()
    vfam, cdr3_len = cdr3_group.split("_")[0], cdr3_group.split("_")[1]

    BINARY_MATRIX_FILENAME = f'{OUTPUTDIR}/{SAMPLENAME}_{vfam}_{cdr3_len}_cdr3.npy'

    subset = unique_vdjs[unique_vdjs['cdr3_group'] == cdr3_group]
    sys.stderr.write(
    	f"Processing VDJ sequence subset: v_family={vfam}, cdr3_length={cdr3_len}, n={len(subset)}\n")

    cdr3_seqs = subset.cdr3.values
    n = len(cdr3_seqs)
    subset_idx = subset.index

    sys.stderr.write("\t computing distance matrix locally...\n")

    D = np.zeros((n,n),np.uint8)

    for i in range(n):
        for j in range(i):
            d = hamming(cdr3_seqs[i], cdr3_seqs[j])
            d = np.uint8(min(d, 255))
            D[i,j] = d
            D[j,i] = d

    np.save(BINARY_MATRIX_FILENAME, D, allow_pickle=False)

    sys.stderr.write(f"\t\t this took {time.time() - new_start} seconds\n")


sys.stderr.write(f'Done! Execution took {(time.time() - start)} seconds\n')
