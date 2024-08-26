import sys
from functools import partial

from Bio import SeqIO, bgzf
import regex

import gzip
import argparse

import time

import pandas as pd
import numpy as np

from pacbio_vdj_utils.blast_utils import *


parser = argparse.ArgumentParser(description='Wrapper for v sequence mutation calling')
parser.add_argument('input_file', help="airr-formatted file")
parser.add_argument('-germline_db', required=True, help="database containing constant region sequences")
parser.add_argument('-outdir', default='.', help="output directory (default: working directory)")
parser.add_argument('-samplename', help="name of new file (default: same as old)")
args = parser.parse_args()

FILENAME = args.input_file
BLAST_DB = args.germline_db
outdir = args.outdir
samplename = args.samplename

if samplename is None:
    samplename=FILENAME.split("/")[-1].split('.tsv.gz')[0]

#########################################################################################################

df = pd.read_table(FILENAME, low_memory=False)

df['sequence_id'] = df['sequence_id'].astype(str)

#use igblast call to determine the constant region start position

df['v_sequence'] = df.apply(lambda x: x.sequence[int(x.v_sequence_start):int(x.cdr3_start)], axis = 1)
v_sequences = df['v_sequence'].unique()
v_sequence_ids = {x:it for it, x in enumerate(v_sequences)}

df['v_seq_id'] = df.v_sequence.map(lambda x: v_sequence_ids[x]).astype(str)

sys.stderr.write("Starting with %d unique sequences...\n" % len(v_sequences))

additional_blastn_options = "-word_size 9 -dust no -penalty -1 -gapopen 3 -gapextend 2 "
#now pipe to blastn
chunk_size = 300

fasta_records = [">{}\n{}".format(seq_id, seq) for seq, seq_id in v_sequence_ids.items()]
blastn_results = []

start = time.time()
sys.stderr.write("Piping sequences to blastn...\n")
for i in range(0, len(fasta_records), chunk_size):
	chunk_begin = i
	chunk_end = chunk_begin + chunk_size

	if chunk_end < len(fasta_records):
		query_string = "\n".join(fasta_records[chunk_begin:chunk_end])
	else:
		query_string = "\n".join(fasta_records[chunk_begin:])

	blastn_out, blastn_err = pipe_to_blastn(query_string, BLAST_DB, 
											evalue="20",
											additional_blastn_options=additional_blastn_options)
	blastn_results.append(return_best_match(blastn_out))

	if i % 5000 == 0:
		sys.stderr.write("[%ds] %d records processed...\n" % (time.time() - start, i))

blastn_results = pd.concat(blastn_results)
blastn_results.sequence_id = blastn_results.sequence_id.astype(str)
sys.stderr.write("Aligned {n} sequences successfully.\n".format(n = blastn_results.shape[0]))

sys.stderr.write("Parsing BTOP strings...\n")
blastn_results['mutations'] = blastn_results.apply(lambda x: parse_btop(x.btop, x.sstart), axis=1)
blastn_results = blastn_results.rename(columns = {"sequence_id" : "v_seq_id",
						  "match"   : "v_db_call",
						  "pident"  : "v_pident",
						  "length"  : "v_match_length",
						  "mismatch": "v_mismatch",
						  "gapopen" : "v_gapopen",
						  "qstart"  : "v_qstart",
						  "qend"    : "v_qend",
						  "sstart"  : "v_sstart",
						  "send"   : "v_send",
						  "evalue"    : "v_evalue",
						  "btop"    : "v_btop",
						  "mutations": "v_mutations"
						  })


df = df.merge(blastn_results, on = "v_seq_id", how = "outer")
sys.stderr.write("[%ds] Done! Saving...\n" % (time.time() - start))
df.to_csv('{}/{}.tsv.gz'.format(outdir,samplename), sep = '\t', index = False)
sys.stderr.write("Output written to: {}/{}.tsv.gz \n".format(outdir, samplename))
