import argparse
import sys
import time

import pandas as pd
import numpy as np

from subprocess import Popen, PIPE

parser = argparse.ArgumentParser(description='Create germline BLAST database')
parser.add_argument('input_file', help="tsv containing information on germline calls")
parser.add_argument('-outdir', default='.', help="output directory (default: working directory)")
parser.add_argument('-samplename', default='default', help='name of sample')
args = parser.parse_args()

FILENAME = args.input_file
outdir = args.outdir
samplename = args.samplename


if samplename == 'default':
    samplename=FILENAME.split("/")[-1].split('_v_germlines.tsv')[0]


OUTFILE="{}/{}_v_germlines.fasta".format(outdir, samplename)

df = pd.read_table(FILENAME)
df = df.sort_values(['db_call','mismatch'], ignore_index=True)

novel_allele_counter=0
with open(OUTFILE,'w') as out_writer:
    for it, row in df.iterrows():
        if row.db_call is None:
            novel_allele_counter += 1
            seq_id = "{}*novel{}".format(row.v_family,novel_allele_counter)
        else:
            seq_id = row.db_call
            if row.mismatch == 0:
                pass
            else:
                seq_id = "{}_{}".format(seq_id,row.mutations)
       
        out_writer.write(">{}\n".format(seq_id))
        out_writer.write("{}\n".format(row.sequence))

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

