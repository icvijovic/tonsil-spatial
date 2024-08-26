import sys
import argparse

import numpy as np
from Levenshtein import hamming, distance
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Calculate distance matrix for strings in fasta file')
parser.add_argument('input', help="input (fasta file)")
parser.add_argument('output', help="output (distance matrix)")
parser.add_argument('-hamming', action='store_true')
args = parser.parse_args()

INPUT_FILENAME = args.input
OUTPUT_FILENAME = args.output
FIXED_LENGTH = args.hamming

sys.stderr.write(f"Working on sequences in {INPUT_FILENAME} \n")
seqs = [str(record.seq) for record in SeqIO.parse(INPUT_FILENAME, "fasta")]

n = len(seqs)

sys.stderr.write(f"This dataset contains {n} sequences \n")

mem = round((n/1024)**2*100)/100
sys.stderr.write(f"Estimated memory requirement:\t {mem}MB\n")

D = np.zeros((n,n), np.uint8)

for i in range(n):
    for j in range(i):
        if FIXED_LENGTH:
            d = hamming(seqs[i], seqs[j])
        else:
            d = distance(seqs[i], seqs[j])
        d = min(d, 255)
        d = np.uint8(d)
        D[i,j] = d
        D[j,i] = d

np.save(f'{OUTPUT_FILENAME}', D, allow_pickle=False)
