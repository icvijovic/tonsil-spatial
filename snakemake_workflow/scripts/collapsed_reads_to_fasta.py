import argparse
import pandas as pd
import sys
############################   PARSER   ###############################

parser = argparse.ArgumentParser()
parser.add_argument('-i')
parser.add_argument('-o', default='reads.fasta')

args = parser.parse_args()

infile = args.i
outfile = args.o

###########################################################################

df = pd.read_csv(infile, usecols=['readId', 'cseq', 'collapsed_read'])

#drop sequences with ambiguous bases
df = df[df['cseq'] == df['collapsed_read']]
print(f"Removed reads with ambiguous bases", file=sys.stderr)

assert df['readId'].nunique() == df.shape[0]
assert df['cseq'].nunique() == df.shape[0]

print(f"All readsId's are unique", file=sys.stderr)

with open(outfile, 'w') as f:
    for i, row in df.iterrows():
        f.write(f">{row['readId']}\n{row['cseq']}\n")

print(f"Done. Wrote {df.shape[0]} reads to {outfile}", file=sys.stderr)
    