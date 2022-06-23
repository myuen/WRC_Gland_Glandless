#!/usr/bin/env python

import argparse
import re

from Bio import SeqIO
from Bio import SeqUtils

# Search for SDR base on motif desc in Lange & Srividya (2019)
# https://doi.org/10.1093/jxb/ery436

parser = argparse.ArgumentParser(description='Simple python script to search for short chain dehydrogenase (SDR) '
                                             'based on 3 defined motif.')

parser.add_argument("-i", "--input", required=True, help="Peptide sequences in fasta format.")
parser.add_argument("-o", "--output", required=True, help="Output file with putative SDR ID")

args = parser.parse_args()
args = vars(args)

out = open(args["output"], "w")


with open(args["input"], "r") as f:
    for s in SeqIO.parse(f, "fasta"):

        seq_str = str(s.seq)

        # Function cannot interpret '*' as stop codon.  Remove '*' from string
        if seq_str.endswith('*'):
            seq_str = seq_str[:-1:]

        # SDR weight between 27-36 kDa.
        mol_wt = SeqUtils.molecular_weight(seq_str, seq_type="protein")

        # Only enforce upper limit for molecular weight in case the sequence is not full-length

        # SDR with NAD+ cofactor
        # Motif 1: GXXXGXG
        # Motif 2: NAG
        # Motif 3: YXXXK
        if re.search("G[A-Z]{3}G[A-Z]G[A-Z]+D[A-Z]+NAG[A-Z]+N[A-Z]+S[A-Z]+Y[A-Z]{3}K", seq_str) and mol_wt < 40000:
            # SeqIO.write(s, out, "fasta-2line")
            out.write(s.id + "\n")

        # SDR with NADPH cofactor
        # Motif 1: GXXXGXG
        # Motif 2: RXXXR
        # Motif 3: NAG
        # Motif 4: [EY]XXXK
        if re.search("G[A-Z]{3}G[A-Z]G[A-Z]+R[A-Z]{3}R[A-Z]+NAG[A-Z]+[EY][A-Z]{3}K", seq_str) and mol_wt < 40000:
            # print(SeqUtils.molecular_weight(s.seq, seq_type="protein"))
            out.write(s.id + "\n")

out.close()
