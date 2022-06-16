#!/usr/bin/env python

import argparse
import re

from Bio import SeqIO

# Search for SDR base on motif desc in Lange & Srividya (2019)
# https://doi.org/10.1093/jxb/ery436

parser = argparse.ArgumentParser(description="Simple python script to search for medium chain dehydrogenase (MDR) "
                                             "based on 3 defined motif.")

parser.add_argument("-i", "--input", required=True, help="Peptide sequences in fasta format.")
parser.add_argument("-o", "--output", required=True, help="Output file with putative MDR ID")

args = parser.parse_args()
args = vars(args)

out = open(args["output"], "w")


with open(args["input"], "r") as f:
    for s in SeqIO.parse(f, "fasta"):

        seq_str = str(s.seq)

        # MDR with NADP as cofactor
        # Motif 1: GHEXXGXXXXXGXXV
        # Motif 2: GXGXXG
        # Motif 3: STSXXK
        if re.search("GHE[A-Z]{2}G[A-Z]{5}G[A-Z]{2}V[A-Z]+G[A-Z]G[A-Z]{2}G[A-Z]+STS[A-Z]{2}K", seq_str):
            # SeqIO.write(s, out, "fasta-2line")
            out.write(s.id + "\n")

        # MDR with NADPH as cofactor
        # Motif 1: GMPGMTA
        # Motif 2: AXXGXXG
        if re.search("GMPGMTA[A-Z]+A[A-Z]{2}G[A-Z]{2}G", seq_str):
            # SeqIO.write(s, out, "fasta-2line")
            out.write(s.id + "\n")

out.close()
