#!/usr/bin/env python
# coding: utf-8

import os
from Bio import SeqIO
from Bio import Entrez


MEP = {"MEP-1-DXPS": "2.2.1.7",
       "MEP-2-DXPR": "1.1.1.267",
       "MEP-3-IspD": "2.7.7.60",
       "MEP-4-IspE": "2.7.1.148",
       "MEP-5-IspF": "4.6.1.12",
       "MEP-6-IspG1": "1.17.7.1",
       "MEP-6-IspG2": "1.17.7.3",
       "MEP-7-IspH": "1.17.7.4"
      }


MEV = {"MEV-1-ACAT": "2.3.1.9",
       "MEV-2-HMGS": "2.3.3.10",
       "MEV-3-HMGR": "1.1.1.34",
       "MEV-4-MK": "2.7.1.36",
       "MEV-5-PMK": "2.7.4.2",
       "MEV-6-MDC": "4.1.1.33"
      }


C1020 = {"C1020-1-IPPI": "5.3.3.2",
         "C1020-2-GPPS": "2.5.1.1",
         "C1020-3-FPPS": "2.5.1.10",
         "C1020-4-GGPPS": "2.5.1.29"
        }


TBB_Biosyn = [MEP, MEV, C1020]


Entrez.email = "myuen@mail.ubc.ca"


for p in TBB_Biosyn:
    for k,v in p.items():
        filename = k + ".fasta"

        if not os.path.isfile(filename):
            print(filename)

            out_handle = open(filename, "w")

            # Downloading...
            query = '({0}[EC/RN Number]) AND Spermatophyta[Organism]'.format(v)

            entrez_handle = Entrez.esearch(db="protein", idtype="acc", term=query, retmax=150)

            record = Entrez.read(entrez_handle)
            
            print(k, "\t", v, "\t", record["Count"])

            accs = record['IdList']

            for acc in accs:
                fasta = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
                
                out_handle.write(fasta.read())

            out_handle.close()
            print("Saved")

            entrez_handle.close()
