#!/usr/bin/python3.5

# This script splits multi alt vcf to single alt per line ## test run on kaviar o/p #
#usage
### arg1=vcf which has multiple alternate alleles,
### arg2=o/p file with one alt per record

import sys
import argparse
import os
import re


#Command line Arguments
tab_file=sys.argv[1]
out_file=sys.argv[2]

out=open(out_file,'w')
#inter_list=[]
i=0
for line in open(tab_file):
    if "#" in line:
        out.write(line)
    if "#" not in line:
        col=line.strip().split("\t")
        chrom=col[0]
        bp=col[1]
        rsid=col[2]
        ref=col[3]
        alt=col[4]
        freq=col[7]
        qual=col[5]
        filters=col[6]
        if "," in alt:
            alt1=alt.strip().split(",")
            freqAF=freq.strip().split(";")[0]
            freqall=freqAF.strip().split("=")[1]
            freq1=freqall.strip().split(",")
            print(alt1)
            print(len(alt1))
            for i in range(len(alt1)):
                out.write("chr"+str(chrom)+"\t"+str(bp)+"\t"+str(rsid)+"\t"+str(ref)+"\t"+alt1[i]+"\t"+qual+"\t"+filters+"\t"+freq1[i]+"\n")
        if "," not in alt:
            freqAF=freq.strip().split(";")[0]
            freqall=freqAF.strip().split("=")[1]
            out.write("chr"+str(chrom)+"\t"+str(bp)+"\t"+str(rsid)+"\t"+str(ref)+"\t"+str(alt)+"\t"+str(qual)+"\t"+str(filters)+"\t"+str(freqall)+"\n")
