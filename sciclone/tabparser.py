#!/usr/bin/python3.5

# this script transforms " ann.all_transcripts.tab" from single sample tumor mutation calling to sciclone accepted VAF.dat format

#usage
### arg1=ann.transcript.tab from Single sample pipeline or T/N pipeline
### arg2= output file

### output file should be in the format with header
#CHROM   POS     TUM1_REF_READS  TUM1_VAR_READS  TUM1_VAF
#chr1    10146   23      3       38.0
#chr1    10146   23      3       38.0
#chr1    10492   2       2       100.0


import sys
import argparse
import os
import re

#Command line Arguments
tab_file=sys.argv[1]
out_file=sys.argv[2]

out=open(out_file,'w')
inter_list=[]
i=0
for line in open(tab_file):
    i=i+1
    if i>1:
        col=line.strip().split("\t")
        chrom=col[0]
        bp=col[1]
        var_reads=(col[29])
        vaf=col[30]
        total_reads=col[41]
        if (var_reads!="."):
            if(total_reads!="."):
                if (str(vaf)!="undef" and int(col[13])>=2):
                    var_reads=(col[29])
                    total_reads=(col[41])
                    vaf1=float(vaf)*100
                    out_line=chrom+"\t"+bp+"\t"+total_reads+"\t"+var_reads+"\t"+str(vaf1)+"\n"
                    out.write(out_line)
                   
