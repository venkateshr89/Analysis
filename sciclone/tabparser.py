#!/usr/bin/python3.5

#usage
### arg1=vcf,
### arg2=target transcript file in gene,transcript format ,
### arg3=tcga file in count,chr,pos

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
                    #inter=chrom+"\t"+bp+"\t"+total_reads+"\t"+var_reads+"\t"+vaf
                    #inter_list.append(inter)
# j=0
# for l in inter_list:
#     j=j+1
#     col=l.strip().split("\t")
#     total=float(col[2])
#     var=float(col[3])
#     ref=col[2]-col[3]
#     print j
#     print(col[0]+"\t"+col[1]+"\t"+ref+"\t"+col[3]+"\t"+col[4])
