#0;136;0c!/usr/bin/python3.5

#usage
### arg1=vcf,
### arg2=target transcript file in gene,transcript format ,
### arg3=tcga file in count,chr,pos

import sys
import argparse
import os
import re


###Exclusions from snpeff annotation
exclusions=["3_prime_UTR_variant", "5_prime_UTR_variant", "downstream_gene_variant", \
            "intron_variant", "non_coding_exon_variant", "synonymous_variant", \
            "upstream_gene_variant", "splice_donor_variant&intron_variant", \
            "splice_region_variant&intron_variant", "splice_region_variant&synonymous_variant"]

#Command line Arguments
tab_file=sys.argv[1]
transcript_file=sys.argv[2]
tn_vcf_file=sys.argv[3]
out_file=sys.argv[4]


canonical_refseq=[]
chrbp_tcga=[]
id_dict={}

for line1 in open(transcript_file):
    word1=line1.strip().split()
    transcriptid=word1[0]
    canonical_refseq.append(transcriptid)

out=open(out_file,'w')
inter_list=[]
i=0

for line in open(tab_file):
    i=i+1
    if i>1:
        col=line.strip().split("\t")
        chrom=col[0]
        bp=col[1]
        chrombp=str(chrom)+"_"+str(bp)
        alt=col[4]
        total_callers=col[13]
        snpeff_impact=col[15]
        transcript_id=col[17].split(".")[0]
        abq=col[7]
        var_reads=col[29]
        vaf=col[30]
        total_reads=col[41]
        iden=str(chrom)+"_"+str(bp)+"_"+str(alt)
        if(int(total_callers)>=2):
            if(str(vaf)!="undef" and var_reads!="." and total_reads!="."): #not in exclusions and int(abq)>=30): #and cosmic!="."): #and str(chrombp) in chrbp_tcga):
                if str(transcript_id) in canonical_refseq:
                    vafr=float(vaf)*100
                    out_line=chrom+"\t"+bp+"\t"+"\t"+total_reads+"\t"+var_reads+"\t"+str(vafr)+"\n"
                    id_dict[iden]=str(out_line)
for kav in open(tn_vcf_file):
    if "#" in kav:
        continue
    else:
        col1=kav.strip().split("\t")
        chromk=col1[0]
        bpk=col1[1]
        altk=col1[4]
        idenk=str(chromk)+"_"+str(bpk)+"_"+str(altk)
        if(idenk in id_dict):
            out.write(id_dict[idenk])
