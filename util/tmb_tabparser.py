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
kaviar_file=sys.argv[3]
out_file=sys.argv[4]
tcga_file=sys.argv[5]

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

for tcgaline in open(tcga_file):
    tcgacol=tcgaline.strip().split()
    chrbp="chr"+str(tcgacol[2])+"_"+str(tcgacol[3])
    chrbp_tcga.append(str(chrbp))
for line in open(tab_file):
    i=i+1
    if i>1:
        col=line.strip().split("\t")
        chrom=col[0]
        bp=col[1]
        chrombp=str(chrom)+"_"+str(bp)
        alt=col[4]
        exac_freq=col[26]
        kg_freq=col[27]
        cosmic=col[28]
        total_callers=col[13]
        snpeff_impact=col[15]
        transcript_id=col[17].split(".")[0]
        abq=col[7]
        iden=str(chrom)+"_"+str(bp)+"_"+str(alt)
        if(int(total_callers)>=2 and abq!="undef"):
            if(snpeff_impact not in exclusions and int(abq)>=30): #and cosmic!="."): #and str(chrombp) in chrbp_tcga):
                if str(transcript_id) in canonical_refseq:
                    out_line=chrom+"\t"+bp+"\t"+str(exac_freq)+"\t"+str(kg_freq)+"\t"+str(total_callers)+"\t"+snpeff_impact+"\t"+transcript_id+"\t"+str(cosmic)+"\t"+str(abq)+"\n"
                    if (exac_freq=="." and kg_freq=="."):
                        id_dict[iden]=str(out_line)
                        #out.write(out_line)
                    if (exac_freq!="." or kg_freq!="."):
                        if (exac_freq=="." and kg_freq!="."):
                            if (float(kg_freq)<=0.01):
                                #print ("kg"+str(kg_freq))
                                #id_dict[iden]=out_line
                                out.write(out_line)
                        if (exac_freq!="." and kg_freq=="."):
                            if (float(exac_freq)<=0.01):
                                #print ("ex"+str(exac_freq))
                                #id_dict[iden]=out_line
                                out.write(out_line)
                    if (exac_freq!="." and kg_freq!="."):
                        if (float(exac_freq)<=0.01 and float(kg_freq)<=0.01):
                            #id_dict[iden]=out_line
                            out.write(out_line)
for kav in open(kaviar_file):
    if "#" in kav:
        continue
    else:
        col1=kav.strip().split("\t")
        chromk=col1[0]
        bpk=col1[1]
        altk=col1[4]
        freqk=col1[7]
        idenk=str(chromk)+"_"+str(bpk)+"_"+str(altk)
        if(idenk in id_dict):
            if (float(freqk)<=0.01):
                out.write(id_dict[idenk])
