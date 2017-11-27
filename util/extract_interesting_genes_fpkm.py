#!/usr/bin/python3.5


import sys
import argparse
import os
import re
import gzip
import scipy.stats as ss

#Command line Arguments
rnaseq_file=sys.argv[1] ##SQL Epitope file
gene_file=sys.argv[2] ## count file with snv and indel caller count information
out_file=sys.argv[3] ## output with Epitope and count information merged

refseq_dict={}
dict_count={}
dict_rna={}
dict_fpkm={}
na=[]

out=open(out_file,'w')
out.write("Gene_Name\tProtein_Name\tRefseq_ID\tsample_genlab_id\tTranscript_ID\tchr:bp\tFPKM_avg\tFPKM_low\tFPKM_high\tStatus\n")

### Mapping Count file to WES epitope file ###

cf=open(gene_file,'r')
for line in cf:
    word=line.strip().split(",")
    refseq=word[2]
    gene_info="\t".join(word[0:3])
    refseq_dict[refseq]=gene_info
for line_ep in open(rnaseq_file,'r'):
    word_ep=line_ep.strip().split("\t")
    if (len(word_ep)!=1):
        refseq_id=word_ep[1]
        data=word_ep[0]+"\t"+word_ep[1]+"\t"+word_ep[7]+"\t"+word_ep[10]+"\t"+word_ep[11]+"\t"+word_ep[12]+"\t"+word_ep[13]+"\n"
        if refseq_id in refseq_dict:
            out.write(refseq_dict[refseq_id]+"\t"+data)


######### Mapping WES to RNA

# x=1
# for line_rna in open(rna_file,'r'):
#     x=x+1
#     if(x==2):
#         word_rna=line_rna.strip().split("\t")
#         for i in range((len(word_rna[3:-7]))):
#             na.append("NA")
#     if(x>=1):
#         word_rna=line_rna.strip().split("\t")
#         var_id_tmp=word_rna[-1].strip().split(" ")
#         var_id_rna=var_id_tmp[0].strip().split("-")[0]
#         all_rna="\t".join(word_rna[3:-7])
#         length_rna_dat=len(word_rna)
#         dict_rna[var_id_rna]=all_rna

# for line_wc in open(wes_out_file,'r'):
#     word_wc=line_wc.strip().split("\t")
#     var_id_wc=word_wc[0]
#     var_id_wc_tmp=var_id_wc.strip().split(":")
#     var_id_wc_chr=var_id_wc_tmp[0]
#     var_id_wc_bp=var_id_wc_tmp[1].strip().split("-")[0]
#     var_id_wc=var_id_wc_chr+":"+var_id_wc_bp
#     all_wc="\t".join(word_wc[0:])
#     if var_id_wc in dict_rna:
#         wes_rna_out.write(all_wc+"\t"+dict_rna[var_id_wc]+"\n")
#     if not var_id_wc in dict_rna:
#         nan="\t".join(na)
#         wes_rna_out.write(all_wc+"\t"+nan+"\n")
# wes_rna_out.close()

# ###Mapping WES_RNA_count file to FPKM ### 
# for line_fpkm in open(fpkm_file,'r'):
#     word_fpkm=line_fpkm.strip().split("\t")
#     transcript_id_fpkm=word_fpkm[0]
#     all_fpkm="\t".join(word_fpkm[0:])
#     dict_fpkm[transcript_id_fpkm]=all_fpkm

# for line_wrc in open (wes_rna_out_file,'r'):
#     word_wrc=line_wrc.strip().split("\t")
#     transcript_id_wrc=word_wrc[9]
# #    print(transcript_id_wrc)
#     all_wrc="\t".join(word_wrc[0:])
#     if transcript_id_wrc in dict_fpkm:
#         wes_rna_fpkm_out.write(dict_fpkm[transcript_id_wrc]+"\t"+all_wrc+"\n")
#     if transcript_id_wrc not in dict_fpkm:
#         wes_rna_fpkm_out.write("NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t"+all_wrc+"\n")
