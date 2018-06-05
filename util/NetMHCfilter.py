#!/usr/bin/python3.5

# Argument 1 -> Epitope file from TN pipeline op
# Argument 2 -> NetMHC binding raw output
# Argument 3 -> Output file

import sys
import argparse
import os
import re
import json
from itertools import chain
from collections import defaultdict

# Positional arguments for filename and path
epitope_file=sys.argv[1]
mhc_file=sys.argv[2]
out_file=sys.argv[3]

# Initializing variables
list1=[]
dict_snp={}
dict_clon={}
dict_tmp={}

#Open epitope file to read and output file to write
epi_f=open(epitope_file,'r')
out_f=open(out_file,'w')
out_f.write("FLAG,POS,HLA,Peptide,Core,Of,Gp,Gl,Tp,Il,Icore,Identity,Score,Aff(nM),%Rank,BindLevel\n")
x=0

#
for line_epi in epi_f:
    col_e=line_epi.strip().split(",")
    snp_id=col_e[1]+","+col_e[8]
    mut_type=col_e[7]
    gene_e=col_e[8]
    aa_change=col_e[12]
    aa_pos=col_e[13]
    rawpeptide=col_e[15]
    mutpeptide=col_e[16]
    #create a dictionary with snp_id as key and gene name, exonic function, aa change, raw peptide seq, mutated peptide seq and aa_pos as value
    dict_snp[snp_id]=str(gene_e+"\t"+mut_type+"\t"+aa_change+"\t"+rawpeptide+"\t"+mutpeptide+"\t"+aa_pos)
#Open MHC file  to read
mhc_f=open(mhc_file,'r')
for line_n in mhc_f:
    # Ignore all empty lines
    if not line_n.strip():
        continue
    # Ignore any line with metadata
    if ("#" in line_n) or  ("---" in line_n) or ("Protein" in line_n) or ("Link" in line_n) or ("Offset" in line_n) or ("NetMHCpan" in line_n) or ("Distance" in line_n) or ("Allele" in line_n) or ("University" in line_n) or ("Pos" in line_n) or ("output" in line_n):
        continue
    # If line has some data then do the following
    else:
        col_n=line_n.strip().split()
        pos=int(col_n[0])
        peptide=col_n[2]
        rank=col_n[13]
        gene_n=col_n[10]
        # if pos =0/1 and lenght of peptide <12
        if(pos==0 or pos==1 and len(peptide)<12):
            for key,value in dict_snp.items():
                ### Check if that gene in NetMHC record is present in epitope file
                if gene_n in value:
                    aa=value.strip().split("\t")[2]
                    aa_pos=int(value.strip().split("\t")[5])
                    exonic_func=value.strip().split("\t")[1]
                    mut_pep=value.strip().split("\t")[4]
                    # if gene is present , check on exonic function and find the position of peptide change
                    if exonic_func == 'nonsynonymous SNV' :
                        aa_position= int(aa[(aa.find('.')+2):len(aa)-1])
                        # check if first 4 AA seq in the NetMHC record is the same as the one in Epitope file and position of peptide change in genome is in <=12 
                        if peptide[:4] in mut_pep and aa_position<=12:
                            out_f.write("MUT AT START OF GENE,"+(",".join(col_n[0:]))+"\n")
                    if exonic_func == 'frameshift deletion' or exonic_func == 'frameshift insertion':
                        aa_position= int(aa[(aa.find('.')+2):len(aa)-2])
                        if(peptide[:4] in mut_pep and aa_position<=12):
                            out_f.write("MUT AT START OF GENE,"+(",".join(col_n[0:]))+"\n")
                    if exonic_func == 'nonframeshift deletion' or exonic_func == 'nonframeshift insertion':
                        if(peptide[:4] in mut_pep and aa_pos<=12):
                            out_f.write("MUT AT START OF GENE,"+(",".join(col_n[0:]))+"\n")
        #If pos=2, peptide length < 11 then ignore since it will not have mutated AA unless it happens in start of the gene/protein
        if(pos==2 and len(peptide)<11):
            for key,value in dict_snp.items():
                if gene_n in value:
                    aa=value.strip().split("\t")[2]
                    aa_pos=int(value.strip().split("\t")[5])
                    exonic_func=value.strip().split("\t")[1]
                    mut_pep=value.strip().split("\t")[4]
                    if peptide[:4] in mut_pep and aa_pos<=12:
                        out_f.write("MUT AT START OF GENE,"+(",".join(col_n[0:]))+"\n")
        #If pos>3 and pos < 10 then write since it will not have mutated AA except start of protein
        if(pos==3 and len(peptide)<10):
            for key,value in dict_snp.items():
                if gene_n in value:
                    aa=value.strip().split("\t")[2]
                    aa_pos=int(value.strip().split("\t")[5])
                    exonic_func=value.strip().split("\t")[1]
                    mut_pep=value.strip().split("\t")[4]
                    if (peptide[:4] in mut_pep and aa_pos<=12):
                        out_f.write("MUT AT START OF GENE,"+(",".join(col_n[0:]))+"\n")
        if(pos==4 and len(peptide)<9):
            for key,value in dict_snp.items():
                if gene_n in value:
                    aa=value.strip().split("\t")[2]
                    aa_pos=int(value.strip().split("\t")[5])
                    exonic_func=value.strip().split("\t")[1]
                    mut_pep=value.strip().split("\t")[4]
                    if peptide[:4] in mut_pep and aa_pos<=12:
                        out_f.write("MUT AT START OF GENE,"+(",".join(col_n[0:]))+"\n")
        if(pos==5 and len(peptide)<8):
            for key,value in dict_snp.items():
                if gene_n in value:
                    aa=value.strip().split("\t")[2]
                    aa_pos=int(value.strip().split("\t")[5])
                    exonic_func=value.strip().split("\t")[1]
                    mut_pep=value.strip().split("\t")[4]
                    if peptide[:4] in mut_pep and aa_pos<=12:
                        out_f.write("MUT AT START OF GENE,"+(",".join(col_n[0:]))+"\n")
        #If pos is greater than 13 (13-real position,12-in file) ignore since it will not have mutated AA except INDELS
        if(pos>12): ##ADD CONDITION SNV TO REMOVE
            for key,value in dict_snp.items():
                ### Check if that gene in NetMHC record is present in epitope file
                if gene_n in value:
                    aa=value.strip().split("\t")[2]
                    aa_pos=int(value.strip().split("\t")[5])
                    exonic_func=value.strip().split("\t")[1]
                    mut_pep=value.strip().split("\t")[4]
                    # if gene is present , check on exonic function and write only if its insertion or deletion
                    if exonic_func == 'nonsynonymous SNV' :
                        if(peptide[:4] in mut_pep):
                            continue
                    if exonic_func == 'frameshift deletion' or exonic_func == 'frameshift insertion':
                        if(peptide[:4] in mut_pep):
                            out_f.write(exonic_func+","+(",".join(col_n[0:]))+"\n")
                    if exonic_func == 'nonframeshift deletion' or exonic_func == 'nonframeshift insertion': 
                        if(peptide[:4] in mut_pep):
                            out_f.write(exonic_func+","+(",".join(col_n[0:]))+"\n")
        else: #write out everything else is not in the conditions above
             out_f.write("-,"+(",".join(col_n[0:]))+"\n")
