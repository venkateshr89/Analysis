Sciclone needs
	 1. vaf.dat file from combined annotated tab file from single sample output (ann.all_transcripts.tab)
	 2. Copy number information from segment-gainloss.transformed.tsv from panel single sample output
	 3. LOH information for the varscan.vcf for corresponding to that tumor sample and cat it together

1. ann.all_transcripts has the following fields
   > CHROMOSOME      POSITION        ID      REFERENCE_ALLELE        ALTERNATE_ALLELES       QUALITY_SCORE   FILTER_LABELS   VARSCAN_ABQ     GT_VARSCAN      GT_GATK GT_SAMTOOLS     GT_FREEBAYES    GT_VARDICT      TOTAL_CALLERS   SNPEFF_ANNOTATED_ALLELE SNPEFF_ANNOTATED_IMPACT SNPEFF_ANNOTATED_GENE   SNPEFF_ANNOTATED_TRANSCRIPT     SNPEFF_ANNOTATED_DNA_CHANGE     SNPEFF_ANNOTATED_PROTEIN_CHANGE ANNOVAR_ANNOTATED_ALLELE        ANNOVAR_ANNOTATED_IMPACT        ANNOVAR_ANNOTATED_GENE  ANNOVAR_ANNOTATED_TRANSCRIPT    ANNOVAR_ANNOTATED_DNA_CHANGE    ANNOVAR_ANNOTATED_PROTEIN_CHANGE        EXAC_FREQ_ESTIMATE      1KG_FREQ_ESTIMATE       COSMIC  TARGET_ALLELE_COUNT     TARGET_ALLELE_FRACTION  A_COUNT A_FRACTION      C_COUNT C_FRACTION      G_COUNT G_FRACTION      T_COUNT T_FRACTION      N_COUNT N_FRACTION      TOTAL_READS     TCGA_RANK  TCGA_COUNT   
   > Use the tabparser.py script to extract specific columns to make vaf.dat needed by sciclone
     ####CHROM   POS     TUM1_REF_READS  TUM1_VAR_READS  TUM1_VAF
     #chr1    10146   23      3       38.0
     #chr1    10146   23      3       38.0
 ##### This is done for all tumor samples

2. CN file 'segment-gainloss.transformed.tsv' has following colums
   > FUBP1   462_313797_8880(FUBP1)_22_1     chr1    78413057        78413377        0.148406        84.8656 0.791023        525     1.108
   > we extract columns chr,start,stop and last field (1.108) to add to CNA file
   > cut -f1-4 cnfile.dat | sort | uniq > final_cn.dat file (input to sciclone does not like multiple CNA reggions) 

3. LOH information is in varscan.vcf as ss=3 (somatic status 3 means LOH variant)
   > Extract LOH by  ``` grep "SS=3" varscan.snp.vcf | awk '{print $1"\t"$2"\t"$2+1}' > exclude1.LOH ```
   > Remember to get LOH for all tumor samples and cat it before you use it for sciclone

Refer to sci.R and sciclone Docker file in repo "pdx_ppmp/pipeline/docker/sciclone" to setup sciclone docker image and run sciclone R script 


ENJOY ANALYSIS !!!!!!
