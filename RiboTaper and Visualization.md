### RiboTaper analysis

RiboTaper runs on Linux, it requires:  
Samtools(0.1.19, must be in the PATH)
Bedtools (v2.17.0)  
**R (>3.0.1 and <4.0)** with 
seqinr_3.1-3  
ade4_1.7-2  
multitaper_1.0-11  
doMC_1.3.3  
iterators_1.0.7  
foreach_1.4.2  
XNomial_1.0.1  
[See Ohler Lab website for detail](https://ohlerlab.mdc-berlin.de/software/RiboTaper_126/) 

**1. Run RiboTaper Annotation**
```
#$RiboTaper_Path is the path to RiboTaper code files
#$GTF path to the GTF files
#$FASTA path to the FASTA files
#$OUTPUT path to the output folder
#$BED path to bedtools v2.17.0

$RiboTaper_Path/create_annotations_files.bash $GTF $FASTA false false $OUTPUT $BED $Taper
```

**2. Run RiboTaper**
```
#$RiboTaper_Path/Ribotaper.sh path to RiboTaper code files
#$RIBO/ribo.bam path to Ribo-seq bam file (STAR output)
#$RNA/RNA.bam path to RNA-seq bam file(STAR output)
#$ANNO path to RiboTaper annotation files
#$BED path to bedtools v2.17.0
#8 is the number of threads used
#24,25,26,27,28 are read length
#8,9,10,11,12 are cut-off and offset (from Ribo-seQC)

$RiboTaper_Path/Ribotaper.sh $RIBO/ribo.bam $RNA/RNA.bam $ANNO 24,25,26,27,28 8,9,10,11,12 $TAPER $BED 8
```
### RiboTaper ORFs_max_filt file output column information:
```
gene_id<-gene id based on the annotation used
gene_symbol<-gene name based on the annotation used
transcript_id<-transcript id based on the annotation used
annotation<-gene biotype based on the annotation used
length<-length of the transcript
strand<-gene in the forward or reverse orientation
n_exons<-n of exons in the transcript
P_sites_sum<-Total sum of P-sites position mapped to the transcript
RNA_sites<-Total sum of RNA-sites position mapped to the transcript
Ribo_cov_aver<-Average coverage by all Ribo-seq reads on the transcript
RNA_cov_aver<-Average coverage by all RNA-seq reads on the transcript
category<-ORF category
ORF_id_tr<-ORF id with transcript coordinates
start_pos<-ORF start position (transcript coordinates)
stop_pos<-ORF stop position (transcript coordinates)
annotated_start<-Annotated ORF start position (transcript coordinates)
annotated_stop<-Annotated ORF stop position (transcript coordinates)
ORF_id_gen<-ORF id with genomic coordinates
ORF_length<-ORF length
ORF_P_sites<-Total sum of P-sites position mapped to the ORF
ORF_Psit_pct_in_frame<-Percentage of P-sites position mapped in frame
ORF_RNA_sites<-Total sum of RNA-sites position mapped to the ORF
ORF_RNAsit_pct_in_frame<-Percentage of RNA-sites position mapped in frame
ORF_pval_multi_ribo<-P-value for the multitaper test on the ORF region using P-sites positions
ORF_pval_multi_rna<-P-value for the multitaper test on the ORF region using RNA-sites positions
ORF_id_tr_annotated<-ORF id with annotated transcript coordinates
n_exons_ORF<-n of exons in the ORF
pct_covered_onlymulti_ribo<-percentage of the ORF region only supported by multimapping Ribo-seq reads
pct_covered_onlymulti_rna<-percentage of the ORF region only supported by multimapping RNA-seq reads
header_tofasta<-header used in the generated protein fasta file
ORF_pept<-peptide sequence of the identified ORF
```
### Step 9: RiboPlotR visualization
For detail information and examples about RiboPlotR please visit: https://github.com/hsinyenwu/RiboPlotR

### Step 10: GWIPS website and TAIR website
[GWIPS](https://gwips.ucc.ie/cgi-bin/hgGateway) 3 Arabidopsis and 1 Maize datasets   

[TAIR](https://www.arabidopsis.org) 1 very deep Arabidopsis seedling dataset
1. Inside a gene page, click fill screen view:
<img width="1056" alt="image" src="https://github.com/kaufm202/Riboseq_ASPB_2024/assets/4383665/07d99a97-21af-485d-a147-66e2c21bad56">
2. Select tracks 
<img width="1045" alt="image" src="https://github.com/kaufm202/Riboseq_ASPB_2024/assets/4383665/8fd62107-b86e-429f-be54-9b0b96b87133">
 
