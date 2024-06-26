## Steps for Ribo-seq preprocessing and quality control
1. [Downloading reads from NCBI SRA (sra-tools)](#1-downloading-reads-from-ncbi-sra-sra-tools)  
2. [Assessing read quality (FastQC/MultiQC)](#2-assessing-read-quality-fastqcmultiqc)  
3. [Adapter/quality trimming (cutadapt)](#3-adapterquality-trimming-cutadapt)  
4. [Remove contaminants (bowtie2)](#4-removing-contaminants-bowtie2)  
5. [Align to transcripts and CDS (STAR)](#5-splice-aware-alignment-to-transcripts-and-cds-star)  
6. [Ribo-seq quality + P-sites (Ribo-seQC)](#6-ribo-seq-quality--p-sites-ribo-seqc)  
7. [Quantification (RSEM)](#7-quantification-rsem)  
### 1. Downloading reads from NCBI SRA (sra-tools)
sra-tools: [GitHub](https://github.com/ncbi/sra-tools), [Documentation](https://github.com/ncbi/sra-tools/wiki)  
Find Ribo-seq and RNA-seq samples on [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/)  

```
### prefetch downloads compressed .sra file from NCBI SRA
SRA="SRR3498209"
prefetch ${SRA} -O ${DIRECTORY}
    # -O output directory path

### fasterq-dump decompressed .sra to output .fastq
fasterq-dump -3 ${DIRECTORY}/${SRA}/${SRA}.sra -O ${DIRECTORY}
    # -3 automatically detects if reads are single or paired end
    # -O output directory path

### Repeat for additional Ribo-seq and RNA-seq samples

### (Optional) Compress files to save space
pigz ${DIRECTORY}/*.fastq
```
Output fastq files include the reads and their quality scores
```
$ head SRR3498209.fastq 
@SRR3498209.1 SN638:766:HC5HHBCXX:1:1101:1200:2077 length=50
ATGGNTGATTTAGCTTCCAAGAACAAGGAGATCGGAAGAGCACACGTCTG
+SRR3498209.1 SN638:766:HC5HHBCXX:1:1101:1200:2077 length=50
DDDD#<DGHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIHIIIIIII
@SRR3498209.2 SN638:766:HC5HHBCXX:1:1101:1222:2081 length=50
AGAGNGGGTGAGAGCCCCGTCGTGCCCGGAGATCGGAAGAGCACACGTCT
+SRR3498209.2 SN638:766:HC5HHBCXX:1:1101:1222:2081 length=50
DDDD#<DGHHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

### 2. Assessing read quality (FastQC/MultiQC)
FastQC: [GitHub](https://github.com/s-andrews/FastQC), [Webpage](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/), [Tutorial video](https://www.youtube.com/watch?v=bz93ReOv87Y)  
MultiQC: [GitHub](https://github.com/MultiQC/MultiQC), [Webpage](https://multiqc.info/), [Documentation](https://multiqc.info/docs/)
```
### Run FastQC for each sample
fastqc -o ${DIRECTORY} ${DIRECTORY}/*.fastq.gz
    # -O output directory path
    # use "*" to process all fastq in the directory

### (Optional) Combine multiple FastQC reports with MultiQC
multiqc ${DIRECTORY}
    # automatically searches for fastqc reports in directory and combines into new .html file
```
Output .html files include many quality metrics, including quality scores, adapter contamination, overrepresented sequences

### 3. Adapter/quality trimming (cutadapt)
cutadapt: [GitHub](https://github.com/marcelm/cutadapt), [Documentation](https://cutadapt.readthedocs.io/en/stable/)
```
### Trimming for Ribo-seq reads
cutadapt -m 20 -M 40 -q 5 --discard-untrimmed -a AGATCGGAAGAGCACACGTCT \
  -o ${DIRECTORY}/${SRA}_trim.fq.gz ${DIRECTORY}/${SRA}.fastq.gz
    # -m minimum read length
    # -M maximum read length
    # -q quality cutoff for 3' end of reads
    # --discard-untrimmed removes any reads that do not contain adapter
    # -a adapter sequence
    # -o output file name

### Single end RNA-seq
cutadapt -m 20 -q 5 -a AGATCGGAAGAGCACACGTCT --poly-a \
  -o ${DIRECTORY}/${SRA}_trim.fq.gz ${DIRECTORY}/${SRA}.fastq.gz
    # --poly-a trims poly-A tails

### Paired end RNA-seq
cutadapt -m 20 -q 5 -a CCCTACACGACGCTCTTCCGATCT -A TTCAGACGTGTGCTCTTCCGATCT --poly-a \
  -o ${DIRECTORY}/${SRA}_1_trim.fq.gz -p ${DIRECTORY}/${SRA}_2_trim.fq.gz \
  ${DIRECTORY}/${SRA}_1.fastq.gz ${DIRECTORY}/${SRA}_2.fastq.gz
    # -o output read 1, -p output read 2
    # -a forward adapter, -A reverse adapter

### Repeat FastQC to check quality and adapter sequences
```

### 4. Removing contaminants (bowtie2)
bowtie2: [GitHub](https://github.com/BenLangmead/bowtie2), [Documentation](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)  
This step requires a .fasta file containing the contaminant sequences (rRNA, tRNA, snRNA, snoRNA). These can be found from a species' genome annotation or by searching sequences on [NCBI](https://www.ncbi.nlm.nih.gov/nuccore). Can also include rRNA/tRNA sequences from closely related species if unable to find species-specific sequences.
```
### Build bowtie index 
bowtie2-build -f ${CONTAM_FASTA} ${INDEX_BASE}
    # -f path to fasta containing contamination sequences
    # INDEX_BASE will be the prefix of the outputted index files

### Align to contaminants and remove aligned reads
## Single end Ribo-seq or RNA-seq
bowtie2 -L 20 -x ${INDEX_BASE} -S ${SRA}_alignments.sam \
  ${DIRECTORY}/${SRA}_trim.fq.gz --un-gz ${DIRECTORY}/${SRA}_rrna_trim.fq.gz
    # -L length of seed substring (less than or equal to minimum read length)
    # -x index, same as INDEX_BASE from bowtie2-build
    # -S outputs alignments to SAM file (do not need this file)
    # --un-gz output unaligned reads to fastq

## Paired end RNA-seq
bowtie2 -L 20 -x $rRNA_ind \
  -1 ${DIRECTORY}/${SRA}_1_trim_paired.fq.gz -2 ${DIRECTORY}/${SRA}_2_trim_paired.fq.gz
  -S ${DIRECTORY}/${SRA}_alignments.sam \
  --un-conc-gz ${DIRECTORY}/${SRA}_%_rrna_trim.fq.gz
    # -1 read 1 input fastq, -2 read 2 input fastq
    # --un-conc-gz output unaligned reads to fastq, writes two files with 1 and 2 in place of %

### Repeat FastQC to check overrepresented sequences
```
If there are still major overrepresented sequences, you can try running [BLASTn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) and if it is a noncoding RNA like rRNA or tRNA, add it to the fasta and run again.

### 5. Splice-aware alignment to transcripts and CDS (STAR)
STAR: [GitHub](https://github.com/alexdobin/STAR), [Documentation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)  

If doing both translation efficiency analysis and ORF discovery, you will need to run STAR twice:  
__For Ribo-seQC, make index using normal GTF with full transcripts.  
For RSEM, [create CDS-only GTF](modifying_gtfs.md#cds-only-gtf-for-rsem), and use that to build index.__  

In addition, you will need to make __separate indexes for Ribo-seq versus RNA-seq__  
This means a total of 4 STAR indexes to complete all analyses: _CDS_Ribo, CDS_RNA, Full_Tx_Ribo, Full_Tx_RNA_
```
### Create STAR index
STAR --runMode genomeGenerate \
  --genomeDir ${INDEX_DIR} \
  --genomeFastaFiles ${GENOME_FASTA} \
  --sjdbGTFfile ${GTF} \
  --sjdbOverhang ${OVERHANG}
    # --genomeDir index output directory
    # --genomeFastaFiles path to genome fasta
    # --sjdbGTFfile path to GTF file
    # --sjdbOverhang max read length minus 1 (usually ~35 for Ribo-seq, longer for RNA-seq)
    # Must make separate index for Ribo-seq and RNA-seq since --sjdbOverhang will differ

### Run STAR alignment for single end (paired end add another file after --readFilesIn)
STAR --genomeDir ${INDEX_DIR} \
  --readFilesCommand zcat \
  --readFilesIn ${DIRECTORY}/${SRA}_rrna_trim.fq.gz \
  --alignIntronMax 5000 \
  --alignIntronMin 15 \
  --outFilterMismatchNmax 1 \
  --outFilterMultimapNmax 20 \
  --outFilterType BySJout \
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 2 \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --outFileNamePrefix "${DIRECTORY}/${SRA}_star_"
    # --readFilesIn input files from bowtie2 output, for paired end put both files
    # see STAR manual for details on other options
```
STAR will output two different alignment files: _Aligned.sortedByCoord.out.bam (genomic alignment) and _Aligned.toTranscriptome.out.bam (transcriptome alignment). BAM is a binary version of [Sequence Alignment Map (SAM)](https://samtools.github.io/hts-specs/SAMv1.pdf) files. BAM files cannot be viewed directly, but can be viewed using samtools:
```
$ samtools view SRR3498209_star_ribo_Aligned.sortedByCoord.out.bam | head
SRR3498209.240605       0       Chr1    558     0       14M6S   *       0       0       TATTCTGAAGTTCTGCAACG    DDDDDIIIIIIIIIIIIIII    NH:i:7  HI:i:1  AS:i:13    nM:i:0
SRR3498209.37112121     0       Chr1    558     0       14M6S   *       0       0       TATTCTGAAGTTCTGCAACG    @@@DD?D::?=C?GFEIEEF    NH:i:7  HI:i:1  AS:i:13    nM:i:0
SRR3498209.42992077     0       Chr1    558     0       14M6S   *       0       0       TATTCTGAAGTTCTGCAACG    CCCFFFFFHHHHHJJJJJJJ    NH:i:7  HI:i:1  AS:i:13    nM:i:0
SRR3498209.58479780     0       Chr1    558     0       14M6S   *       0       0       TATTCTGAAGTTCTGCAACG    DDDDDIIHHIIIIIIIIIII    NH:i:7  HI:i:1  AS:i:13    nM:i:0
SRR3498209.112905847    0       Chr1    558     0       14M6S   *       0       0       TATTCTGAAGTTCTGCAACG    DDDDDIIIHIIIIIGIIIHI    NH:i:7  HI:i:1  AS:i:13    nM:i:0
SRR3498209.47658999     0       Chr1    1061    255     7S20M   *       0       0       GCGGGGTCAACTCCCCCCACCTCCCCC     ;<99@(;(8@)).=@(8=(<=1:>(==     NH:i:1     HI:i:1  AS:i:17 nM:i:1
SRR3498209.38461983     0       Chr1    1068    255     7S17M   *       0       0       GTTTGAGCCCCACCTCCCCCCCCC        ;;;(()()<@((((((((;@/0<'        NH:i:1     HI:i:1  AS:i:16 nM:i:0
SRR3498209.48318447     16      Chr1    1081    3       20M6S   *       0       0       CCCCCCCCCCACCACCCAAAAGATGC      70&1&8-(8(2(8=</;)@=)(5;<;      NH:i:2     HI:i:1  AS:i:17 nM:i:1
SRR3498209.30620583     16      Chr1    3609    1       18M6S   *       0       0       ACTTCACTGTCTTCCACCTATCCT        IIIHF<<IIIIIIIIIIIIDDDDD        NH:i:3     HI:i:1  AS:i:15 nM:i:1
SRR3498209.947027       0       Chr1    3856    255     24M     *       0       0       GTTGAAGTAGCCATCAGCGAGGTC        DDDDDIIIIIIIIIIIIIIIIIII        NH:i:1     HI:i:1  AS:i:23 nM:i:0
```

### 6. Ribo-seq quality + P-sites (Ribo-seQC)
Ribo-seQC (R package): [Original GitHub](https://github.com/ohlerlab/RiboseQC), [Documentation](https://github.com/lcalviell/Ribo-seQC/blob/master/RiboseQC-manual.pdf)  

Use "_Aligned.sortedByCoord.out.bam" files from full transcript aligned STAR. Can merge or subset reads using [samtools](https://www.htslib.org/)

```
# in bash
### Merge samples with samtools
samtools merge ${BAM_OUT} ${BAM_IN_LIST}
samtools index ${BAM_OUT}

### Subset aligned reads
samtools view -s .1 -b -o ${SUBSET_BAM_OUT} ${BAM_INPUT}
    # -s subsample fraction (i.e., .1 = 10%)
    # -b output in bam format
    # -o output file name
samtools index ${SUBSET_BAM_OUT}
```
If using R version 4.0 or above, we recommend installing a slightly modified version of Ribo-seQC ([Updated GitHub](https://github.com/hsinyenwu/RiboseQC_R4.2.1))
```
# in R
### Install updated Ribo-seQC
library("devtools")
install_github(repo = "hsinyenwu/RiboseQC_R4.2.1")
library("RiboseQC")

prepare_annotation_files(annotation_directory = OUT_DIR,
                         genome_seq = GENOME_FASTA, gtf_file = GTF, 
                         export_bed_tables_TxDb = F, forge_BSgenome = F, create_TxDb = T,
                         circ_chroms = CIRCULAR_SEQUENCES)
    # annotation_directory path to directory to output annotation files
    # genome_seq path to genome fasta file
    # gtf_file path to GTF
    # circ_chroms list of circular chromosomes
    #     i.e. CIRCULAR_SEQUENCES <- c("ChrC", "ChrM") for Arabidopsis TAIR10
    #     remove circ_chroms if chloroplast/mitochondria sequences are not known

RiboseQC_analysis(annotation_file = RANNOT,
                  bam_files = BAM_IN,
                  report_file = OUT_FILE,
                  write_tmp_files = F)
    # annotation_file path to file generated by prepare_annotation_files ending with "_Rannot"
    # bam_files input bam file
    # report_file output html file prefix
    # If 3-nt is too low, Ribo-seQC will not calculate P-site offsets,
    #     add the option readlength_choice_method = "all" to force P-site calcs
```
The main report from Ribo-seQC will be an HTML file, which also contains the P-site offset values. In addition, Ribo-seQC outputs P-site positions that can be used for RiboPlotR with a bit of reformatting:
```
# in R
in_plus <- "path/to/file_P_sites_plus.bedgraph"
in_minus <- "path/to/file_P_sites_minus.bedgraph"
out_file <- "path/to/desired_out_file.riboplotr"

plus <- read.table(in_plus, sep="\t")
minus <- read.table(in_minus, sep="\t")

plus$strand <- "+"
minus$strand <- "-"

minus <- minus[,c(1, 3, 4, 5)]
plus <- plus[,c(1, 3, 4, 5)]

names <- c("chr", "start", "count", "strand")
colnames(plus) <- names
colnames(minus) <- names

comb <- rbind(plus, minus)
comb <- comb[,c(3, 1, 2, 4)]
write.table(comb, file=out_file, col.names = FALSE,
            row.names = FALSE, quote = FALSE, sep="\t")
```
The output file has four columns: count, chromosome, position, strand. This file can be used as the P-site input for RiboPlotR.
```
$ head mergedBamFile.bam_P_sites.riboplotr
21      Chr1    3741    +
1       Chr1    4744    +
13      Chr1    4769    +
8       Chr1    4838    +
18      Chr1    4854    +
9       Chr1    5080    +
18      Chr1    5240    +
26      Chr1    5261    +
26      Chr1    23585   +
2       Chr1    23609   +
```
To use RiboPlotR, you will also need the GTF and RNA-seq BAM file and index, which you can merge and subset to a manageable file size as described earlier with samtools.

### 7. Quantification (RSEM)
RSEM: [GitHub](https://github.com/deweylab/RSEM), [Documentation](https://www.encodeproject.org/documents/0c78ea4b-9392-421b-a6f3-6c858b6002aa/@@download/attachment/RSEM_Documentation.pdf)  
RSEM uses the "_Aligned.toTranscriptome.out.bam" outputs from CDS-aligned STAR  
Here you again need to make a __separate index for Ribo-seq versus RNA-seq__
```
### Build index
rsem-prepare-reference --gtf ${CDS_GTF} ${GENOME_FASTA} ${INDEX_DIR}
    # RSEM will build off the STAR, so ${INDEX_DIR} MUST be the same one used for STAR
    # Run separately for Ribo-seq and RNA-seq
```
To run the quantification, you will need the read length mean standard deviation for single end reads. You can either estimate from the FastQC reports or calculate the mean and SD for each sample using:
```
zcat ${DIRECTORY}/${SRA}_rrna_trim.fq.gz | awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("total=%d \t avg=%f \t stddev=%f\n",n,m,sq/n-m*m);}' -
```
The mean and SD will be different for Ribo-seq and RNA-seq, but the values should be similar between samples of the same type, so you can use the same values within the sample type.
```
### Single end RNA-seq and Ribo-seq
rsem-calculate-expression \
  --fragment-length-mean ${AVE} \
  --fragment-length-sd ${SD} \
  --bam --no-bam-output \
  --forward-prob=1 \
  --seed-length 21 \
  --alignments ${DIRECTORY}/${SRA}_Aligned.toTranscriptome.out.bam \
  ${STAR_CDS_INDEX_DIR} ${DIRECTORY}/${SRA}_rsem
    # change mean, SD, and index for RNA-seq vs. Ribo-seq
    # change --forward-prob depending on library strandedness

### Paired end RNA-seq
rsem-calculate-expression \
  --paired-end --bam --no-bam-output \
  --forward-prob=0 \
  --alignments ${rna_in}/${sample}_star_rna_Aligned.toTranscriptome.out.bam \
  ${STAR_CDS_INDEX_DIR}/rna ${DIRECTORY}/${SRA}_rsem
    # paired end does not need mean and SD
    # change --forward-prob depending on library strandedness
```
RSEM will output two tab-delimited files with abundances at the gene level and transcript level (.genes.results and .isoforms.results) which can be used for downstream translation efficiency and differential expression analysis.
```
$ head SRR3498209_ribo.isoforms.results
transcript_id   gene_id length  effective_length        expected_count  TPM     FPKM    IsoPct
AT1G01010.1     AT1G01010       1290    1260.51 96.00   2.64    3.00    100.00
AT1G01020.1     AT1G01020       738     708.51  0.00    0.00    0.00    0.00
AT1G01020.2     AT1G01020       576     546.51  0.00    0.00    0.00    0.00
AT1G01020.3     AT1G01020       711     681.51  0.00    0.00    0.00    0.00
AT1G01020.4     AT1G01020       711     681.51  0.00    0.00    0.00    0.00
AT1G01020.5     AT1G01020       597     567.51  84.00   5.12    5.83    100.00
AT1G01020.6     AT1G01020       315     285.51  0.00    0.00    0.00    0.00
AT1G01030.1     AT1G01030       1077    1047.51 0.00    0.00    0.00    0.00
AT1G01030.2     AT1G01030       1008    978.51  72.00   2.55    2.90    100.00
```