
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

### 5. Align to genome (STAR)
STAR: [GitHub](https://github.com/alexdobin/STAR), [Documentation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)  
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

### Run STAR alignment
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
STAR will output two different alignment files: _Aligned.sortedByCoord.out.bam (genomic alignment) and _Aligned.toTranscriptome.out.bam (transcriptome alignment)

### 6. Ribo-seq quality + P-sites (Ribo-seQC)
Ribo-seQC (R package): [Original GitHub](https://github.com/ohlerlab/RiboseQC), [Documentation](https://github.com/lcalviell/Ribo-seQC/blob/master/RiboseQC-manual.pdf)  

<ins>Use "_Aligned.sortedByCoord.out.bam" files from STAR.</ins> Can merge or subset reads using [samtools](https://www.htslib.org/)

```
# in bash
### Merge samples with samtools
samtools merge ${BAM_OUT} ${BAM_IN_LIST}

### Subset aligned reads
samtools view -s .1 -b -o ${BAM_SUBSET} ${BAM_INPUT}
    # -s subsample fraction (i.e., .1 = 10%)
    # -b output in bam format
    # -o output file name
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
```
The main report from Ribo-seQC will be an HTML file. In addition, Ribo-seQC outputs P-sites that can be used for RiboPlotR with a bit of reformatting:
```
in_plus <- "path/to/file_coverage_plus.bedgraph"
in_minus <- "path/to/file_coverage_minus.bedgraph"
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

