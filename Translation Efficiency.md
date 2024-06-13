### Analyzing differential translation with Xtail

Since the differential translation analysis needs to consider the variability for both RNA-seq and Ribo-seq, it is not as accurate as the RNA-seq (differential expression) analysis.

[Original paper](https://www.nature.com/articles/ncomms11194)  
[More examples](https://rdrr.io/github/xryanglab/xtail/f/vignettes/xtail.Rmd)  
Studying the changes in translation efficiency in two conditions. You need to 
1. Make a gtf/gff3 file for only CDS region of all genes.
2. Map the Ribo-seq and RNA-seq reads to the CDS only gtf/gff regions of the genome with STAR.
3. Quantify Ribo-seq and RNA-seq reads with the STAR output bam files with RSEM.
4. For step 2 or 3, you can also use Salmon or Kallisto for quantification. tximport (see below) can also input the output from Salmon or Kallisto.

#### Install Xtail in R.
```
#Need to download gfortran to install Xtail in R4.3 or above.
#For example: in https://cran.r-project.org/bin/macosx/tools/ download gfortran-12.2-universal.pkg
#Install Xtail
library(devtools)
install_github("xryanglab/xtail")
library(xtail)
```

#### Example code:
```
######################
library(xtail)
library(tximport)
#D and T are two different treatments
#RNA (the input files here are RSEM quantification files for each sample (biological replicate).)
dir <- "/path/to/CDS_quant_RSEM/RNA/"
sampletxt <-data.frame(sample=c("D1","D2","D3","T1","T2","T3")) 
files1 <- file.path(dir,sampletxt$sample,paste0(sampletxt$sample,".genes.results"))
names(files1) <- c(paste0("D_RNA", 1:3),paste0("T_RNA", 1:3))

#ribo (same, the input files here are RSEM quantification files for each sample (biological replicate).)
dir2 <- "/path/to/CDS_quant_RSEM/ribo/"
sampletxt2 <-data.frame(sample=c("D1","D2","D3","T1","T2","T3"))
files2 <- file.path(dir2,sampletxt2$sample,paste0(sampletxt2$sample,".genes.results"))
names(files2) <- c(paste0("D_ribo", 1:3),paste0("T_ribo", 1:3))

files <- c(files1,files2)

txi.RSEM <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
counts <- txi.RSEM$counts

#Remove those genes with low RNA-seq read counts
keep <- rowSums(counts[,1:6]) >= 100
counts <- counts[keep,] 

rpf <- round(counts[,7:12],0)
rownames(rpf) <- rownames(counts)
mrna <- round(counts[,1:6],0)
rownames(mrna) <- rownames(counts)
condition <- c("D","D","D","T","T","T")

# https://rdrr.io/github/xryanglab/xtail/f/vignettes/xtail.Rmd
test.results <- xtail(mrna,rpf,condition,bins=1000) #bins set to 1000 to allow quick sample run, use default 10000 for better results
head(test.results$resultsTable)
# log2 fold changes in two conditions
plotFCs(test.results)
# RPF-to-mRNA ratios in two conditions
plotRs(test.results)

#test.results contain all statistical analysis results for each gene
#test.results contain all statistical analysis results for each gene
DF <- as.data.frame(test.results$resultsTable)
DF$Gene <- rownames(test.results$resultsTable)
head(test.results$resultsTable)
head(DF)

############################################
# adj p < 0.01, 1.5 fold induced filtering #
############################################
sig.results2 <- DF %>% filter(pvalue.adjust < 0.01) %>% arrange(pvalue.adjust)
sig.resultsP <- sig.results2 %>% filter(log2FC_TE_v1 > log2(1.5)) %>% arrange(pvalue.adjust)
nrow(sig.resultsP)
head(sig.resultsP)

```
