### Analyzing translation efficiency with Xtail
[Original paper](https://www.nature.com/articles/ncomms11194)
```
######################
library(xtail)
library(tximport)
#D and A are two different treatments
#RNA
dir <- "/path/to/CDS_quant_RSEM/RNA/"
sampletxt <-data.frame(sample=c("D1","D2","D3","T1","T2","T3")) 
files1 <- file.path(dir,sampletxt$sample,paste0(sampletxt$sample,".genes.results"))
names(files1) <- c(paste0("D_RNA", 1:3),paste0("T_RNA", 1:3))

#ribo
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
test.results <- xtail(mrna,rpf,condition,bins=1000) #bins set to 1000 to allow quick sample run
head(test.results$resultsTable)
# log2 fold changes in two conditions
plotFCs(test.results)
# RPF-to-mRNA ratios in two conditions
plotRs(test.results)

```
