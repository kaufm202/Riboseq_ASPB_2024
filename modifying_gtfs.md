# What is a GFF/GTF?
Ensembl provides in-depth description [here](https://useast.ensembl.org/info/website/upload/gff.html).
GFF stands for "General Feature Format". GFF version 2 is referred to as GTF. GFF version 3 is referred to as GFF3 or sometimes just GFF. The two versions only differ in how the attribute column us formatted and generally GTF uses gene_id and transcript_id, while GFF3 uses Parent, Name, ID.  

Many tools will only accept GTF format, so <ins>__it is generally better to download the GTF format if possible.__</ins>

GFF/GTF are tab-delimited files with 1 row per feature and 9 columns:
1. seqid - name of chromosome/scaffold that feature is on
2. source - data source (program, database, or project name)
3. type - type of feature (i.e. gene, mRNA, exon)
4. start - start position of feature
5. end - end position of feature
6. score - floating point value, generally not useful
7. strand - forward (+) or reverse (-) strand
8. phase - indicates phase of first base for coding sequence
9. attributes - list of values with additional information about feature  

The beginning of a GFF/GTF may also contain header rows that are commented out with '##'
```
### Example gtf

$ head arabidopsis_araport11_tx.gtf
##gff-version 2
##source-version rtracklayer 1.58.0
##date 2024-05-07
Chr1    Araport11       gene    3631    5899    .       +       .       gene_id "AT1G01010"; transcript_id "AT1G01010"; ID "AT1G01010"; Name "AT1G01010"
Chr1    Araport11       transcript      3631    5899    .       +       .       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; ID "AT1G01010.1"; Parent "AT1G01010"; Name "AT1G01010.1"
Chr1    Araport11       CDS     3760    3913    .       +       0       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; ID "AT1G01010:CDS:1"; Parent "AT1G01010.1"; Name "NAC001:CDS:1"
Chr1    Araport11       CDS     3996    4276    .       +       2       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; ID "AT1G01010:CDS:2"; Parent "AT1G01010.1"; Name "NAC001:CDS:2"
Chr1    Araport11       CDS     4486    4605    .       +       0       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; ID "AT1G01010:CDS:3"; Parent "AT1G01010.1"; Name "NAC001:CDS:3"
Chr1    Araport11       CDS     4706    5095    .       +       0       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; ID "AT1G01010:CDS:4"; Parent "AT1G01010.1"; Name "NAC001:CDS:4"
Chr1    Araport11       CDS     5174    5326    .       +       0       gene_id "AT1G01010"; transcript_id "AT1G01010.1"; ID "AT1G01010:CDS:5"; Parent "AT1G01010.1"; Name "NAC001:CDS:5"
```

# GFF/GTF errors
One issue with the GFF/GTF format is that different programs expect different formatting (some can only use GTF), require different attributes to be present, or cannot handle certain types of transcripts or features. This leads to frequent and sometimes quite cryptic error messages. Here I will provide some tools that have helped me to avoid and overcome these problems.

# Converting and filtering GFF/GTF
## Command line tools
There are multiple command line tools like gffread and AGAT that can interconvert GFF and GTF. I have found that different tools work for different applications, but have not found a tool that works well for every species and for every program, and these tools typically do not allow much customization.  
See summary of command line tools [here](https://agat.readthedocs.io/en/latest/gff_to_gtf.html)

## rtracklayer
I currently use rtracklayer for all GFF/GTF modifications/conversions since it allows much more customization and allows you to manipulate the GTF like a dataframe in R.  
rtracklayer: [GitHub](https://github.com/lawremi/rtracklayer), [Documentation](https://bioconductor.org/packages/release/bioc/vignettes/rtracklayer/inst/doc/rtracklayer.pdf)  

Below, I will show a few methods I use to resolve various errors that I have encountered while analyzing different species. Keep in mind that <ins>some annotations may not require all of these modifications, and some may require additional modifications.</ins> The nice thing about rtracklayer is that you can make precise modifications as needed when you get errors or notice problems.

Import GFF/GTF into rtracklayer in R:
```
library(rtracklayer)
library(GenomicFeatures)

### Importing GFF3/GTF
in_file <- "path/to/gff3_file.gff3
in_gtf <- rtracklayer::import(in_file, format="gff3")
new_gtf <- in_gtf
```
Once imported, the gtf is now a [GRanges](https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomicRanges/html/GRanges-class.html) object, and can manipulated much like an R dataframe.

```
> head(in_gtf, n=3)

GRanges object with 3 ranges and 19 metadata columns:
      seqnames    ranges strand |    source     type     score     phase
         <Rle> <IRanges>  <Rle> |  <factor> <factor> <numeric> <integer>
  [1]     Chr1 3631-5899      + | Araport11     gene        NA      <NA>
  [2]     Chr1 3631-5899      + | Araport11     mRNA        NA      <NA>
  [3]     Chr1 3760-3913      + | Araport11     CDS         NA         0
                   ID         Name                   Note      symbol
          <character>  <character>        <CharacterList> <character>
  [1]       AT1G01010    AT1G01010 NAC domain containin..      NAC001
  [2]     AT1G01010.1  AT1G01010.1 NAC domain containin..      NAC001
  [3] AT1G01010:CDS:1 NAC001:CDS:1 NAC domain containin..        <NA>
                   full_name computational_description           locus
             <CharacterList>           <CharacterList> <CharacterList>
  [1] NAC domain containin..    NAC domain containin..         2200935
  [2] NAC domain containin..    NAC domain containin..                
  [3]                           NAC domain containin..                
          locus_type          Parent  conf_class conf_rating
         <character> <CharacterList> <character> <character>
  [1] protein_coding                        <NA>        <NA>
  [2]           <NA>       AT1G01010           2        ****
  [3]           <NA>     AT1G01010.1        <NA>        <NA>
                        gene        curator_summary Derives_from
             <CharacterList>        <CharacterList>  <character>
  [1]                                                       <NA>
  [2] 2200934,UniProt=Q0WV96 Member of the NAC do..         <NA>
  [3]                        Member of the NAC do..         <NA>
               Dbxref
      <CharacterList>
  [1]                
  [2]                
  [3]                
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

### Remove unneeded feature types
GFFs often contain features that may not be necessary for your application
```
> table(in_gtf$type)

                    gene                     mRNA                      CDS 
                      35                       77                      367 
                    exon           five_prime_UTR                  protein 
                     418                      116                       77 
         three_prime_UTR                    miRNA miRNA_primary_transcript 
                      80                        3                        2 
                    tRNA 
                       1
```
For our analysis, we may be only interested in mRNA transcripts and their coding sequences:
```
keep_types <- c("gene", "mRNA", "exon", "CDS")
new_gtf <- new_gtf[new_gtf$type %in% keep_types,]
```
And with that filter, we remove many lines that may have caused errors:
```
> message("Removed ", length(in_gtf) - length(new_gtf), " out of ", length(in_gtf), " lines from GFF")
Removed 279 out of 1176 lines from GFF
```
### Removing unneeded attributes
In addition to unneeded features, there will likely be many attributes that you do not need that are in the GRanges object as metadata columns:
```
> colnames(mcols(in_gtf))
 [1] "source"                    "type"                     
 [3] "score"                     "phase"                    
 [5] "ID"                        "Name"                     
 [7] "Note"                      "symbol"                   
 [9] "full_name"                 "computational_description"
[11] "locus"                     "locus_type"               
[13] "Parent"                    "conf_class"               
[15] "conf_rating"               "gene"                     
[17] "curator_summary"           "Derives_from"             
[19] "Dbxref"
```
Sometimes metadata columns contain characters or strings that cause errors, so can remove most of those attributes by specifying the ones we do need. 
```
### Be sure to always keep source, type, score, and phase, since they are required in GTF
keep_cols <- c("source", "type", "score", "phase", "ID", "Parent", "Name", "locus_type")
new_gtf <- new_gtf[,keep_cols]
new_gtf$Parent <- as.character(new_gtf$Parent)
```
This now greatly simplifies the GTF:
```
> colnames(mcols(new_gtf))
[1] "source"     "type"       "score"      "phase"      "ID"        
[6] "Parent"     "Name"       "locus_type"

> head(new_gtf, n=3)
GRanges object with 3 ranges and 8 metadata columns:
      seqnames    ranges strand |    source     type     score     phase
         <Rle> <IRanges>  <Rle> |  <factor> <factor> <numeric> <integer>
  [1]     Chr1 3631-5899      + | Araport11     gene        NA      <NA>
  [2]     Chr1 3631-5899      + | Araport11     mRNA        NA      <NA>
  [3]     Chr1 3760-3913      + | Araport11     CDS         NA         0
                   ID          Parent         Name     locus_type
          <character> <CharacterList>  <character>    <character>
  [1]       AT1G01010                    AT1G01010 protein_coding
  [2]     AT1G01010.1       AT1G01010  AT1G01010.1           <NA>
  [3] AT1G01010:CDS:1     AT1G01010.1 NAC001:CDS:1           <NA>
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
We also need to fix a couple other formatting issues before continuing:
```
### Some attributes are imported as "CharacterList" instead of "character", so we can fix that like so:
new_gtf$Parent <- as.character(new_gtf$Parent)

### GFFs will either use "mRNA" or "transcript" as the type for mRNA transcripts.
### "transcript" is generally the preferred type for GTFs, so we can convert "mRNA" to "transcript"
new_gtf$type <- gsub("mRNA", "transcript", new_gtf$type)

```
### Adding transcript_id and gene_id
Some GFF files do not contain the attributes "gene_id" and "transcript_id", but many programs require these attributes for the GTF format, so we must add it manually:
```
### Start by inserting transcript_id and gene_id as empty columns
### for some programs, gene_id and transcript_id must be in the 5th and 6th position
mcols(new_gtf) <- cbind(mcols(new_gtf)[,1:4], gene_id = "", transcript_id = "", 
                          mcols(new_gtf)[5:ncol(mcols(new_gtf))])

### Fill in the empty columns for each feature type from Name, ID, and Parent columns
# May require some modification depending on whether the gene_id or transcript_id needs to come from Name or ID
new_gtf[new_gtf$type == "gene"]$transcript_id <- new_gtf[new_gtf$type == "gene"]$Name
new_gtf[new_gtf$type == "gene"]$gene_id <- new_gtf[new_gtf$type == "gene"]$Name

new_gtf[new_gtf$type == "transcript"]$transcript_id <- new_gtf[new_gtf$type == "transcript"]$Name
new_gtf[new_gtf$type == "transcript"]$gene_id <- new_gtf$Name[match(new_gtf[new_gtf$type == "transcript"]$Parent, new_gtf$ID)]

new_gtf[new_gtf$type == "exon"]$transcript_id <- new_gtf$Name[match(new_gtf[new_gtf$type == "exon"]$Parent, new_gtf$ID)]
new_gtf[new_gtf$type == "exon"]$gene_id <- new_gtf$gene_id[match(new_gtf[new_gtf$type == "exon"]$Parent, new_gtf$ID)]

new_gtf[new_gtf$type == "CDS"]$transcript_id <- new_gtf$Name[match(new_gtf[new_gtf$type == "CDS"]$Parent, new_gtf$ID)]
new_gtf[new_gtf$type == "CDS"]$gene_id <- new_gtf$gene_id[match(new_gtf[new_gtf$type == "CDS"]$Parent, new_gtf$ID)]
```
Once finished, every feature should have gene_id and transcript_id
```
> head(new_gtf)

GRanges object with 6 ranges and 10 metadata columns:
      seqnames    ranges strand |    source        type     score     phase
         <Rle> <IRanges>  <Rle> |  <factor> <character> <numeric> <integer>
  [1]     Chr1 3631-5899      + | Araport11        gene        NA      <NA>
  [2]     Chr1 3631-5899      + | Araport11  transcript        NA      <NA>
  [3]     Chr1 3760-3913      + | Araport11         CDS        NA         0
  [4]     Chr1 3996-4276      + | Araport11         CDS        NA         2
  [5]     Chr1 4486-4605      + | Araport11         CDS        NA         0
  [6]     Chr1 4706-5095      + | Araport11         CDS        NA         0
          gene_id transcript_id              ID      Parent         Name
      <character>   <character>     <character> <character>  <character>
  [1]   AT1G01010     AT1G01010       AT1G01010        <NA>    AT1G01010
  [2]   AT1G01010   AT1G01010.1     AT1G01010.1   AT1G01010  AT1G01010.1
  [3]   AT1G01010   AT1G01010.1 AT1G01010:CDS:1 AT1G01010.1 NAC001:CDS:1
  [4]   AT1G01010   AT1G01010.1 AT1G01010:CDS:2 AT1G01010.1 NAC001:CDS:2
  [5]   AT1G01010   AT1G01010.1 AT1G01010:CDS:3 AT1G01010.1 NAC001:CDS:3
  [6]   AT1G01010   AT1G01010.1 AT1G01010:CDS:4 AT1G01010.1 NAC001:CDS:4
          locus_type
         <character>
  [1] protein_coding
  [2]           <NA>
  [3]           <NA>
  [4]           <NA>
  [5]           <NA>
  [6]           <NA>
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

### Filter to transcripts that GenomicFeatures thinks are good
GenomicFeatures does a good job of filtering out transcripts that have various issues like overlapping exons, transcripts with exons from different strands, transcripts with multiple CDS, etc. If we only keep genes/transcripts that GenomicFeatures keeps when making a TxDb, then we are less likely to encounter errors down the line.

```
txdb <- GenomicFeatures::makeTxDbFromGFF(in_file, format="gff3")
# you may receive some warning messages telling you which transcripts were removed and why

gn_in_txdb <- unique(na.omit(transcriptLengths(txdb)$gene_id))  # get genes that are in TxDb
tx_in_txdb <- unique(na.omit(transcriptLengths(txdb)$tx_name))  # get transcripts that are in TxDb
tx_in_txdb <- append(tx_in_txdb, gn_in_txdb)  # need to include the genes so that they don't get removed in the transcript filtering

new_gtf <- new_gtf[new_gtf$gene_id %in% gn_in_txdb]
new_gtf <- new_gtf[new_gtf$transcript_id %in% tx_in_txdb]
```

### Remove duplicated transcripts
Remove any transcripts where there are different transcripts with the same name:
```
tx_gtf <- new_gtf[new_gtf$type == "transcript"]
tx_dups <- unique(tx_gtf$transcript_id[duplicated(tx_gtf$transcript_id)])
new_gtf <- new_gtf[!new_gtf$transcript_id %in% tx_dups]
```
### Remove individual problematic transcripts
It is not always clear why a certain transcript causes an error, but sometimes it is easier to simply remove the causal transcript:
```
bad_tx <- c("AT1G01305.1", "AT1G01240.5")  # these transcripts are fine, this is just an example to show this method
new_gtf <- new_gtf[!new_gtf$transcript_id %in% bad_tx]
```
### Remove genes that have no transcripts
Some programs expect at least one transcript from genes in the GFF, so we must remove any genes that do not have an associated transcript:
```
tx_genes <- unique(new_gtf[new_gtf$type == "transcript"]$gene_id)
new_gtf <- new_gtf[new_gtf$gene_id %in% tx_genes,]
```
### CDS only GTF for RSEM
If using Ribo-seq and RNA-seq to quantify translation with [RSEM](https://github.com/deweylab/RSEM), then you will need to align to only the coding sequences, which can be done by keeping only genes, mRNA, and CDS, then converting CDS to exons:
```
keep_types <- c("gene", "mRNA", "transcript", "CDS")
new_gtf <- in_gtf[in_gtf$type %in% keep_types,]

new_gtf$type <- gsub("CDS", "exon", new_gtf$type)
```
### Export to GTF
Once all necessary modifications are made, export modified GFF to GTF format
```
out_file <- "path/to/gtf_output.gtf"
rtracklayer::export(new_gtf, out_file, format="gtf")
```
# Try a different source
If nothing seems to fix your GFF/GTF problems, try downloading the annotation from a different source. Many genomes and annotations are available through multiple sources that have slight differences in the formatting and features they include. For example, Arabidopsis thaliana Araport 11 is available through [Phytozome](https://phytozome-next.jgi.doe.gov/), [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/), [EnsemblPlants](https://plants.ensembl.org/index.html), or [TAIR](https://www.arabidopsis.org/).