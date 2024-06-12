# Ribo-seq analysis preprocessing and quality control
## Steps for Ribo-seq preprocessing and quality control
1. [Downloading reads from NCBI SRA (sra-tools)](pipeline#1-downloading-reads-from-ncbi-sra-sra-tools)  
2. [Assessing read quality (FastQC/MultiQC)](pipeline#2-assessing-read-quality-fastqcmultiqc)  
3. [Adapter/quality trimming (cutadapt)](pipeline#3-adapterquality-trimming-cutadapt)  
4. [Remove contaminants (bowtie2)](pipeline#4-removing-contaminants-bowtie2)  
5. [Align to transcripts and CDS (STAR)](pipeline#5-splice-aware-alignment-to-transcripts-and-cds-star)  
6. [Ribo-seq quality + P-sites (Ribo-seQC)](pipeline#6-ribo-seq-quality--p-sites-ribo-seqc)  
7. [Quantification (RSEM)](pipeline#7-quantification-rsem)  
<img src="https://github.com/kaufm202/Riboseq_ASPB_2024/assets/113535253/8aa554d5-6099-4964-8143-4e74bc346fb6" width=50% height=50%>

## Setup and general recommendations
Steps 1-5 and 7 are run with bash scripts using packages available through [miniconda](https://docs.anaconda.com/free/miniconda/)

Step 6 (Ribo-seQC) is run in [R](https://cran.r-project.org/bin/windows/base/)  

It is best to run this pipeline on a computing cluster/server. Many of the steps require a large amount of disk space for reads and alignment, large amounts of RAM for species with large genomes, or take a long time to run on a personal computer.

If running multiple samples, it is more efficient to run scripts as loops in bash:
```
sra_list="SRR3498209 SRR3498210 SRR3498211 SRR3498215 SRR3498216 SRR3498217"
type_list="ribo ribo ribo rna rna rna"

sra_arr=($sra_list)
type_arr=($type_list)

for i in "${!sra_arr[@]}"
do
	sra="${sra_arr[$i]}"
	read_type="${type_arr[$i]}"

	if [ "$read_type" == "ribo" ]; then
        commands for Ribo-seq here
	elif [ "$read_type" == "rna" ]; then
        command for RNA-seq here
    fi
done
```
