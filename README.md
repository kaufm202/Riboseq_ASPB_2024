## Ribo-seq analysis preprocessing and quality control
<img src="https://github.com/kaufm202/Riboseq_ASPB_2024/assets/113535253/35185277-46a1-4db9-bba4-32a8da6bdf4c" width=65% height=65%>

### Setup and general recommendations
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
### RiboTaper for ORF identification
<img alt="image" src="https://github.com/kaufm202/Riboseq_ASPB_2024/assets/4383665/74942c9d-9734-452a-91d6-d803a9b137bb" width=65% height=65%>  
[Instruction and example code](https://github.com/kaufm202/Riboseq_ASPB_2024/blob/main/RiboTaper%20and%20Visualization.md)  

### RiboPlotR for ORF visualization
<img alt="image" src="https://github.com/kaufm202/Riboseq_ASPB_2024/assets/4383665/a35752ab-1992-46ff-9160-4545be3c97f3" width=65% height=65%>  
[Instruction and example code](https://github.com/hsinyenwu/RiboPlotR)  
[Make gtf files for uORFs](https://github.com/hsinyenwu/RiboPlotR_addition)

