# Analysis of Last Exon variation between Pbodies and Cytoplasmic samples

## Requirements

### Packages and dependencies
- Snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- Samtools
- Bedtools

It is advised to install these 3 in a single conda environment.

- R Packages : ggplot2,reshape2,dplyr,biomaRt,changepoint,argparse,seqinr

### Datas
You will need a proper STAR index of Human Genome 38

You also need the fastq files of the 3 Pbodies and 3 Presort samples

All other required datas are available in the "data" folder


### Run
#### Step 1 : Complete the config.json files with following parameters:

- **fastq_dir**: Location of FASTQ files.
- **star_index**: Location of STAR HG38 index.
- **sample_names**: List of sample names. Same name as your fastq files, without extension (.fq.gz). Default : ["sample1","sample2","sample3","sample4","sample5","sample6"]
- **data_dir**: Location of joined "data" folder.
- **output_dir**: Location where results will be written.

#### Step 2 : Run
```
snakemake
```
-It is strongly advise, if possible, to add "--cores 6" to speed up the process.
