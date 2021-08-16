import os
import csv
import gzip
import datetime

def current_date():
    return datetime.datetime.now().strftime("%d %B %Y, %H:%M:%S")

def start_log(log_file, rule_name):
    with open(log_file, "w") as f:
        f.write("******\n")
        f.write("start of rule " + rule_name + " : " + current_date() + "\n")

def end_log(log_file, rule_name):
    with open(log_file, "a") as f:
        f.write("\nend of rule " + rule_name + " : " + current_date() + "\n")
        f.write("******\n");

configfile: "config.json"

OUTPUT_DIR = config['output_dir'] if 'output_dir' in config else "."
SAMPLE_NAMES = config['sample_names'] if 'sample_names' in config else ["sample1","sample2","sample3","sample4","sample5","sample6"]
FASTQ_SUFFIX = ".fastq.gz"
FASTQ_DIR = config['fastq_dir'] if 'fastq_dir' in config else "/store/EQUIPES/SSFA/antoine.laine/PBodies_06_02_20/data"
STAR_INDEX = config['star_index'] if 'star_index' in config else "/store/EQUIPES/SSFA/MEMBERS/yunfeng.wang/dkpl-anno/index_star_hg38/star"
DATA = config['data_dir'] if 'data_dir' in config else "/store/EQUIPES/SSFA/antoine.laine/PBodies_06_02_20/STAR_Analysis"

STAR_DIR = OUTPUT_DIR + "/alignments_STAR"
DEPTH_DIR = OUTPUT_DIR + "/depth"
RESULTS_DIR = OUTPUT_DIR + "/results"
LOGS = OUTPUT_DIR + "/logs"

rule all:
    input: RESULTS_DIR + "/Final_Last_Exon.tsv"

rule star_align:
    input:
        fastq = FASTQ_DIR + "/{sample}" + FASTQ_SUFFIX,
        index = STAR_INDEX
    output:
        bam  = STAR_DIR + "/{sample}.bam"
    params:
        output_dir = STAR_DIR + "/{sample}"
    log : LOGS + "/{sample}-star_align.log"
    threads: 4
    run:
        start_log(log[0],"star_align")
        shell("STAR --genomeDir {input.index} --runThreadN 4 --readFilesIn {input.fastq} --readFilesCommand gunzip -c --alignIntronMax 50000 --alignSJoverhangMin 10 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.output_dir}")
        end_log(log[0],"star_align")


rule samtools_index:
    input:
        bam = STAR_DIR + "/{sample}.bam"
    output:
        bai = STAR_DIR + "/{sample}.bam.bai"
    log : LOGS + "/{sample}-bai.log"
    threads: 1
    run:
        start_log(log[0],"bai")
        shell("samtools index {input.bam}")
        end_log(log[0],"bai")



rule samtools_depth:
    input:
        bam = STAR_DIR + "/{sample}.bam",
        bai = STAR_DIR + "/{sample}.bam.bai",
        bed = DATA + "/lastExon_MeanSup10Final.bed"
    output:
        depth = DEPTH_DIR + "/{sample}.depth"
    log : LOGS + "/{sample}-depth.log"
    threads: 1
    run:
        start_log(log[0],"depth")
        shell("samtools depth -b {input.bed} -aa {input.bam} > {output.depth} 2>> {log}")
        end_log(log[0],"depth")


rule join_depth:
    input:
        depth = expand("{depth}/{sample}.depth",depth=DEPTH_DIR,sample=SAMPLE_NAMES)
    output:
        join_depth = DEPTH_DIR + "/Join_Depth1to6"
    log : LOGS + "/join_depth.log"
    threads: 1
    run:
        start_log(log[0],"join_depth")
        shell("paste {input.depth[0]} <(awk '{{print $3}}' {input.depth[1]}) <(awk '{{print $3}}' {input.depth[2]}) <(awk '{{print $3}}' {input.depth[3]}) <(awk '{{print $3}}' {input.depth[4]}) <(awk '{{print $3}}' {input.depth[5]}) > {output.join_depth} 2>> {log}")
        end_log(log[0],"join_depth")

rule mean_depth:
    input:
        join_depth = DEPTH_DIR + "/Join_Depth1to6",
        bed = DATA + "/lastExon_MeanSup10Final.bed"
    output:
        split_depth = DEPTH_DIR + "/split_depth.txt",
        vector_PB = DEPTH_DIR + "/vector_PB.txt",
        vector_CYTO = DEPTH_DIR + "/vector_CYTO.txt"
    log : LOGS + "/mean_depth.log"
    threads:1
    run:
        start_log(log[0],"mean_depth")
        shell("scripts/get_mean_depth_vectors.sh {input.bed} {input.join_depth} {output.split_depth} {output.vector_PB} {output.vector_CYTO}")
        end_log(log[0],"mean_depth")

rule r_analysis:
    input:
        vector_PB = DEPTH_DIR + "/vector_PB.txt",
        vector_CYTO = DEPTH_DIR + "/vector_CYTO.txt",
        complete_bed = DATA + "/last_exon_final.bed",
        fasta = DATA + "/signif_last_exon.fa",
    output:
        final_tsv = RESULTS_DIR + "/Final_Last_Exon.tsv"
    params:
        output_dir = RESULTS_DIR
    log : LOGS + "/rscript.log"
    threads:1
    run:
        start_log(log[0],"rscript")
        shell("Rscript scripts/last-exon2.R --vectorPB {input.vector_PB} --vectorCYTO {input.vector_CYTO} --bed {input.complete_bed} --fasta {input.fasta} --output {params.output_dir} --log {log} 2>> {log}")
        end_log(log[0],"rscript")
