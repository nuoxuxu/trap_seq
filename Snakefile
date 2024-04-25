import os
from pathlib import Path

flist = [path for path in Path("data/WEN26595.20240222/240214_A02072_0139_AH2M75DSXC").iterdir() if path.name.endswith("gz")]
sample_list = list(set([str(path).rsplit("_", 3)[0].split("/")[-1] for path in flist]))

localrules: generate_input_text_file

rule all:
    input:
        "proc/rMATS_results",
        "proc/tmp"

rule star_index:
    input:
        genome_fasta_files=f"{os.environ['GENOMIC_DATA_DIR']}Ensembl/Human/Release_104/Raw/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        annotations_gtf=f"{os.environ['GENOMIC_DATA_DIR']}Ensembl/Human/Release_104/Raw/Homo_sapiens.GRCh38.104.gtf"
    output:
        directory(f"{os.environ['GENOMIC_DATA_DIR']}Ensembl/Human/Release_104/STAR_genomeDir")
    conda: "trap_seq"
    resources:
        runtime="1h"
    shell: "STAR --runThreadN 80 --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.genome_fasta_files} --sjdbGTFfile {input.annotations_gtf} --sjdbOverhang 100"

def get_fq1(wildcards):
    return f"data/WEN26595.20240222/240214_A02072_0139_AH2M75DSXC/{wildcards.sample}_L004_R1_001.fastq.gz"

def get_fq2(wildcards):
    return f"data/WEN26595.20240222/240214_A02072_0139_AH2M75DSXC/{wildcards.sample}_L004_R2_001.fastq.gz"

rule align:
    input:
        index=f"{os.environ['GENOMIC_DATA_DIR']}Ensembl/Human/Release_104/STAR_genomeDir",
        fq1=get_fq1,
        fq2=get_fq2
    params:
        annotations_gtf=f"{os.environ['GENOMIC_DATA_DIR']}Ensembl/Human/Release_104/Raw/Homo_sapiens.GRCh38.104.gtf"
    output:
        directory('proc/STAR_results/{sample}')
    shell: "STAR --runThreadN 80 --genomeDir {input.index} --readFilesIn {input.fq1} {input.fq2} --sjdbGTFfile {params.annotations_gtf} --readFilesCommand zcat --outFileNamePrefix {output}/ --outSAMtype BAM SortedByCoordinate"

IP_list = [s for s in sample_list if "IP" in s]
KO_list = [s for s in IP_list if s.count("_") == 3]
WT_list = [s for s in IP_list if s.count("_") == 2]

rule generate_input_text_file:
    input: 
        B1=expand('/scratch/s/shreejoy/nxu/trap_seq/proc/STAR_results/{sample}/Aligned.sortedByCoord.out.bam', sample=KO_list),
        B2=expand('/scratch/s/shreejoy/nxu/trap_seq/proc/STAR_results/{sample}/Aligned.sortedByCoord.out.bam', sample=WT_list)
    output:
        "proc/B1.txt",
        "proc/B2.txt"
    run:
        with open(output[0], "w") as f:
            f.write(",".join(input.B1))
        with open(output[1], "w") as f:
            f.write(",".join(input.B2))

rule rMATS:
    input:
        gtf=f"{os.environ['GENOMIC_DATA_DIR']}Ensembl/Human/Release_104/Raw/Homo_sapiens.GRCh38.104.gtf",
        B1="proc/B1.txt",
        B2="proc/B2.txt"
    output:
        directory("proc/rMATS_results"),
        directory("proc/tmp")
    conda: "trap_seq"
    resources:
        runtime="5h"
    shell: "rmats.py --gtf {input.gtf} --b1 {input.B1} --b2 {input.B2} --readLength 151 --od {output[0]} --tmp {output[1]} --nthread 80"