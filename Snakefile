import os
from pathlib import Path

flist = [path for path in Path("data/240214_A02072_0139_AH2M75DSXC").iterdir() if path.name.endswith("gz")]
sample_list = list(set([str(path).rsplit("_", 3)[0].split("/")[-1] for path in flist]))
localrules: generate_input_text_file, generate_dds, get_interaction_results, create_juncfiles, shinyapp

# rule all:
#     input:
#         "results/RE_quant/rodriguez_results.csv"

# rule all:
#     input:
#         expand("proc/STAR_results/{sample}", sample=sample_list)

# rule all:
#     input:
#         expand("proc/leafcutter/junc_files/{sample}.junc", sample=sample_list)

rule all:
    input:
        expand("proc/leafcutter/ds_results/IP_only/leafcutter_ds{suffix}", suffix=["_cluster_significance.txt", "_effect_sizes.txt", "_results.rds"])

###################################
#           STAR                  #
###################################

rule star_index:
    input:
        genome_fasta_files=f"{os.environ['GENOMIC_DATA_DIR']}Ensembl/Human/Release_104/Raw/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        annotations_gtf=f"{os.environ['GENOMIC_DATA_DIR']}Ensembl/Human/Release_104/Raw/Homo_sapiens.GRCh38.104.gtf"
    output:
        directory(f"{os.environ['GENOMIC_DATA_DIR']}Ensembl/Human/Release_104/STAR_genomeDir")
    conda: "rMATS"
    resources:
        runtime="1h"
    shell: "STAR --runThreadN 80 --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.genome_fasta_files} --sjdbGTFfile {input.annotations_gtf} --sjdbOverhang 100"

def get_fq1(wildcards):
    return f"data/240214_A02072_0139_AH2M75DSXC/{wildcards.sample}_L004_R1_001.fastq.gz"

def get_fq2(wildcards):
    return f"data/240214_A02072_0139_AH2M75DSXC/{wildcards.sample}_L004_R2_001.fastq.gz"

rule align:
    input:
        index=f"{os.environ['GENOMIC_DATA_DIR']}Ensembl/Human/Release_104/STAR_genomeDir",
        fq1=get_fq1,
        fq2=get_fq2
    params:
        annotations_gtf=f"{os.environ['GENOMIC_DATA_DIR']}Ensembl/Human/Release_104/Raw/Homo_sapiens.GRCh38.104.gtf"
    resources:
        runtime="1h"
    output:
        directory('proc/STAR_results/{sample}')
    conda:
        "rMATS"
    shell: "STAR --runThreadN 80 --genomeDir {input.index} --readFilesIn {input.fq1} {input.fq2} --sjdbGTFfile {params.annotations_gtf} --outSAMstrandField intronMotif --readFilesCommand zcat --outFileNamePrefix {output}/ --outSAMtype BAM SortedByCoordinate"

###############################
#        Leafcutter           #
###############################

rule run_regtools:
    input:
        bamfile="proc/STAR_results/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        junc="proc/leafcutter/junc_files/{sample}.junc"
    conda: "trap_seq"
    resources:
        runtime="1h"
    envmodules:
        "CCEnv"
    shell:
        """
        samtools index --threads 80 {input.bamfile}
        regtools-0.4.2 junctions extract -a 8 {input.bamfile} -o {output.junc}
        """

rule create_juncfiles:
    input:
        "proc/leafcutter/junc_files"
    output:
        "proc/leafcutter/juncfiles.txt"
    shell:
        "find {input} -type f > {output}"

rule intron_clustering:
    input:
        juncfiles="proc/leafcutter/juncfiles.txt"
    params:
        prefix="proc/leafcutter/cluster_results/trap_seq"
    output:
        expand("proc/leafcutter/cluster_results/trap_seq{suffix}", suffix=["_perind_numers.counts.gz", "_perind.counts.gz", "_pooled", "_refined", "_sortedlibs"])
    resources:
        runtime="1h"
    conda: "leafcutter"
    shell:
        "python2 /scratch/s/shreejoy/nxu/leafcutter/clustering/leafcutter_cluster_regtools.py -j {input.juncfiles} -m 50 -o {params.prefix} -l 500000 -r proc/leafcutter/run_dir"

# type is "IP_only" or "all"
rule diff_spl:
    input:
        perind_numers="proc/leafcutter/cluster_results/trap_seq_perind_numers.counts.gz",
        groups_file="proc/leafcutter/{type}.txt"
    params:
        prefix="proc/leafcutter/ds_results/{type}/leafcutter_ds",
        min_samples_per_intron=lambda wildcards: {"IP_only": 3, "all": 5}[wildcards.type]
    output:
        expand("proc/leafcutter/ds_results/{type}/leafcutter_ds{suffix}", suffix=["_cluster_significance.txt", "_effect_sizes.txt", "_results.rds"], allow_missing=True)
    resources:
        runtime="1h"
    conda: "leafcutter"
    shell:
        "/scratch/s/shreejoy/nxu/leafcutter/scripts/leafcutter_ds.R --num_threads 80 {input.perind_numers} {input.groups_file} -o {params.prefix} -i {params.min_samples_per_intron}"

rule annotations:
    input: "/scratch/s/shreejoy/nxu/Genomic_references/hg38/Raw/Homo_sapiens.GRCh38.104.gtf"
    output:
        expand("proc/leafcutter/visualization/trap_seq{suffix}", suffix=["_all_exons.txt.gz", "_all_introns.bed.gz", "_fiveprime.bed.gz", "_threeprime.bed.gz"])
    params:
        prefix="proc/leafcutter/visualization/trap_seq"
    shell:
        "/scratch/s/shreejoy/nxu/leafcutter/leafviz/gtf2leafcutter.pl -o {params.prefix} {input}"

rule shinyapp:
    input:
        groups_file="proc/leafcutter/{type}.txt",
        perind_numers="proc/leafcutter/cluster_results/trap_seq_perind_numers.counts.gz",
    params:
        prefix="proc/leafcutter/visualization/trap_seq"
    output:
        "proc/leafcutter/visualization/leafviz_{type}.RData"
    conda: "leafcutter"
    shell:
        "/scratch/s/shreejoy/nxu/leafcutter/leafviz/prepare_results.R -m {input.groups_file} {input.perind_numers} proc/leafcutter/ds_results/{wildcards.type}/leafcutter_ds_cluster_significance.txt proc/leafcutter/ds_results/{wildcards.type}/leafcutter_ds_effect_sizes.txt {params.prefix} -o {output}"

# cd /scratch/s/shreejoy/nxu/leafcutter/leafviz
# ./run_leafviz.R /scratch/s/shreejoy/nxu/trap_seq/proc/leafcutter/visualization/leafviz_IP_only.RData
    
###############################
#             rMATs           #
###############################

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
    conda: "rMATS"
    resources:
        runtime="5h"
    shell: "rmats.py --gtf {input.gtf} --b1 {input.B1} --b2 {input.B2} --readLength 151 --od {output[0]} --tmp {output[1]} --nthread 80"

###########################
#         DESeq2          #
###########################

rule generate_dds:
    input: 
        gtf_path=Path(os.environ["GENOMIC_DATA_DIR"]).joinpath("Ensembl/Human/Release_104/Raw/Homo_sapiens.GRCh38.104.gtf")
    output: 
        "proc/dds_wendy.rds",
        "proc/txi_wendy.rds",
        "proc/dds_nuo.rds",
        "proc/txi_nuo.rds",
    conda: "trap_seq"
    script: "scripts/01_generate_dds.R"

rule get_interaction_results:
    input: "proc/dds_wendy.rds"
    output: 
        "results/interaction/raw_interaction_results.csv",
        "results/interaction/shrink_interaction_results.csv"
    conda: "trap_seq"
    script: "scripts/03_interaction.R"