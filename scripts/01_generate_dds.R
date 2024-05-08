library(tximport)
library(DESeq2)
library(GenomicFeatures)
library(stringr)
library(here)
library(tidyverse)

save_txi <- function(files, tx2gene, output_path) {
    txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
    saveRDS(txi, output_path)
}
save_dds <- function(files, tx2gene, output_path) {
    txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
    samples <- factor(names(files), levels = c("4_7", "4_9", "95_3", "87", "90", "91"))
    assay <- factor(str_extract(files, "IP|Input"), levels = c("IP", "Input"))
    condition <- factor(ifelse(grepl("_", samples), "PTEN_KO", "WT"), levels = c("WT", "PTEN_KO"))
    coldata <- data.frame(samples = samples, assay = assay, condition = condition)
    dds <- DESeqDataSetFromTximport(txi, coldata, design = ~ samples)
    dds <- DESeq2::estimateSizeFactors(dds)
    saveRDS(dds, output_path)
}

txdb <- makeTxDbFromGFF(snakemake@input[["gtf_path"]])
k <- keys(txdb, keytype = "TXNAME")

# wendy
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene <- bind_rows(tx2gene, data.frame(TXNAME = c("EGFP"), GENEID = c("EGFP")))
files <- list.files("data/salmon_with_eGFP/compiled_quants_egfp", full.names = TRUE)
names(files) <- stringr::str_extract(files, "\\d{1,2}(_\\d)?")
save_dds(files, tx2gene, snakemake@output[[1]])

files <- list.files("data/salmon_with_eGFP/compiled_quants_egfp", full.names = TRUE)
names(files) <- stringr::str_extract(files, "\\d+(_\\d+)?_(Input|IP)")
save_txi(files, tx2gene, snakemake@output[[2]])

# nuo
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene <- bind_rows(tx2gene, data.frame(TXNAME = c("eGFP"), GENEID = c("eGFP")))
files <- file.path(here("data/salmon_results"), list.files(here("data/salmon_results")), "quant.sf")
names(files) <- str_extract(list.files(here("data/salmon_results")), "(\\d+(_\\d+)?)")
save_dds(files, tx2gene, snakemake@output[[3]])

files <- list.files("data/salmon_with_eGFP/compiled_quants_egfp", full.names = TRUE)
names(files) <- stringr::str_extract(files, "(\\d+(_\\d+)?)_(Input|IP)")
save_txi(files, tx2gene, snakemake@output[[4]])