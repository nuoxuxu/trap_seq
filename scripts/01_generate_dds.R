library(tximport)
library(dplyr)
library(DESeq2)
library(GenomicFeatures)
library(stringr)

# get salmon results
txdb <- makeTxDbFromGFF(snakemake@input[["gtf_path"]])
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene <- bind_rows(tx2gene, data.frame(TXNAME = c("EGFP"), GENEID = c("EGFP")))
files <- list.files(snakemake@input[["salmon_results"]], full.names = TRUE)
names(files) <- stringr::str_extract(files, "\\d{1,2}(_\\d)?")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# filter
filtered_count_matrix <- txi$counts[rowSums(txi$counts) != 0, ]
filtered_count_matrix <- filtered_count_matrix[rowSums(filtered_count_matrix) >= 20, ]

# create DESEq2DataSet
samples <- names(files)
assay <- str_extract(files, "IP|Input")
condition <- ifelse(grepl("_", samples), "PTEN_KO", "WT")
coldata <- data.frame(samples = samples, assay = assay, condition = condition)
dds <- DESeqDataSetFromTximport(txi, coldata, design = ~ samples)
dds <- DESeq2::estimateSizeFactors(dds)

# Write DESeq2DataSet to disk
saveRDS(dds, snakemake@output[[1]])