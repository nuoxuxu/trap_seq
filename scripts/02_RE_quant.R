library(DESeq2)
library(preprocessCore)

dds <- readRDS(snakemake@input[[1]])
dds <- readRDS("proc/dds.rds")

up_vec <- counts(dds, normalized=TRUE)[, colData(dds)$assay == "IP"]
down_vec <- counts(dds, normalized=TRUE)[, colData(dds)$assay == "Input"]
RE <- up_vec / down_vec
RE_quant_normalized <- preprocessCore::normalize.quantiles(RE, keep.names=TRUE)

if (!dir.exists("results/RE_quant")) {
  dir.create("results/RE_quant", recursive = TRUE)
}

write.csv(RE, snakemake@output[["RE_unormalized"]])
write.csv(RE_quant_normalized, snakemake@output[["RE_normalized"]])
