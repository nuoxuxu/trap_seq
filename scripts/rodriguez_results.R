library(stringr)
library(dplyr)

gtf_path <- "/scratch/s/shreejoy/nxu/Genomic_references/hg38/Raw/Homo_sapiens.GRCh38.111.gtf"
annotation_from_gtf <- rtracklayer::import(gtf_path) %>%
    dplyr::as_tibble() %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name, gene_biotype)

RE <- read.csv(snakemake@input[[1]], row.names = 1)
RE_KO <- RE[, str_detect(colnames(RE), "4_7|4_9|95_3")]
RE_WT <- RE[, str_detect(colnames(RE), "87|90|91")]

log2_RE_FC <- log(rowMeans(RE_KO, na.rm = TRUE) / rowMeans(RE_WT, na.rm = TRUE), 2)
Z_score <- rowMeans(RE_KO, na.rm = TRUE) - rowMeans(RE_WT, na.rm = TRUE) / sqrt(apply(RE_KO, 1, sd, na.rm = TRUE)**2 + apply(RE_WT, 1, sd, na.rm = TRUE)**2)

results <- data.frame(gene_id = names(log2_RE_FC), log2_RE_FC = log2_RE_FC, Z_score = Z_score)

# add gene symbols to the data frame
results <- results %>%
    mutate(gene_id = rownames(results)) %>%
    left_join(annotation_from_gtf, join_by(gene_id))

# Selecting genes with significant shifts in RE. This is the same filtering threshold as in Rodriguez et al.
results <- results %>%
    filter_if(is.numeric, all_vars(!is.na(.) & !is.infinite(.))) %>%
    filter(abs(Z_score) > 2 & abs(log2_RE_FC) > 0.2) %>%
    arrange(desc(Z_score))

write.csv(results, snakemake@output[[1]])