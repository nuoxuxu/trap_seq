# Preparation

## Download libraries

install.packages("dplyr")
install.packages("stringr")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("preprocessCore", version="1.62")
BiocManager::install("DESeq2")

## Loading libraries

library(dplyr)
library(preprocessCore)
library(stringr)
library(DESeq2)

## Define functions

quantile_normalize <- function(dds) {
    up_vec <- counts(dds, normalized=TRUE)[, colData(dds)$assay == "IP"]
    down_vec <- counts(dds, normalized=TRUE)[, colData(dds)$assay == "Input"]
    RE <- up_vec / down_vec
    normalize.quantiles(RE, keep.names=TRUE)
}

get_results <- function(RE_normalized) {
    RE_KO <- RE_normalized[, str_detect(colnames(RE_normalized), "4_7|4_9|95_3")]
    RE_WT <- RE_normalized[, str_detect(colnames(RE_normalized), "87|90|91")]

    log2_RE_FC <- log(rowMeans(RE_KO, na.rm = TRUE) / rowMeans(RE_WT, na.rm = TRUE), 2)
    Z_score <- rowMeans(RE_KO, na.rm = TRUE) - rowMeans(RE_WT, na.rm = TRUE) / sqrt(apply(RE_KO, 1, sd, na.rm = TRUE)**2 + apply(RE_WT, 1, sd, na.rm = TRUE)**2)

    data.frame(gene_id = names(log2_RE_FC), log2_RE_FC = log2_RE_FC, Z_score = Z_score)
}

get_dim <- function(results) {
    results %>%
        filter_if(is.numeric, all_vars(!is.na(.) & !is.infinite(.))) %>%
        filter(abs(Z_score) > 2 & abs(log2_RE_FC) > 0.2) %>%
        dim()
}

## Load data

dds_nuo <- readRDS("proc/dds_nuo.rds")
dds_wendy <- readRDS("proc/dds_wendy.rds")

# Calculate RE

nuo_RE_normalized <- quantile_normalize(dds_nuo)
nuo_rodriguez_results <- get_results(nuo_RE_normalized)
nuo_rodriguez_results %>% get_dim()