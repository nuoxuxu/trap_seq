# Preparation

## Loading libraries

library(dplyr)
library(preprocessCore)
library(stringr)

## Define functions

quantile_normalize <- function(df) {
    df %>%
        data.matrix() %>%
        preprocessCore::normalize.quantiles() %>%
        as.data.frame()
}
get_results <- function(RE_normalized) {
    RE_KO <- RE[, str_detect(colnames(RE), "4_7|4_9|95_3")]
    RE_WT <- RE[, str_detect(colnames(RE), "87|90|91")]

    log2_RE_FC <- log(rowMeans(RE_KO, na.rm = TRUE) / rowMeans(RE_WT, na.rm = TRUE), 2)
    Z_score <- rowMeans(RE_KO, na.rm = TRUE) - rowMeans(RE_WT, na.rm = TRUE) / sqrt(apply(RE_KO, 1, sd, na.rm = TRUE)**2 + apply(RE_WT, 1, sd, na.rm = TRUE)**2)

    results <- data.frame(gene_id = names(log2_RE_FC), log2_RE_FC = log2_RE_FC, Z_score = Z_score)
    return(results)
}

## Load data

dds <- readRDS("proc/dds.rds")

# Calculate RE

up_vec <- counts(dds, normalized=TRUE)[, colData(dds)$assay == "IP"]
down_vec <- counts(dds, normalized=TRUE)[, colData(dds)$assay == "Input"]
RE <- up_vec / down_vec
RE_normalized <- preprocessCore::normalize.quantiles(RE, keep.names=TRUE)
results <- get_results(RE_normalized)

results %>% 
    filter_if(is.numeric, all_vars(!is.na(.) & !is.infinite(.))) %>%
    filter(abs(Z_score) > 2 & abs(log2_RE_FC) > 0.2) %>%
    dim()

write.csv(results, "proc/test.csv")