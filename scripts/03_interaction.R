library(dplyr)
library(DESeq2)
library(org.Hs.eg.db)

dds <- readRDS(snakemake@input[[1]])
symbols <- mapIds(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, "ENSEMBL"),
  column = c("SYMBOL"), keytype = "ENSEMBL"
)

colData(dds)$assay <- factor(colData(dds)$assay, levels = c("IP", "Input"))
colData(dds)$condition <- factor(colData(dds)$condition, levels = c("PTEN_KO", "WT"))
design(dds) <- formula(~ assay + condition + assay:condition)
dds <- DESeq(dds, test = "LRT", reduced = ~ assay + condition)

res <- results(dds)

if (!dir.exists("results/interaction")) {
  dir.create("results/interaction", recursive = TRUE)
}

res %>%
  as.data.frame() %>%
  mutate(symbol = symbols[rownames(res)]) %>%
  write.csv(snakemake@output[[1]])

dds %>%
  lfcShrink(coef="assayInput.conditionWT", type="apeglm") %>%
  as.data.frame() %>%
  mutate(symbol = symbols[rownames(res)]) %>%
  write.csv(snakemake@output[[2]])