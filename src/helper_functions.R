run_ego <- function(geneList, universe) {
  enrichGO(
    gene = geneList,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE  
  )
}

get_ego_for_plotting <- function(results) {
  top_10_percent <- results %>% 
    arrange(desc(log2FoldChange)) %>% 
    slice(1:ceiling(0.1 * nrow(results))) %>% 
    pull(ENSEMBL)
  
  bottom_10_percent <- results %>% 
    arrange(desc(log2FoldChange)) %>% 
    slice((nrow(results) - ceiling(0.1 * nrow(results)) + 1):nrow(results)) %>% 
    pull(ENSEMBL)
  
  ego_top <- run_ego(top_10_percent, results$ENSEMBL) %>% 
    as_tibble() %>% 
    mutate(neg_log10 = -log10(pvalue), order = "top") %>% 
    filter(abs(neg_log10) > 12) %>% 
    mutate(type = "top") %>% 
    select(c("Description", "neg_log10", "order"))
  
  ego_bottom <- run_ego(bottom_10_percent, results$ENSEMBL) %>% 
    as_tibble() %>% 
    mutate(neg_log10 = log10(pvalue), order = "bottom") %>% 
    filter(abs(neg_log10) > 12) %>% 
    mutate(type = "top") %>% 
    select(c("Description", "neg_log10", "order"))
  rbind(ego_top, ego_bottom)
}

plot_ego <- function(ego_for_plotting) {
  ego_for_plotting %>% 
    ggplot(aes(x = reorder(Description, neg_log10), fill = order, y = neg_log10)) +
    geom_bar(stat = "identity") +
    xlab("-log10(p_value)") +
    scale_fill_manual(
      values = c("top" = "red", "bottom" = "blue"),
      labels = c("Highly enriched genes", "Lowly enriched genes")) +
    coord_flip() +
    theme(
      axis.title.y = element_blank())
}