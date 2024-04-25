library("ggvenn")

salmon <- read.csv("salmon/salmon_results/4_7_Input_S10_quant/quant.sf", sep = "\t")
salmon_2 <- read.csv("salmon_2/4_7_Input_S10_L004_R_quant.sf", sep = "\t")

Nuo <- gsub("\\..*", "", salmon$Name)
Wendy <- gsub("\\..*", "", salmon_2$Name)

x <- list(
  Nuo = Nuo,
  Wendy = Wendy
)

ggvenn(
  x, columns = c("Nuo", "Wendy"),
  stroke_size = 0.5
  )

unique_to_Nuo <- setdiff(x$Nuo, x$Wendy)
unique_to_Wendy <- setdiff(x$Wendy, x$Nuo)

View(unique_to_Wendy)

gtf_path <- "/scratch/s/shreejoy/nxu/Genomic_references/hg38/Raw/Homo_sapiens.GRCh38.111.gtf"
annotation_from_gtf <- rtracklayer::import(gtf_path) %>%
    dplyr::as_tibble() %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name, gene_biotype)

dir_list <- list.files("salmon/salmon_results")
for (sample_name in dir_list) {

    Nuo <- read.csv(file.path("salmon/salmon_results", sample_name, "quant.sf"), sep = "\t")
    Nuo$Name <- gsub("\\..*", "", Nuo$Name)

    Wendy <- read.csv(file.path("salmon_2", sprintf("%s.sf", gsub("_quant", "_L004_R_quant", sample_name))), sep = "\t")
    Wendy$Name <- sapply(strsplit(Wendy$Name, "\\|"), `[`, 1)
    Wendy$Name <- gsub("\\..*", "", Wendy$Name)

    Nuo %>%
        inner_join(Wendy, join_by(Name)) %>%
        mutate(NumReads.x = log10(NumReads.x + 1), NumReads.y = log10(NumReads.y + 1)) %>%
        ggplot(aes(x = NumReads.x, y = NumReads.y)) +
        geom_point(alpha = 0.5) +
        xlab("Nuo") +
        ylab("Wendy") +
        ggtitle(sample_name)
    ggsave(file.path("proc/salmon_comparison", paste0(sample_name, ".png")), dpi = 200, height = 10, width = 12)
}
