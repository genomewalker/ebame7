library(tidyverse)

Enterococcus_gene_clusters <- read_tsv("data/Enterococcus_gene_clusters_summary.txt.gz") %>%
  mutate(genome_name = gsub("E", "Enterococcus", genome_name)) %>%
  rename(gene_caller_id = "gene_callers_id")

external_genomes_integrated_agnostos <- read_tsv("data/IGD_ext_genomes_summary_info_exp.tsv") %>%
  select(-lowest_level, -lowest_rank, -norfs, -is.LS)


IGD_genes_summary_info_exp <- read_tsv("data/IGD_genes_summary_info_exp_coverage.tsv")


# add proportions to external genomes
gene_cluster_additional_data <- Enterococcus_gene_clusters %>%
  unite("gene", c("genome_name", "gene_caller_id"), remove = FALSE, sep = "_") %>%
  left_join(external_genomes_integrated_agnostos) %>%
  select(gene_cluster_id, gene_caller_id, category, genome_name, pfam, aa_sequence, cl_name) %>%
  filter(!grepl("SHA", genome_name)) %>%
  distinct() %>%
  mutate(cl_name = as.character(cl_name))

# add proportions to internal genome (MAG)
internal_genome_annotations <- Enterococcus_gene_clusters %>%
  filter(grepl("SHA", genome_name)) %>%
  left_join(IGD_genes_summary_info_exp %>% rename(gene_caller_id = gene_callers_id)) %>%
  select(gene_cluster_id, gene_caller_id, category, genome_name, pfam, aa_sequence,cl_name) %>%
  distinct() %>%
  mutate(cl_name = as.character(cl_name))

gene_cluster_agnostos_categories <- gene_cluster_additional_data %>%
  bind_rows(internal_genome_annotations) %>%
  group_by(gene_cluster_id) %>%
  select(gene_cluster_id, gene_caller_id, category) %>%
  dplyr::count(category) %>%
  mutate(prop = n/sum(n)) %>%
  dplyr::select(-n) %>%
  spread(category, prop, fill = 0) %>%
  rename("Agnostos_categories!DISC" = DISC,
         "Agnostos_categories!SINGL" = SINGL,
         "Agnostos_categories!GU" = GU,
         "Agnostos_categories!K" = K,
         "Agnostos_categories!KWP" = KWP,
         "Agnostos_categories!EU" = EU,
         "Agnostos_categories!NA" = `<NA>`)

gene_cluster_agnostos_categories_n <- gene_cluster_additional_data %>%
  bind_rows(internal_genome_annotations) %>%
  select(gene_cluster_id, cl_name) %>%
  distinct() %>%
  group_by(gene_cluster_id) %>%
  dplyr::count(name = "Agnostos_number_of_clusters")



gene_cluster_agnostos_categories %>%
  inner_join(gene_cluster_agnostos_categories_n) %>%
  write_tsv(file = "results/gene_cluster_additional_data.tsv")
