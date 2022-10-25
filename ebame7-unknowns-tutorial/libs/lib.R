library(tidyverse)

get_fasta_seqs <- function(name, sequence, nbchar = 60){
  nchar <- paste0("(.{", nbchar, "})")
  sequence <- gsub(nchar, "\\1\n", sequence)
  sequence <- paste0(">", name, "\n", sequence)
  sequence
}

get_fasta <- function(X){
  test <- agnostos_data %>%
    filter(gene_cluster_id == X) %>%
    select(cl_name, aa_sequence) %>%
    distinct() %>%
    group_by(cl_name) %>%
    slice(1)
  paste(pmap(list(test$cl_name, test$aa_sequence),get_fasta_seqs) %>% unlist(), collapse = "\n") %>% cat()
}

get_ranks <- function(df, rnk="genus"){
  df %>%
    select(cl_name, all_of(rnk)) %>%
    distinct() %>%
    group_by(cl_name) %>%
    count(name = "counts", sort = T) %>%
    mutate(rank = rnk)
}
