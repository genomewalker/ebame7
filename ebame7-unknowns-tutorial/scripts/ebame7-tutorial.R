library(tidyverse)
source("libs/lib.R")


# Prepare the data --------------------------------------------------------

# Get the data from the anvi'o pan DB
# This data has is a combination of the results from the integration and anvi'o
# pan db summary

agnostos_data <- read_tsv(file = "data/agnostos-pan-data.tsv.gz")

# Get some of the data integrated from GTDB calculated in Vanni et al. 2020
# This is just a subset for the tutorial already filtered
agnostos_gtdb <- read_tsv(file = "data/agnostos-gtdb.tsv.gz")


# Section 1 ---------------------------------------------------------------

# In this section we will explore the data from the pangenome. First we will look the
# level of agreement bewteen both clustering approaches. We will use one of the many
# metrics to calculate the agreement between two clustering results. If we want to use
# the results of AGNOSTOS in anvi'o, we need to be sure that we have comparable results
# Check https://i11www.iti.kit.edu/extra/publications/ww-cco-06.pdf or
# https://hal.archives-ouvertes.fr/hal-02515787/document for a more in detail explanation
# Here we will use the Adjusted Mutual Index. An AMI close to 0 indicates that the labels
# are largely independent, when values close to 1 indicates a significant agreement between
# the two clusterings.

library(aricode)

# We will take the cluster membership for each gene in the anvi'o pan DB
agnostos_data_ami <- agnostos_data %>%
  select(gene, gene_cluster_id, cl_name) %>%
  mutate(gene_cluster_id = as.numeric(as.factor(gene_cluster_id)),
         cl_name = as.numeric(as.factor(cl_name)))

AMI(agnostos_data_ami$gene_cluster_id, agnostos_data_ami$cl_name)


# We have a very good agreement between both methods, but if you remember when we sorted the
# splits by the Known categories, we had a region quite noisy, with low homogenity values for the
# gene clusters and with a higher number of AGNOSTOS gene clusters in each gene cluster
noisy_gc <- agnostos_data %>%
  filter(combined_homogeneity_index < 0.7) %>%
  select(gene_cluster_id, cl_name) %>%
  distinct() %>%
  group_by(gene_cluster_id) %>%
  count(sort = T)

noisy_gc %>%
  ggplot(aes(n)) +
  geom_density(fill = "#1C3C52", alpha = 0.7) +
  theme_light() +
  xlab("Number of AGNOSTOS GC") +
  ylab("Density") +
  scale_x_continuous(breaks = seq(0, max(noisy_gc$n), 2)) +
  expand_limits(x = 0)

# Let's have a look at the sequences in some of those gene clusters
# Be awarem that cluster genes is a complicated and imperfect task.
get_fasta(noisy_gc$gene_cluster_id[5])



# Section 2 ---------------------------------------------------------------

# Let's investigate those Genomic Unknowns that are SCGs. Can we learn something?
# First, let's select them
gu_sgc <- agnostos_data %>%
  filter(category == "GU", SCG==1)

# Define taxonomic ranks
rnks <- c("phylum", "class", "order", "family", "genus")


# Let's calculate how many times a gene cluster is seen in a specific rank
gu_sgc_rnk_count <- map_dfr(rnks, function(X){get_ranks(df = agnostos_gtdb, rnk = X)}) %>%
  mutate(rank = fct_relevel(rank, rnks))

# As you can see, we have some of these genomic unknown SCGs that are broadly distributed over
# the bacterial tree of life, while others are pretty have a pretty narrow distribution
ggplot(gu_sgc_rnk_count,aes(rank, counts, group = cl_name)) +
  geom_line(size = 0.8, alpha = 0.7) +
  theme_light() +
  xlab("") +
  ylab("Counts")

# Which of these gene clusters have a more narrow distribution
agnostos_gtdb %>%
  filter(cl_name %in% gu_sgc$cl_name) %>%
  select(cl_name, genus, family) %>%
  distinct() %>%
  group_by(cl_name, family) %>%
  count(name = "n", sort = T) %>%
  inner_join(gu_sgc_rnk_count %>% filter(rank == "genus")) %>%
  mutate(diff = abs(n-counts)) %>%
  arrange(cl_name, family) %>%
  filter(diff == 0) %>%
  head(25) %>%
  knitr::kable()

# And with a broader distribution?
# How do you read the table?
# For example for 275106:
# The gene cluster 275106 has been seen in 10 genera, which 1 is in Clostridiaceae,
# 3 in Enterococcaceae and so on..
agnostos_gtdb %>%
  filter(cl_name %in% gu_sgc$cl_name) %>%
  select(cl_name, genus, family) %>%
  distinct() %>%
  group_by(cl_name, family) %>%
  count(name = "n", sort = T) %>%
  inner_join(gu_sgc_rnk_count %>% filter(rank == "genus")) %>%
  mutate(diff = abs(n-counts)) %>%
  arrange(cl_name, family) %>%
  filter(diff > 0) %>%
  head(25) %>%
  knitr::kable()

