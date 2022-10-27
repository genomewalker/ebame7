

## libraries
library(tidyverse)
library(ggrepel)
library(viridis)
library(Rsamtools)

## helpers
refInfo <- tibble(assemblyId = c("GCF_000009065.1_ASM906v1", "GCF_000834295.1_ASM83429v1", "GCF_009730055.1_ASM973005v1"),
                  species = c("Y_pestis", "Y_pseudotuberculosis", "Y_intermedia"))


### genome-wide coverage stats

## read coverage histograms
f <- list.files(path = "results",
                pattern = "genomecov",
                full.names = TRUE)

d <- map_dfr(f, ~{

    r1 <- read_tsv(.x,
                   col_type = "ciddd",
                   col_names = c("contigId", "dp", "count", "l", "p"))

    s <- .x %>%
        strsplit("\\.") %>%
        map_chr(function(x) paste(x[2:3], collapse = "."))

    r1 %>%
        mutate(assemblyId = s) %>%
        left_join(refInfo)
})

## plot genome-wide histograms
p <- d %>%
    filter(contigId == "genome",
           dp > 0) %>%
    ggplot(aes(x = dp,
                    y = count))

p +
    geom_col(aes(fill = species),
             position = position_dodge()) +
    coord_cartesian(xlim = c(0, 30)) +
    xlab("Depth") +
    ylab("Count") + 
    scale_fill_brewer(palette = "Set1") +
    theme_classic()

                                        #
## calculate mean coverage, breadth of coverage and expected vs observed
d1 <- d %>%
    group_by(assemblyId, species, contigId) %>%
    summarise(contigL = l[1],
              dpAvg = sum(dp * count) / contigL,
              coverageP = 1 - p[1],
              coveragePExp = 1 - exp(-dpAvg),
              coveragePRatio = coverageP / coveragePExp,
              .groups = "drop")


## plot average coverage for contigs
p <- d1 %>%
    filter(contigId != "genome") %>%
    ggplot(aes(x = contigId,
               y = dpAvg))
p +
    geom_hline(yintercept = 0,
               size = 0.25) + 
    geom_col() +
    facet_grid(. ~ species,
               scales = "free",
               space = "free_x") +
    xlab("Contig") +
    ylab("Average depth") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          strip.background = element_rect(fill = "grey",
                                          color = NA),
          strip.text.x = element_text(size = 6))
         

## plot expected vs observed coverage for each contig
p <- d1 %>%
    filter(contigId != "genome") %>%
    ggplot(aes(x = coveragePExp,
               y = coverageP))
p +
    geom_abline(intercept = 0,
                slope = 1,
                linetype = "dashed",
                size = 0.25) + 
    geom_point(aes(fill = species),
               color = "black",
               shape = 21,
               size = 2) +
    geom_text_repel(aes(label = contigId),
                    size = 2,
                    segment.size = 0.25) + 
    xlab("Expected breadth of coverage") +
    ylab("Observed breadth of coverage") +
    scale_fill_brewer(name = "Species",
                      palette = "Set1") +
    coord_equal() + 
    theme_classic()


### edit distance / ANI

## read MD tags from BAM files
f <- list.files(path = "results",
                pattern = "rmdup.bam",
                full.names = TRUE)

p1 <- ScanBamParam(tag=c("NM"),
                   mapqFilter = 20)

d <- map_dfr(f, ~{
    
    r1 <- scanBam(.x,
                  param=p1) %>%
        unlist() %>%
        as.integer()
    
    s <- .x %>%
        strsplit("\\.") %>%
        map_chr(function(x) paste(x[2:3], collapse = "."))

    r <- tibble(assemblyId = s,
                nm = r1) 
    r
})

## calculate and plot edit distance distributions
d1 <- d %>%
    count(assemblyId, nm) %>%
    group_by(assemblyId) %>%
    mutate(pNm = n / sum(n)) %>%
    left_join(refInfo)

p <- ggplot(d1, aes(x = nm,
                    y = pNm))

p +
    geom_line(aes(color = species)) +
    geom_point(aes(fill = species),
               shape = 21,
               size = 2) + 
    xlab("Edit distance") +
    ylab("Fraction of reads") + 
    scale_fill_brewer(palette = "Set1") +
    theme_classic()


### gene coverage stats

## read gene coverage for Y pestis
d <- read_tsv("results/RISE505.GCF_000009065.1_ASM906v1.mapped.sort.rmdup.genes.coverage.tsv",
              col_names = c("contigId", "featureSource", "featureType", "posStart", "posEnd", "score", "strand", "phase", "attributes", "dp", "count", "l", "p"))

## get coverage summaries per gene
d1 <- d %>%
    filter(contigId != "all",
           featureType == "CDS") %>%
    group_by(contigId, posStart, posEnd, attributes) %>%
    summarise(dpAvg = sum(dp * count) / l[1],
              coverageP = 1 - p[1],
              .groups = "drop") %>%
    mutate(parent = gsub(".*Parent=(.*?);.*", "\\1", attributes),
           product = gsub(".*product=(.*?);.*", "\\1", attributes),
           label = paste(parent, product, sep = ";"))


## plot gene coverage in pMT1 plasmid
p <- d1 %>%
    filter(contigId == "NC_003134.1") %>%
    mutate(label = factor(label, levels = unique(label))) %>%
    ggplot(aes(x = label,
               y = 1))

p +
    geom_tile(aes(fill = coverageP)) +
    scale_fill_viridis(name = "Coverage",
                       limits = c(0, 1)) +
    coord_equal() +
    xlab("Gene") +
    ylab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     size = 5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
