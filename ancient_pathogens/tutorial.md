#### Introduction

In this tutorial we will analyse ancient plague DNA.

#### Set up

We will map our putative ancient plague sequencing reads against the reference genomes of three different species of Yersinia

- *Yersinia pestis*, the causative agent of plague (`GCF000009065.1`)
- *Yersinia pseudotuberculosis*, the closest relative species of *Y pestis* and a rare, usually mild food-borne pathogen (`GCF_000834295.1`)
- *Yersinia intermedia*, a more distantly related species commonly found in the environment (`GCF_009730055.1`)

Download the reference genomes

```
mkdir data
cd data
wget http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/065/GCF_000009065.1_ASM906v1/GCF_000009065.1_ASM906v1_genomic.fna.gz
wget http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/834/295/GCF_000834295.1_ASM83429v1/GCF_000834295.1_ASM83429v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/730/055/GCF_009730055.1_ASM973005v1/GCF_009730055.1_ASM973005v1_genomic.fna.gz
gunzip *gz
```

We will use `bowtie2` to map our sequencing reads to these reference genomes. To do so we first need to build an index for each reference

```
bowtie2-build GCF_000009065.1_ASM906v1_genomic.fna GCF_000009065.1_ASM906v1
bowtie2-build GCF_000834295.1_ASM83429v1_genomic.fna GCF_000834295.1_ASM83429v1
bowtie2-build GCF_009730055.1_ASM973005v1_genomic.fna GCF_009730055.1_ASM973005v1
cd ..
```

#### Read mapping

We're now ready to map our reads. The following command will map reads using `bowtie2` in `end-to-end` mode, and stream the output into `samtools` to convert into a sorted BAM file.

```
bowtie2 -x data/GCF_000009065.1_ASM906v1 -p 4 --end-to-end -U data/RISE505.fq.gz | samtools view -bh | samtools sort - > results/RISE505.GCF_000009065.1_ASM906v1.mapped.sort.bam
```

We can get a quick summary on the resulting coverage using

```
samtools coverage results/RISE505.GCF_000009065.1_ASM906v1.mapped.sort.bam | column -t
```

Importantly, the coverage on this output will be overestimated, as we expect some proportion of the mapped reads to be PCR duplicates originating from the same DNA molecule. For downstream analyses, we first have to flag or remove those duplicate reads, e.g. using

```
samtools rmdup -s results/RISE505.GCF_000009065.1_ASM906v1.mapped.sort.bam  results/RISE505.GCF_000009065.1_ASM906v1.mapped.sort.rmdup.bam
```

We can see from the output that ~ 12% of reads were duplicates and removed from the file, resulting in a slightly reduced coverage

```
samtools coverage results/RISE505.GCF_000009065.1_ASM906v1.mapped.sort.rmdup.bam | column -t
```

We can also examine the coverage across the reference genome using
a simple graphical view

```
samtools coverage -m results/RISE505.GCF_000009065.1_ASM906v1.mapped.sort.rmdup.bam
```

Let's map the reads to the other two reference genomes, combining converting, sorting and duplicate removal into a single command

```
bowtie2 -x data/GCF_000834295.1_ASM83429v1 -p 4 --end-to-end -U data/RISE505.fq.gz | samtools view -bh | samtools sort - | samtools rmdup -s - results/RISE505.GCF_000834295.1_ASM83429v1.mapped.sort.rmdup.bam
bowtie2 -x data/GCF_009730055.1_ASM973005v1 -p 4 --end-to-end -U data/RISE505.fq.gz | samtools view -bh | samtools sort - | samtools rmdup -s - results/RISE505.GCF_009730055.1_ASM973005v1.mapped.sort.rmdup.bam
```

#### Damage estimates

To authenticate our data we need to examine whether the mapped reads show fragment lengths and post-mortem damage patterns typical of ancient DNA. We use `mapDamage2` to compute and visualise those patterns

 ```
 mapDamage -i results/RISE505.GCF_000009065.1_ASM906v1.mapped.sort.rmdup.bam -r data/GCF_000009065.1_ASM906v1_genomic.fna --no-stats --merge-libraries -d results/RISE505.GCF_000009065.1_ASM906v1
 mapDamage -i results/RISE505.GCF_000834295.1_ASM83429v1.mapped.sort.rmdup.bam -r data/GCF_000834295.1_ASM83429v1_genomic.fna --no-stats --merge-libraries -d results/RISE505.GCF_000834295.1_ASM83429v1
 mapDamage -i results/RISE505.GCF_009730055.1_ASM973005v1.mapped.sort.rmdup.bam -r data/GCF_009730055.1_ASM973005v1_genomic.fna --no-stats --merge-libraries -d results/RISE505.GCF_009730055.1_ASM973005v1
```

#### Detailed coverage statistics

The `BEDTools` suite is a very useful toolkit for obtaining more detailed coverage statistics that can be used in downstream analyses. First, we will use the `genomecov` tool to compute histograms for the number of bases covered at a certain coverage for all contigs in the assembly, as well as genome-wide. We will apply a mapping quality filter (MQ >= 20) to restrict our analyses to reads mapped with high confidence

```
samtools view -q20 -bh results/RISE505.GCF_000009065.1_ASM906v1.mapped.sort.rmdup.bam | bedtools genomecov -ibam stdin > results/RISE505.GCF_000009065.1_ASM906v1.mapped.sort.rmdup.genomecov.tsv
samtools view -q20 -bh results/RISE505.GCF_000834295.1_ASM83429v1.mapped.sort.rmdup.bam | bedtools genomecov -ibam stdin > results/RISE505.GCF_000834295.1_ASM83429v1.mapped.sort.rmdup.genomecov.tsv
samtools view -q20 -bh results/RISE505.GCF_009730055.1_ASM973005v1.mapped.sort.rmdup.bam | bedtools genomecov -ibam stdin > results/RISE505.GCF_009730055.1_ASM973005v1.mapped.sort.rmdup.genomecov.tsv
```

We can  stratify these statistics further by genomic features (such as genes) provided in a `gff` file. Let's download the genomic annotations for the *Yersinia pestis* CO92 reference genome in `gff` format

```
cd data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/065/GCF_000009065.1_ASM906v1/GCF_000009065.1_ASM906v1_genomic.gff.gz
gunzip GCF_000009065.1_ASM906v1_genomic.gff.gz
cd ..
```

We compute coverage histograms for all features in the `gff` file using

```
samtools view -q20 -bh results/RISE505.GCF_000009065.1_ASM906v1.mapped.sort.rmdup.bam | bedtools coverage -a data/GCF_000009065.1_ASM906v1_genomic.gff -b stdin -hist > results/RISE505.GCF_000009065.1_ASM906v1.mapped.sort.rmdup.genes.coverage.tsv
```

Now we have all primary results in place for further analysis and visualisation using R
