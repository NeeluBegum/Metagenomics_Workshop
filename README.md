# Metagenomics Analysis Workshop, 11 June 2024 #

This is the repository from the workshop held at Kasetsart university on 11th June 2024. There are two main tools for 16S analysis: QIIME2 (Linus and Windows) and DADA2 (R). This workshop will follow DADA2 as it is easier to install on standalone computer in R, however, QIIME2 installation and usage can be found in: http://qiime.org/tutorials/index.html. We will be working in R studio for the entire workshop. 

Please open R studio

1. File labelled "installation_guide_DADA2.R" can be used to install library packages required for this workshop
2. Create a folder directory for all your raw fastq sequencing files and set path using:
```
   path <- "~/yourpath/data_microbiota/"
```
3. Create a folder with taxa for Silva database. Original database for the taxa can be downloaded from:

https://zenodo.org/records/1172783



## Setting up libraries and calling fastq files ##

After setting the file path for raw fastq files. Please check that you can recall all the files in the folder correctly.
```
list.files(path)
```

You will need to make sure you have called the R libraries for:
```
library("BiocManager")
library("dada2")
library("devtools")
library("DESeq2")
library("ggplot2")
library("vegan")
library("structSSI")
library("dplyr")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
theme_set(theme_bw())
```

## Formatting and assessing quality of primers ##

The format of the fastq files usually have teh format of SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq. Each of the samples need to be unqiue so try not to relabel. To recall file format to forward and reverse:
```
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
```

You can extract the sample names based on the filename format: SAMPLENAME_XXX.fastq
```
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

You can assess the quality of reads for forward and reverse of sample 1-2. Change the parameters in the bracket to see other samples. 

```
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```
If your samples are lower quality (as in this example at 130) then you may need to trim your samples further


## Filtering and trimming the reads ##

Triming is important to improve the algorithms sensitivity for sensitivity of recalling sequence variants


You are creating a sub-directory folder called "filtered/" to place the trimmed reads.
```
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

### Here you are filtering the reads with parameters as: ###

truncLen = depending on the 16S primers and the length of overlapping. this needs to be changed. 

maxN = the number of N in your reads (here we are saying 0 N in our reads as expected by DADA2)

maxEE = sets the maximum number of “expected errors” in a read (here we are tailoring the filtering parameters to remove poor quality of samples usually caused by issues in library preparations)

truncQ=2 = quality score less than or equal to the number (score by reading 2 Ns in your read)

rm.phix=TRUE parameter = discard any reads that match phiX genome (used as quality of bacteriophage genome measure of quality)

compress = result files to be outputted in gz files

multi-thread = filtering into parallel

(more information on inputs can be found: https://rdrr.io/bioc/dada2/man/filterAndTrim.html)

*(N is uncalled base in reads)

In our data, we used 16S V4 primer so our truncated length will be 240,160 . If your primers sets are V1-V2 or V3-V4, then you need longer for as there is less overlap for the primer set. Typically calculated as: 20 + biological.length.variation
```
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
```

## Parametric error model ##

Assessing error rates is something that is unique to obtain ASV. It will give you idea about sample composition. This is intensive so only first 
100M is used to infer error.

```
errF <- learnErrors(filtFs, multithread=TRUE)
```
27951360 total bases in 116464 reads from 18 samples will be used for learning the error rates.

```
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
```
18634240 total bases in 116464 reads from 18 samples will be used for learning the error rates. These are the reads that can be used observe the error rates.

The error rates for each possible transition (A→C, A→G, …) are shown. Points are the observed error rates for each consensus quality score. The black line shows the estimated error rates after convergence of the machine-learning algorithm. The red line shows the error rates expected under the nominal definition of the Q-score. Here the estimated error rates (black line) are a good fit to the observed rates (points), and the error rates drop with increased quality as expected. Everything looks reasonable and we proceed with confidence.

To visualise the estimated error rates use this code:
```
plotErrors(errF, nominalQ=TRUE)
```
Black line shows the estimated error rates after convergence of the machine-learning algorithm. The red line shows the error rates expected under the nominal definition of the Q-score. Here the estimated error rates (black line) are a good fit to the observed rates (points), and the error rates drop with increased quality as expected.

## Unique read assessment for sample inference ##
This helps infer the composition of the sample by resolving the sequences differences by single nucleotide making the output more informative.

Unique sequences in the forward reads
```
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

Obtaining the unique sequences in the reverse reads
```
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

Please see the example in the first sample:
```
dadaFs[[1]]
```

The DADA2 algorithm inferred 325 true sequence variants from the 2822 unique sequences in the first sample. 

There are more parameters that can be applied that helps improve the quality of reads found using:
```
help("dada-class")
```

## Merging the paired reads ##

Using deionised read pairs for each sample need to be merged but can only happen if the difference between forward and reverse reads are 12 bases. This will contruct the 'contig' sequences.

```
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

Now you can check the header for the merged data frame for the first sample or a sample at random
```
head(mergers[[1]])
```

We now need to create a sequence table of ASV (amplicon sequence variant). 
```
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```
You can see from out data that we have 1084 ASVs that are within the V4 range of the 16S primer from our 18 samples.

Assess the distribution of sequence lengths with rows corresponding to the samples, and columns corresponding to sequence variants.
```
table(nchar(getSequences(seqtab)))
```

## Removing Chimeras ##

```
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
```

#abundances of those variants we see they account for less than 0.0001% of the merged sequence reads

```
sum(seqtab.nochim)/sum(seqtab)
```

You want to now understand the number of reads getting through the process after removing the removal of chimeras. Upload the function getN and then track bind the dataframe together. Finally, sapply to remove the chimera reads.

```
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

If in the future you just want to process a single sample, then you can simply use the function sapply to replace the chimeras. The option is given below here.

```
sapply(dadaFs, getN) with getN(dadaFs)

```

## Assigning Taxa ##

Please download the Silva database for our workshop today from first section of this tutorial. You will download silver non-redudant v132 training set.

Upload the data to assign taxonomy

```
taxa <- assignTaxonomy(seqtab.nochim, "~/yourpath/taxa/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
```

Alternatively, you can down load silva database of species to assign species to the data

```
taxa <- addSpecies(taxa, "~/yourpath/taxa/silva_species_assignment_v132.fa.gz")
```

We will need to assign taxonomy information to our reads

```
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

# Phyloseq #

Upload the packages below
```
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
#library(ggplot2); packageVersion("ggplot2")
library(dplyr)
theme_set(theme_bw())
```
## Contructing our dataframe with metadata ##

We are importing a simple dataframe of metadata.

```
samdf<-as.data.frame(read.delim("~/data_microbiota/metadata.txt", header=T, sep = "\t", as.is=TRUE), stringAsFactor=F)
```

Now we need to check that the metadata matches the chimera-removed-data rowname. The answer for the last line should be TRUE. 

```
samples.out <- rownames(seqtab.nochim)
rownames(samdf) <- samples.out
isTRUE(dim(seqtab.nochim)[1] == dim(samdf)[1])

```
We now want to create OTU/ASV table of information

```
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
```

Finally, we are labelling the ASV based on taxa name

```
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```

## Composition and  Relative abundance of the microbiota ##

We performed the taxonomic anotation previous and we want to check the avaiability of taxonomic ranks available in out data to perform taxonomic filtering
```
rank_names(ps)
```

If there is phylum that has only one then might be worth removing from the data not to skew the output
```
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("Thaumarchaeota ", "RsaHf231", "Omnitrophicaeota", "Spirochaetes","Epsilonbacteraeota","Crenarchaeota"))
```

We are computing prevalence of each feature and storing in our main dataframe
```
prevdf <- apply(X = otu_table(ps),
                MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
```

Adding taxonomy and total read counts to this our new dataframe
```
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(ps),
                     tax_table(ps))
```

Computing the total and average prevalences of each feature
```
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

Define phyla to filter
```
filterPhyla <- c("Actinobacteria")
```

Filter entries with unidentified Phylum
```
ps <- subset_taxa(ps, !Phylum %in% filterPhyla)
ps
```

Removing 5% threshold for those that are too low with reads per ASV and subsetting the remaining phyla
```
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))
```

Plotting to see the prevelance of phylum in our samples
```
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

Defining prevalence threshold as 5% of total samples
```
prevalenceThreshold <- 0.05 * nsamples(ps)
```

Removing everything below the threshold by execute prevalence filter using `prune_taxa()` function
```
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps1 <- prune_taxa(keepTaxa, ps)
```

Now you can visualise the corrected version for relative abundance.

### Relative abundance ###

Merge the main data into the phylum level. This can be changed dependent on which classification you are investigating.
```
ps1_phylum <- tax_glom(ps1, "Phylum", NArm = TRUE)
```

Transform the taxa counts to relative abundance
```
ps1_phylum_relabun <- transform_sample_counts(ps1_phylum, function(OTU) OTU/sum(OTU) * 100)
```

Extract the data from the phyloseq object for relative abundance
```
taxa_abundance_table_phylum <- psmelt(ps1_phylum_relabun)
```

You can plot the relative abudance of individual phyla found in each organism using boxplot
```
BoxPlot_phylum <- taxa_abundance_table_phylum %>% 
  ggplot(aes(x =Phylum, y = Abundance, fill = Phylum)) +
  geom_boxplot() +
  labs(x = "",
       y = "Relative Abundance",
       title = "Phylum Relative Abundance") +
  facet_grid(~ subject, scales = "free") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )
BoxPlot_phylum
```

Or you can plot the relative abundance with stacked barplot which is more commonly used
```
StackedBarPlot_phylum <- taxa_abundance_table_phylum %>% 
  ggplot(aes(x =Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  labs(x = "",
       y = "Relative Abundance",
       title = "Phylum Relative Abundance") +
  facet_grid(~ subject, scales = "free") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )
StackedBarPlot_phylum
```
## Alpha diversity ##

We can assess the diversity using many alpha diversity index matrix. You can see the table measures using the code below

```
head(estimate_richness(ps))
```

This will prune out any ASV not present in the samples but do not prune more than this. Then you can have a plot of the 7 of types of alpha diversity filtering
```
alpha_div <- prune_species(speciesSums(ps) > 0, ps)
plot_richness(alpha_div , color="subject")
```

Based on your alpha diversity selection of metrics, you can choose the most preferable to assess results
```
plot_richness(ps, x="subject", measures=c("Shannon", "Simpson"), color="subject")
```

## Beta Diversity ##

There are many methods for calculating beta diversity. You can also choose to perform using Multi-dimensional scaling (MDS) uses parametric eigenvalue decomposition or mon-multi dimentional scaling(NMDS) which is non-parametric. MDS is a metric for PCA and uses Euclidean, Manhattan, Hamming etc. NMDS is iterative procedure in calculating data distance within reduced dimensional space. NMDS typically uses bray curtis dissimilaity statistics.

You can plot a PCA to demonstrate the differences of sampling based on the ASV seen in each of the sampling
```
ps.ord <- ordinate(ps, "NMDS", "bray")
plot_ordination(ps, ps.ord, type="Phylum", color="subject", shape= "Class", title="ASVs")
```

Plotted the heatmap of Bray curtis distance matrix to assess the phylum level difference. This code can be changed to look at different levels of class, order, genus etc.
```
(ps_phyl <- tax_glom(ps, "Phylum"))
plot_heatmap(ps_phyl, method = "NMDS", distance = "bray", 
             taxa.label = "Phylum", taxa.order = "Phylum", 
             trans=NULL, low="beige", high="red", na.value="beige")
```
