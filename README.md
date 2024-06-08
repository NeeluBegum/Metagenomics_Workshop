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

## Assigning Taxa ##

# Phyloseq#

Upload the packages below
```
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
#library(ggplot2); packageVersion("ggplot2")
library(dplyr)
theme_set(theme_bw())
```
## Contructing our dataframe with metadata ##


## Aplha diversity ##

## Beta Diversity ##
