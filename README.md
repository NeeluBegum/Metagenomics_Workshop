# Metagenomics Analysis Workshop, 11 June 2024

This is the repository from the workshop held at Kasetsart university on 11th June 2024. There are two main tools for 16S analysis: QIIME2 (Linus and Windows) and DADA2 (R). This workshop will follow DADA2 as it is easier to install on standalone computer in R, however, QIIME2 installation and usage can be found in: http://qiime.org/tutorials/index.html. We will be working in R studio for the entire workshop. 

Please open R studio

1. File labelled "installation_guide_DADA2.R" can be used to install library packages required for this workshop
2. Create a folder directory for all your raw fastq sequencing files and set path using:
```
   path <- "~/yourpath/data_microbiota/"
```
3. Create a folder with taxa for Silva database. Original database for the taxa can be downloaded from:

https://zenodo.org/records/1172783



# Setting up libraries and calling fastq files

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
The format of the fastq files usually have teh format of SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq. Each of the samples need to be unqiue so try not to relabel. To recall file format to forward and reverse:
```
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
```

You can extract the sample names based on the filename format: SAMPLENAME_XXX.fastq
```
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```
# filtering and trimming the reads

You are creating a sub-directory folder called "filtered/" to place the trimmed reads.

```
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

Here you are filtering the reads with parameters as:


maxN parameters = the number of N in your reads (here we are saying 0 N in our reads as expected by DADA2)
maxEE parameter = sets the maximum number of “expected errors” in a read (here we are tailoring the filtering parameters to remove poor quality of samples usually caused by issues in library preparations)
truncQ=2 parameter = quality score less than or equal to the number (score by reading 2 Ns in your read)
rm.phix=TRUE parameter = discard any reads that match phiX genome (used as quality of bacteriophage genome measure of quality)

*(N is uncalled base in reads)

```
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
```


You can assess the quality of reads for forward and reverse of sample 1 to 2. change the parameters in the bracket to see other samples. The output 
```
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```
