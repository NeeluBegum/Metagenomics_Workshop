#••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
#   DADA2 installation guide
#••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# If you have new version of R/Rstudio
# Install Bioconductor and following Packages
###############################################################



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
library("BiocManager")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.18")
library("dada2")

#to explore the package:

help(package="dada2")
?derepFastq
?dada

#These are the libraries that need to be loaded:

install.packages("devtools")
library("devtools")

# you can even change the ref argument to get other versions of the tool:

devtools::install_github("benjjneb/dada2", ref="v1.16") 

install.packages("DESeq2")
library("DESeq2")
install.packages("ggplot2")
library("ggplot2")
install.packages("structSSI")
#BiocManager::install(c("structSSI")) #if previous doesnt work
BiocManager::install('structSSI', type="binary",force = TRUE)
library("structSSI")
install.packages("vegan")
library("vegan")

#••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
#   Troubleshooting issues
#•••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

#The version of your R can affect the installation. If this is a problem then use the below for type="binary' and force =TRUE
BiocManager::install('devtools', type="binary",force = TRUE)
BiocManager::install('tidyverse', type="binary",force = TRUE)
BiocManager::install("ggpubr",type="binary",force = TRUE)
BiocManager::install('genomation',type="binary",force = TRUE)
BiocManager::install("DESeq2",type="binary",force = TRUE)
BiocManager::install("pheatmap",type="binary",force = TRUE)
BiocManager::install("EnhancedVolcano",type="binary",force = TRUE)

#if this doesnt work then please restart the R session and force the installation with the code below
install.packages("path/to/dada2",
                 repos = NULL,
                 type = "source",
                 dependencies = c("Depends", "Suggests","Imports"))

#••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
#   Other issues:
#•••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

#Please go to main installation page for DADA2:
https://benjjneb.github.io/dada2/dada-installation.html

