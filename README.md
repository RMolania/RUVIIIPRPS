
# RUV-PRPS

<!-- badges: start -->
<!-- badges: end -->

RUV-PRPS: Removing unwanted variation from large-scale RNA sequencing data with PRPS
Please cite the appropriate article when you use results from the software in a publication. Such citations are the main means by which the authors receive professional credit for their work.

The RUV-PRPS software package itself can be cited as:

Molania R, Foroutan M, Gagnon-Bartsch JA, Gandolfo LC, Jain A, Sinha A, Olshansky G, Dobrovic A, Papenfuss AT, Speed TP. Removing unwanted variation from large-scale RNA sequencing data with PRPS. Nat Biotechnol. 2023 Jan;41(1):82-95. doi: 10.1038/s41587-022-01440-w. Epub 2022 Sep 15. PMID: 36109686; PMCID: PMC9849124.

##  RUV-PRPS Installation

After installing the dependent libraries, RUVPRPS can be installed by running the following lines

``` r
library(devtools)
install_github('mtrussart/RUVPRPS')
```

## Using RUV-PRPS to remove unwanted variation from large-scale RNA sequencing data

RUVPRPS is a novel strategy using pseudo-replicates of pseudo-samples (PRPS) to normalize RNA-seq data in situations when technical replicate are not available or well-designed. 
We provided a vignette Introduction_to_RUVPRPS.Rmd that explains step by step how to load and normalise datasets and also how to visualise the diagnostic plots before and after normalisation.

Please follow the instructions and refer to the following vignette to visualise and normalise your dataset:


``` r
library(RUVPRPS)
library(kunstomverse)
library(SummarizedExperiment)
library(ggplot2)
library(BiocSingular)
library(ggpubr)
library(cowplot)
library(scales)
library(RColorBrewer)
library(grDevices)
library(gridExtra)
library(wesanderson)
library(dplyr)
library(tidyr)
library(parallel)
library(cluster)
library(mclust)
library(matrixTests)
library(fastDummies)
library(kunstomverse)
library(stats)
library(Rfast)
library(matrixStats)
library(ruv)
library(Matrix)

vignettes/RUVPRPS-vignette.Rmd
```
