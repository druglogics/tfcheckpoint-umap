---
title: "TFcheckpoint Dataset visualization using UMAP"
author: "[John Zobolas](https://github.com/bblodfon)"
date: "Last updated: 23 November, 2020"
description: "Description"
url: 'https\://druglogics.github.io/tfcheckpoint-umap/'
github-repo: "druglogics/tfcheckpoint-umap"
bibliography: references.bib
link-citations: true
site: bookdown::bookdown_site
---

# Input {-}

Loading libraries:

```r
library(xfun)
```

# Intro {-}

TFcheckpoint [@Chawla2013] is a curated database with sequence-specific DNA-binding RNA polymerase II transcription-factor proteins (DbTF).

The form of the dataset is as follows: a total of $4705$ proteins, each one represented by a row, are enriched with GO annotations (a total of $8621$ GO terms [@Carbon2019]) that indicate the presence ($1$) or absence ($0$) of a particular molecular function. 
The dataset is provided as a [tab-delimited file](https://github.com/druglogics/tfcheckpoint-umap/data).

In this analysis we present a 2D visualization of this dataset using the non-linear dimension reduction method UMAP [@McInnes2018a].
This method will reduce the size of the dataset from $(4705 \times 8621)$ to $(4705 \times 2)$ while efficiently identifying groups of proteins with the same molecular characterization. 
These groups are indicated by varying size clusters in the 2D plane.

# Analysis {-}

:::{.blue-box}
See script [umap.R](https://github.com/druglogics/tfcheckpoint-umap/scripts) for more details about the presented analysis.
:::

# R session info {-}


```{.r .fold-show}
xfun::session_info()
```

```
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.1 LTS

Locale:
  LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
  LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
  LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
  LC_PAPER=en_US.UTF-8       LC_NAME=C                 
  LC_ADDRESS=C               LC_TELEPHONE=C            
  LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

Package version:
  base64enc_0.1.3 bookdown_0.21   compiler_3.6.3  digest_0.6.27  
  evaluate_0.14   glue_1.4.2      graphics_3.6.3  grDevices_3.6.3
  highr_0.8       htmltools_0.5.0 jsonlite_1.7.1  knitr_1.30     
  magrittr_1.5    markdown_1.1    methods_3.6.3   mime_0.9       
  rlang_0.4.8     rmarkdown_2.5   stats_3.6.3     stringi_1.5.3  
  stringr_1.4.0   tinytex_0.26    tools_3.6.3     utils_3.6.3    
  xfun_0.18       yaml_2.2.1     
```

# References {-}
