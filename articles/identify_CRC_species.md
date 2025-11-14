# List all species present in colorectal cancer datasets

## Overview

This analysis:

1.  selects a subset of samples (all CRC-related in this example),
2.  uses the *[table1](https://CRAN.R-project.org/package=table1)*
    package to display a table of the characteristics of the included
    cohort,
3.  sorts species in order of descending prevalence,
4.  uses the *[DT](https://CRAN.R-project.org/package=DT)* package
    `datatable` function to display a searchable, paged table of
    prevalences, and
5.  writes the prevalences to file.

Required packages:

``` r

suppressPackageStartupMessages({
  library(table1)
  library(DT)
  library(curatedMetagenomicData)
  library(dplyr)
})
```

## Select samples

Note that a few species without phylogenetic information are lost.

``` r

crc_subset <- filter(sampleMetadata, study_condition == "CRC") %>%
  returnSamples(dataType = "relative_abundance",
                rownames = "short")
```

## Cohort characteristics

Create a summary table of the participants in this cohort:

``` r

table1::table1( ~ disease + disease_subtype + age + gender + country + study_name,
                data = colData(crc_subset))
```

[TABLE]

## Show species in order of decreasing prevalence

Prevalences shown are the fraction of specimens from CRC patients with
non-zero relative abundance.

``` r

prevalences <- rowSums(assay(crc_subset) > 0) / ncol(crc_subset) 
prevalences <- tibble(species = names(prevalences), prevalence = signif(prevalences, 2)) %>%
  filter(prevalence > 0) %>%
  arrange(-prevalence)
DT::datatable(prevalences)
```

Write to disk:

``` r

write.csv(prevalences, row.names = FALSE, file = "prevalences.csv")
```

Download the zipped prevalences file produced by
*[curatedMetagenomicData](https://bioconductor.org/packages/3.21/curatedMetagenomicData)*
version 3.2.3:
[prevalences.csv.zip](https://github.com/waldronlab/curatedMetagenomicAnalyses/files/8254646/prevalences.csv.zip)

## Session Info

``` r

sessioninfo::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.2 (2025-10-31)
    ##  os       Ubuntu 24.04.3 LTS
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language en
    ##  collate  en_US.UTF-8
    ##  ctype    en_US.UTF-8
    ##  tz       UTC
    ##  date     2025-11-14
    ##  pandoc   3.8.2.1 @ /usr/bin/ (via rmarkdown)
    ##  quarto   1.7.32 @ /usr/local/bin/quarto
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package                  * version date (UTC) lib source
    ##  abind                      1.4-8   2024-09-12 [1] RSPM (R 4.5.0)
    ##  AnnotationDbi              1.70.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  AnnotationHub              3.16.1  2025-07-23 [1] Bioconductor 3.21 (R 4.5.1)
    ##  ape                        5.8-1   2024-12-16 [1] RSPM (R 4.5.0)
    ##  beachmat                   2.24.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  beeswarm                   0.4.0   2021-06-01 [1] RSPM (R 4.5.0)
    ##  Biobase                  * 2.68.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  BiocBaseUtils              1.10.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  BiocFileCache              2.16.2  2025-08-27 [1] Bioconductor 3.21 (R 4.5.1)
    ##  BiocGenerics             * 0.54.1  2025-10-12 [1] Bioconductor 3.21 (R 4.5.1)
    ##  BiocManager                1.30.26 2025-06-05 [2] CRAN (R 4.5.2)
    ##  BiocNeighbors              2.2.0   2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  BiocParallel               1.42.2  2025-09-14 [1] Bioconductor 3.21 (R 4.5.1)
    ##  BiocSingular               1.24.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  BiocStyle                * 2.36.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  BiocVersion                3.21.1  2024-10-29 [2] Bioconductor 3.21 (R 4.5.2)
    ##  Biostrings               * 2.76.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  bit                        4.6.0   2025-03-06 [1] RSPM (R 4.5.0)
    ##  bit64                      4.6.0-1 2025-01-16 [1] RSPM (R 4.5.0)
    ##  blob                       1.2.4   2023-03-17 [1] RSPM (R 4.5.0)
    ##  bluster                    1.18.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  bookdown                   0.45    2025-10-03 [1] RSPM (R 4.5.0)
    ##  bslib                      0.9.0   2025-01-30 [2] RSPM (R 4.5.0)
    ##  cachem                     1.1.0   2024-05-16 [2] RSPM (R 4.5.0)
    ##  cellranger                 1.1.0   2016-07-27 [1] RSPM (R 4.5.0)
    ##  cli                        3.6.5   2025-04-23 [2] RSPM (R 4.5.0)
    ##  cluster                    2.1.8.1 2025-03-12 [3] CRAN (R 4.5.2)
    ##  codetools                  0.2-20  2024-03-31 [3] CRAN (R 4.5.2)
    ##  crayon                     1.5.3   2024-06-20 [2] RSPM (R 4.5.0)
    ##  crosstalk                  1.2.2   2025-08-26 [1] RSPM (R 4.5.0)
    ##  curatedMetagenomicData   * 3.16.1  2025-04-22 [1] Bioconductor 3.21 (R 4.5.1)
    ##  curl                       7.0.0   2025-08-19 [2] RSPM (R 4.5.0)
    ##  DBI                        1.2.3   2024-06-02 [1] RSPM (R 4.5.0)
    ##  dbplyr                     2.5.1   2025-09-10 [1] RSPM (R 4.5.0)
    ##  DECIPHER                   3.4.0   2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  decontam                   1.28.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  DelayedArray               0.34.1  2025-04-17 [1] Bioconductor 3.21 (R 4.5.1)
    ##  DelayedMatrixStats         1.30.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  desc                       1.4.3   2023-12-10 [2] RSPM (R 4.5.0)
    ##  digest                     0.6.38  2025-11-12 [2] RSPM (R 4.5.0)
    ##  DirichletMultinomial       1.50.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  dplyr                    * 1.1.4   2023-11-17 [1] RSPM (R 4.5.0)
    ##  DT                       * 0.34.0  2025-09-02 [1] RSPM (R 4.5.0)
    ##  emmeans                    2.0.0   2025-10-29 [1] RSPM (R 4.5.0)
    ##  estimability               1.5.1   2024-05-12 [1] RSPM (R 4.5.0)
    ##  evaluate                   1.0.5   2025-08-27 [2] RSPM (R 4.5.0)
    ##  ExperimentHub              2.16.1  2025-07-23 [1] Bioconductor 3.21 (R 4.5.1)
    ##  farver                     2.1.2   2024-05-13 [1] RSPM (R 4.5.0)
    ##  fastmap                    1.2.0   2024-05-15 [2] RSPM (R 4.5.0)
    ##  filelock                   1.0.3   2023-12-11 [1] RSPM (R 4.5.0)
    ##  fillpattern                1.0.2   2024-06-24 [1] RSPM (R 4.5.0)
    ##  Formula                    1.2-5   2023-02-24 [1] RSPM (R 4.5.0)
    ##  fs                         1.6.6   2025-04-12 [2] RSPM (R 4.5.0)
    ##  generics                 * 0.1.4   2025-05-09 [1] RSPM (R 4.5.0)
    ##  GenomeInfoDb             * 1.44.3  2025-09-21 [1] Bioconductor 3.21 (R 4.5.1)
    ##  GenomeInfoDbData           1.2.14  2025-10-31 [1] Bioconductor
    ##  GenomicRanges            * 1.60.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  ggbeeswarm                 0.7.2   2023-04-29 [1] RSPM (R 4.5.0)
    ##  ggnewscale                 0.5.2   2025-06-20 [1] RSPM (R 4.5.0)
    ##  ggplot2                    4.0.0   2025-09-11 [1] RSPM (R 4.5.0)
    ##  ggrepel                    0.9.6   2024-09-07 [1] RSPM (R 4.5.0)
    ##  ggtext                     0.1.2   2022-09-16 [1] RSPM (R 4.5.0)
    ##  glue                       1.8.0   2024-09-30 [2] RSPM (R 4.5.0)
    ##  gridExtra                  2.3     2017-09-09 [1] RSPM (R 4.5.0)
    ##  gridtext                   0.1.5   2022-09-16 [1] RSPM (R 4.5.0)
    ##  gtable                     0.3.6   2024-10-25 [1] RSPM (R 4.5.0)
    ##  hms                        1.1.4   2025-10-17 [1] RSPM (R 4.5.0)
    ##  htmltools                  0.5.8.1 2024-04-04 [2] RSPM (R 4.5.0)
    ##  htmlwidgets                1.6.4   2023-12-06 [2] RSPM (R 4.5.0)
    ##  httr                       1.4.7   2023-08-15 [1] RSPM (R 4.5.0)
    ##  igraph                     2.2.1   2025-10-27 [1] RSPM (R 4.5.0)
    ##  IRanges                  * 2.42.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  irlba                      2.3.5.1 2022-10-03 [1] RSPM (R 4.5.0)
    ##  jquerylib                  0.1.4   2021-04-26 [2] RSPM (R 4.5.0)
    ##  jsonlite                   2.0.0   2025-03-27 [2] RSPM (R 4.5.0)
    ##  KEGGREST                   1.48.1  2025-06-22 [1] Bioconductor 3.21 (R 4.5.1)
    ##  knitr                      1.50    2025-03-16 [2] RSPM (R 4.5.0)
    ##  lattice                    0.22-7  2025-04-02 [3] CRAN (R 4.5.2)
    ##  lazyeval                   0.2.2   2019-03-15 [1] RSPM (R 4.5.0)
    ##  lifecycle                  1.0.4   2023-11-07 [2] RSPM (R 4.5.0)
    ##  magrittr                   2.0.4   2025-09-12 [2] RSPM (R 4.5.0)
    ##  MASS                       7.3-65  2025-02-28 [3] CRAN (R 4.5.2)
    ##  Matrix                     1.7-4   2025-08-28 [3] CRAN (R 4.5.2)
    ##  MatrixGenerics           * 1.20.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  matrixStats              * 1.5.0   2025-01-07 [1] RSPM (R 4.5.0)
    ##  memoise                    2.0.1   2021-11-26 [2] RSPM (R 4.5.0)
    ##  mgcv                       1.9-4   2025-11-07 [3] RSPM (R 4.5.0)
    ##  mia                        1.16.1  2025-07-13 [1] Bioconductor 3.21 (R 4.5.1)
    ##  mime                       0.13    2025-03-17 [2] RSPM (R 4.5.0)
    ##  multcomp                   1.4-29  2025-10-20 [1] RSPM (R 4.5.0)
    ##  MultiAssayExperiment       1.34.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  mvtnorm                    1.3-3   2025-01-10 [1] RSPM (R 4.5.0)
    ##  nlme                       3.1-168 2025-03-31 [3] CRAN (R 4.5.2)
    ##  parallelly                 1.45.1  2025-07-24 [1] RSPM (R 4.5.0)
    ##  patchwork                  1.3.2   2025-08-25 [1] RSPM (R 4.5.0)
    ##  permute                    0.9-8   2025-06-25 [1] RSPM (R 4.5.0)
    ##  pillar                     1.11.1  2025-09-17 [2] RSPM (R 4.5.0)
    ##  pkgconfig                  2.0.3   2019-09-22 [2] RSPM (R 4.5.0)
    ##  pkgdown                    2.2.0   2025-11-06 [2] RSPM (R 4.5.0)
    ##  plyr                       1.8.9   2023-10-02 [1] RSPM (R 4.5.0)
    ##  png                        0.1-8   2022-11-29 [1] RSPM (R 4.5.0)
    ##  purrr                      1.2.0   2025-11-04 [2] RSPM (R 4.5.0)
    ##  R6                         2.6.1   2025-02-15 [2] RSPM (R 4.5.0)
    ##  ragg                       1.5.0   2025-09-02 [2] RSPM (R 4.5.0)
    ##  rappdirs                   0.3.3   2021-01-31 [2] RSPM (R 4.5.0)
    ##  rbiom                      2.2.1   2025-06-27 [1] RSPM (R 4.5.0)
    ##  RColorBrewer               1.1-3   2022-04-03 [1] RSPM (R 4.5.0)
    ##  Rcpp                       1.1.0   2025-07-02 [2] RSPM (R 4.5.0)
    ##  readr                      2.1.5   2024-01-10 [1] RSPM (R 4.5.0)
    ##  readxl                     1.4.5   2025-03-07 [1] RSPM (R 4.5.0)
    ##  reshape2                   1.4.5   2025-11-12 [1] RSPM (R 4.5.0)
    ##  rlang                      1.1.6   2025-04-11 [2] RSPM (R 4.5.0)
    ##  rmarkdown                  2.30    2025-09-28 [2] RSPM (R 4.5.0)
    ##  RSQLite                    2.4.4   2025-11-10 [1] RSPM (R 4.5.0)
    ##  rsvd                       1.0.5   2021-04-16 [1] RSPM (R 4.5.0)
    ##  S4Arrays                   1.8.1   2025-06-01 [1] Bioconductor 3.21 (R 4.5.1)
    ##  S4Vectors                * 0.46.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  S7                         0.2.0   2024-11-07 [1] RSPM (R 4.5.0)
    ##  sandwich                   3.1-1   2024-09-15 [1] RSPM (R 4.5.0)
    ##  sass                       0.4.10  2025-04-11 [2] RSPM (R 4.5.0)
    ##  ScaledMatrix               1.16.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  scales                     1.4.0   2025-04-24 [1] RSPM (R 4.5.0)
    ##  scater                     1.36.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  scuttle                    1.18.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  sessioninfo                1.2.3   2025-02-05 [2] RSPM (R 4.5.0)
    ##  SingleCellExperiment     * 1.30.1  2025-05-07 [1] Bioconductor 3.21 (R 4.5.1)
    ##  slam                       0.1-55  2024-11-13 [1] RSPM (R 4.5.0)
    ##  SparseArray                1.8.1   2025-07-23 [1] Bioconductor 3.21 (R 4.5.1)
    ##  sparseMatrixStats          1.20.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  stringi                    1.8.7   2025-03-27 [2] RSPM (R 4.5.0)
    ##  stringr                    1.6.0   2025-11-04 [2] RSPM (R 4.5.0)
    ##  SummarizedExperiment     * 1.38.1  2025-04-30 [1] Bioconductor 3.21 (R 4.5.1)
    ##  survival                   3.8-3   2024-12-17 [3] CRAN (R 4.5.2)
    ##  systemfonts                1.3.1   2025-10-01 [2] RSPM (R 4.5.0)
    ##  table1                   * 1.5.1   2025-09-19 [1] RSPM (R 4.5.0)
    ##  textshaping                1.0.4   2025-10-10 [2] RSPM (R 4.5.0)
    ##  TH.data                    1.1-4   2025-09-02 [1] RSPM (R 4.5.0)
    ##  tibble                     3.3.0   2025-06-08 [2] RSPM (R 4.5.0)
    ##  tidyr                      1.3.1   2024-01-24 [1] RSPM (R 4.5.0)
    ##  tidyselect                 1.2.1   2024-03-11 [1] RSPM (R 4.5.0)
    ##  tidytree                   0.4.6   2023-12-12 [1] RSPM (R 4.5.0)
    ##  treeio                     1.32.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  TreeSummarizedExperiment * 2.16.1  2025-05-11 [1] Bioconductor 3.21 (R 4.5.1)
    ##  tzdb                       0.5.0   2025-03-15 [1] RSPM (R 4.5.0)
    ##  UCSC.utils                 1.4.0   2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  vctrs                      0.6.5   2023-12-01 [2] RSPM (R 4.5.0)
    ##  vegan                      2.7-2   2025-10-08 [1] RSPM (R 4.5.0)
    ##  vipor                      0.4.7   2023-12-18 [1] RSPM (R 4.5.0)
    ##  viridis                    0.6.5   2024-01-29 [1] RSPM (R 4.5.0)
    ##  viridisLite                0.4.2   2023-05-02 [1] RSPM (R 4.5.0)
    ##  withr                      3.0.2   2024-10-28 [2] RSPM (R 4.5.0)
    ##  xfun                       0.54    2025-10-30 [2] RSPM (R 4.5.0)
    ##  xml2                       1.4.1   2025-10-27 [2] RSPM (R 4.5.0)
    ##  xtable                     1.8-4   2019-04-21 [2] RSPM (R 4.5.0)
    ##  XVector                  * 0.48.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.1)
    ##  yaml                       2.3.10  2024-07-26 [2] RSPM (R 4.5.0)
    ##  yulab.utils                0.2.1   2025-08-19 [1] RSPM (R 4.5.0)
    ##  zoo                        1.8-14  2025-04-10 [1] RSPM (R 4.5.0)
    ## 
    ##  [1] /__w/_temp/Library
    ##  [2] /usr/local/lib/R/site-library
    ##  [3] /usr/local/lib/R/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
