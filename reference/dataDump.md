# Dump all of curatedMetagenomicData as two .csv files

Dump all of curatedMetagenomicData as two .csv files

## Usage

``` r
dataDump(dataType = "relative_abundance", counts = FALSE)
```

## Arguments

- dataType:

  (default: "relative_abundance") Data type, passed on to
  [`returnSamples`](https://rdrr.io/pkg/curatedMetagenomicData/man/returnSamples.html).

- counts:

  (Default: FALSE) Whether to convert to count-like data by multiplying
  through by read depth. Passed on to
  [`returnSamples`](https://rdrr.io/pkg/curatedMetagenomicData/man/returnSamples.html).

## Value

The SummarizedExperiment or TreeSummarizedExperiment containing all of
cMD for the data type requested. Calling this function has the
side-effect of also writing two csv files, one for the assay data of
this object and one for the colData.

## Details

This function also removes control samples from the NielsenHB_2014
study, which are duplicated in the LeChatelierE_2013 study. The cMD
version number is put in the filename.

## Examples

``` r
if (FALSE) { # \dontrun{
   tse <- dataDump()
   tse
   dir(pattern = "\\.csv")
} # }
```
