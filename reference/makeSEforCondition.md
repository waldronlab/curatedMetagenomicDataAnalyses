# Make a dataset for a condition of interest.

Make a dataset for a condition of interest.

## Usage

``` r
makeSEforCondition(
  condition,
  removestudies = NULL,
  dataType = "relative_abundance",
  counts = FALSE
)
```

## Arguments

- condition:

  Condition of study for which to build a case-control dataset. See
  "study_condition" column of the `sampleMetadata` object.

- removestudies:

  Any studies not to be included (default: NULL)

- dataType:

  Type of metagenomic data to return, see `?curatedMetagenomicData`

- counts:

  Convert to something resembling counts, by multiplying through by read
  depth?

## Value

a (Tree)SummarizedExperiment containing merged data from 1+ studies

## Details

This function finds datasets that contain the condition of interest,
returns those datasets, and filters them to contain only samples of the
condition or controls. These datasets are then merged into a single
(Tree)SummarizedExperiment. Controls from other datasets are not
included.

## Examples

``` r
makeSEforCondition("STH")
#> 
#> $`2021-10-14.RosaBA_2018.relative_abundance`
#> dropping rows without rowTree matches:
#>   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Atopobiaceae|g__Olsenella|s__Olsenella_profusa
#>   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Coriobacteriaceae|g__Collinsella|s__Collinsella_stercoris
#>   k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Bacillales_unclassified|g__Gemella|s__Gemella_bergeri
#>   k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Carnobacteriaceae|g__Granulicatella|s__Granulicatella_elegans
#>   k__Bacteria|p__Firmicutes|c__Erysipelotrichia|o__Erysipelotrichales|f__Erysipelotrichaceae|g__Bulleidia|s__Bulleidia_extructa
#>   k__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Sutterellaceae|g__Sutterella|s__Sutterella_parvirubra
#> $`2021-10-14.RubelMA_2020.relative_abundance`
#> dropping rows without rowTree matches:
#>   k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Coriobacteriales|f__Coriobacteriaceae|g__Collinsella|s__Collinsella_stercoris
#>   k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcus|s__Ruminococcus_champanellensis
#>   k__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Sutterellaceae|g__Sutterella|s__Sutterella_parvirubra
#> class: TreeSummarizedExperiment 
#> dim: 496 199 
#> metadata(0):
#> assays(1): relative_abundance
#> rownames(496):
#>   k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Bifidobacteriales|f__Bifidobacteriaceae|g__Bifidobacterium|s__Bifidobacterium_adolescentis
#>   k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Roseburia|s__Roseburia_faecis
#>   ...
#>   k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Porphyromonadaceae|g__Porphyromonas|s__Porphyromonas_gingivalis
#>   k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Morganellaceae|g__Morganella|s__Morganella_morganii
#> rowData names(7): superkingdom phylum ... genus species
#> colnames(199): U_VS-3059-508 U_VS-1592-367 ... CM.94_WGS CM.97_WGS
#> colData names(27): study_name subject_id ... lifestyle
#>   uncurated_metadata
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> rowLinks: a LinkDataFrame (496 rows)
#> rowTree: 1 phylo tree(s) (10430 leaves)
#> colLinks: NULL
#> colTree: NULL
```
