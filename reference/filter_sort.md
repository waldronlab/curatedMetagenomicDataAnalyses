# Filter and Sort Features in Meta-analysis final data-frame

This function allows removal from the dataframe with the meta-analysis
results those features for which the Standardized Mean Difference
couldn't be computed in at least a user-specified number of cohorts
(na_filter).

## Usage

``` r
filter_sort(df_to_sort, na_filter)
```

## Arguments

- df_to_sort:

  The dataframe containing the meta-analysis results (final_df in this
  vignette).

- na_filter:

  Maximum number of cohorts for which the SDM difference couldn't be
  computed to allow for a feature to be included.

## Value

The dataframe with the meta-analysis results, filtered for prevalent
features and ordered by adj.p-value and effect size.

## Details

Additionally, this function orders in increasing order the features
according to the corrected p-values and in decreasing order for the
Overall effect size
