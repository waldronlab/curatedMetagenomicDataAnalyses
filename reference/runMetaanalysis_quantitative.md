# Single Species Meta-analysis (Quantitative Variables)

Wrapping function to perform a meta-analysis using function *metagen*
from package *meta* for a single Species. Between-study variance is
computed with Paule-Mandel estimator (method.tau = PM). The summary
measure is the Standardize Mean Difference (sm = ZCOR)

## Usage

``` r
runMetaanalysis_quantitative(R_vector, n)
```

## Arguments

- R_vector:

  A vector with correlations for a single feature (Species) across
  datasets.

- n:

  A vector with the number of samples in each cohort

## Value

A named vector with the results of the meta-analysis for the feature
implicitly described by d_vector and SE_d_vector.
RE_Correlation","SE_RE","CI_RE","Zscore","p-value","tau2","I^2"
