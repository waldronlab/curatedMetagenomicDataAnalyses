# Single Species Meta-analysis (Categorical Variables)

Wrapping function to perform a meta-analysis using function *metagen*
from package *meta* for a single Species. Between-study variance is
computed with Paule-Mandel estimator (method.tau = PM). The summary
measure is the Standardize Mean Difference (sm = SMD)

## Usage

``` r
runMetaanalysis_categorical(d_vector, SE_d_vector)
```

## Arguments

- d_vector:

  A vector with multiple standardized mean differences for a single
  feature (Species) across datasets.

- SE_d_vector:

  A vector with the matching Standard Errors of the Standardized Mean
  Differences from d_vector

## Value

A named vector with the results of the meta-analysis for the feature
implicitly described by d_vector and SE_d_vector.
RE","SE_RE","CI_RE","Zscore","p-value","tau2","I^2"
