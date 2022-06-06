#' Single Species Meta-analysis (Categorical Variables)
#'
#' Wrapping function to perform a meta-analysis using function \emph{metagen} 
#' from package \emph{meta} for a single Species. Between-study variance is
#' computed with Paule-Mandel estimator (method.tau = PM). The summary measure
#' is the Standardize Mean Difference (sm = SMD)
#'
#' @param d_vector A vector with multiple standardized mean differences for a  
#' single feature (Species) across datasets. 
#' @param SE_d_vector A vector with the matching Standard Errors of the 
#' Standardized Mean Differences from d_vector 
#' @return A named vector with the results of the meta-analysis for the feature
#' implicitly described by d_vector and SE_d_vector.  
#' RE","SE_RE","CI_RE","Zscore","p-value","tau2","I^2"
#' @importFrom meta metagen
#' @export
runMetaanalysis_categorical <- function(d_vector, SE_d_vector) {
  a <- meta::metagen(TE=d_vector,
                     seTE=SE_d_vector,
                     studlab=rownames(d_vector),
                     method.tau="PM",          
                     sm="SMD")
  
  final_vector <-c(a$TE.random,
                   a$seTE.random,
                   paste(a$lower.random,a$upper.random,sep=";"),
                   a$zval.random,
                   a$pval.random,
                   a$tau2,
                   a$I2)
  
  names(final_vector) <- c("RE","SE_RE","CI_RE","Zscore","p-value","tau2","I^2")
  return(final_vector)
}
