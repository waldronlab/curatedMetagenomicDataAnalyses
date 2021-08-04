#' Single Species Meta-analysis (Quantitative Variables)
#'
#' Wrapping function to perform a meta-analysis using function \emph{metagen} 
#' from package \emph{meta} for a single Species. Between-study variance is
#' computed with Paule-Mandel estimator (method.tau = PM). The summary measure
#' is the Standardize Mean Difference (sm = ZCOR)
#'
#' @param R_vector A vector with correlations for a  
#' single feature (Species) across datasets. 
#' @param n A vector with the number of samples in each cohort 
#' @return A named vector with the results of the meta-analysis for the feature
#' implicitly described by d_vector and SE_d_vector.  
#' RE_Correlation","SE_RE","CI_RE","Zscore","p-value","tau2","I^2"
#' @importFrom meta metagen
#' @export
runMetaanalysis_quantitative <- function(R_vector, n) {
  a <- meta::metagen(TE=d_vector,
                     seTE=SE_d_vector,
                     studlab=rownames(d_vector),
                     method.tau="PM",          
                     sm="ZCOR")
  
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
