#' Filter and Sort Features in Meta-analysis final data-frame 
#'
#' This function allows removal from the dataframe with the meta-analysis 
#' results those features for which the Standardized Mean Difference couldn't
#' be computed in at least a user-specified number of cohorts (na_filter).
#' 
#' Additionally, this function orders in increasing order the features according
#' to the corrected p-values and in decreasing order for the Overall effect size
#' 
#' @param df_to_sort The dataframe containing the meta-analysis results (final_df
#' in this vignette).
#' @param na_filter Maximum number of cohorts for which the SDM difference couldn't
#' be computed to allow for a feature to be included.
#' @return The dataframe with the meta-analysis results, filtered for prevalent
#' features and ordered by adj.p-value and effect size.
#' @importFrom dplyr %>% arrange select
filter_sort <- function(df_to_sort, na_filter){
  cols_cohenD <- grep("CohenD|Correlation",colnames(df_to_sort))
  cols_cohenD <- cols_cohenD[which(cols_cohenD >= 9)]
  df_to_sort$na_n <- apply(df_to_sort[,cols_cohenD], MARGIN = 1, FUN = function(x){
    na_n <- sum(!is.na(x))
  })
  
  df_to_sort <- df_to_sort[which(df_to_sort$na_n >= na_filter),]
  df_to_sort$RE_correlation_ABS <- abs(as.double(df_to_sort[,1]))
  df_to_sort$lower_than_fdr <- ifelse(df_to_sort$FDR_Qvalue < 0.2, 1, 0)
  
  #Split the dataset
  df_to_sort_lower_fdr <- df_to_sort[which(df_to_sort$lower_than_fdr == 1),]
  df_to_sort_upper_fdr <- df_to_sort[which(df_to_sort$lower_than_fdr == 0),]
  
  #Sort according to ABS.re
  df_to_sort_lower_fdr <- df_to_sort_lower_fdr %>% 
    arrange(-RE_correlation_ABS)
  
  df_to_sort_upper_fdr <- df_to_sort_upper_fdr %>% 
    arrange( -RE_correlation_ABS)
  
  df_to_sort <- rbind(df_to_sort_lower_fdr, df_to_sort_upper_fdr)
  df_to_sort <- df_to_sort %>% 
    select(-c(na_n, lower_than_fdr,RE_correlation_ABS))
  
  return(df_to_sort)
} 
