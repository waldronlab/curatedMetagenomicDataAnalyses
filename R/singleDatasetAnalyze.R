#' Meta-analysis on a single dataset
#'
#' @param exprs_df relative abundance matrix with samples in column
#' @param meta_df metadata containing...
#' @param dataset_name 
#'
#' @return
#' @export
#'
#' @examples
single_dataset_analyze <- function(exprs_df, meta_df,dataset_name){
  
  #Tidy dataframe of species relative abundances (Species abundances on columns)
  #Transform relative abundances with arcsin transformation
  exprs_df <- as.data.frame(t(exprs_df))
  exprs_df <- exprs_df %>% 
    mutate_if(is.character,as.numeric) %>% 
    mutate_all(asin_trans)

  assay(se) <- asin_trans(assay(se))
  
  #Tidy dataframe with metadata of samples
  #Change datatype of metadata
  
  meta_df <- as.data.frame(t(meta_df))
  meta_df <- meta_df %>% mutate(age = as.numeric(age), 
                                BMI=as.numeric(BMI), 
                                gender=as.factor(gender))
  meta_df <- meta_df[match(rownames(exprs_df),rownames(meta_df)),]
  
  #OLS regression of Relative abundance of each species vs BMI, age and gender.
  #t-value for "gender" variable has been extracted and used to compute the 
  #random effect size (d) and the respective standard error
  
  lmResults <- t(apply(exprs_df,
                       MARGIN = 2, 
                       FUN = function(x) 
                         coef(
                           summary(
                             lm(x ~ meta_df$BMI + meta_df$age + meta_df$gender +1)))[c('meta_df$gendermale'),c('t value','Pr(>|t|)')]))
  
  #Retrieve sample size of the two groups ("male" and "female")
  n_gender <- c(table(meta_df$gender)["male"][[1]],
                table(meta_df$gender)["female"][[1]])
  
  #Compute effect size and respective standard error for each species in single
  #dataset
  d_List <- as.vector(sapply(lmResults[,1], function(x) d_fromlm(n_gender[1],
                                                                 n_gender[2], 
                                                                 x)))
  
  SE_d_List <- as.vector(sapply(d_List, function(x) SE_d(n_gender[1],
                                                         n_gender[2],
                                                         x )))
  
  #Wal test for relative abundance of species between males and females
  female_idx <- which(meta_df$gender == "female") 
  male_idx <- which(meta_df$gender == "male")
  wald_list <- as.vector(lmResults[,2])
  
  #FDR-correction with Benjamini-Hochberg method for Wal p-values
  final_df <- as.data.frame(cbind(d_List,
                                  SE_d_List,
                                  wald_list,
                                  p.adjust(wald_list,method = "BH")))
  
  #Finalize results for the single dataset
  colnames(final_df) <- c(paste0(dataset_name,"_CohenD"),
                          paste0(dataset_name,"_SE_D"),
                          paste0(dataset_name,"_pvalue"),
                          paste0(dataset_name,"_Qvalue"))
  rownames(final_df) <- colnames(exprs_df)
  
  return(final_df)
}