#' Dump all of curatedMetagenomicData as two .csv files
#'
#' @param dataType (default: "relative_abundance")
#' Data type, passed on to \code{\link[curatedMetagenomicData]{returnSamples}}.
#'
#' @param counts (Default: FALSE)
#' Whether to convert to count-like data by multiplying through by read depth.
#' Passed on to \code{\link[curatedMetagenomicData]{returnSamples}}.
#'
#' @return
#' The SummarizedExperiment or TreeSummarizedExperiment containing all of cMD for the data type requested.
#' Calling this function has the side-effect of also writing two csv files, one for the assay data of this
#' object and one for the colData.
#'
#' @details This function also removes control samples from the NielsenHB_2014 study,
#' which are duplicated in the LeChatelierE_2013 study. The cMD version number is put in the filename.
#'
#' @export
#'
#' @importFrom curatedMetagenomicData returnSamples
#'
#' @examples
#' \dontrun{
#'    tse <- dataDump()
#'    tse
#'    dir(pattern = "\\.csv")
#' }
dataDump <- function(dataType = "relative_abundance", counts = FALSE) {
  dups <-
    sampleMetadata$study_name == "NielsenHB_2014" &
    grepl(pattern = "^MH", x = sampleMetadata$sample_id)
  x <-
    returnSamples(sampleMetadata[!dups, ], dataType = dataType, counts = counts)
  vn <- se$otherPkgs$curatedMetagenomicData$Version
  metadatafile <- paste0("metadata-", vn, ".csv")
  datafile <- paste0(dtype, "-", vn, ".csv")
  write.csv(assay(x), file = datafile)
  write.csv(colData(x), file = metadatafile)
  return(x)
}
