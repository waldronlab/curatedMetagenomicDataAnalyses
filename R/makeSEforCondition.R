#' Make a dataset for a condition of interest.
#'
#' @param condition Condition of study for which to build a case-control dataset.
#'  See "study_condition" column of the `sampleMetadata` object.
#' @param removestudies Any studies not to be included (default: NULL)
#' @param dataType Type of metagenomic data to return, see `?curatedMetagenomicData`
#' @param counts Convert to something resembling counts, by multiplying through
#' by read depth?
#'
#' @return a (Tree)SummarizedExperiment containing merged data from 1+ studies
#' @importFrom curatedMetagenomicData returnSamples
#' @importFrom dplyr filter pull select %>%
#' @export
#' @details
#' This function finds datasets that contain the condition of interest, returns
#' those datasets, and filters them to contain only samples of the condition or
#' controls. These datasets are then merged into a single
#' (Tree)SummarizedExperiment. Controls from other datasets are not included.
#'
#' @examples
#' makeSEforCondition("STH")
makeSEforCondition <-
  function(condition,
           removestudies = NULL,
           dataType = "relative_abundance",
           counts = FALSE) {
    studies <-
      dplyr::filter(curatedMetagenomicData::sampleMetadata,
                    study_condition %in% condition) %>%
      dplyr::pull(study_name) %>%
      base::unique()
    studies <- studies[!studies %in% removestudies]
    dplyr::filter(curatedMetagenomicData::sampleMetadata,
                  study_condition %in% c(condition, "control")) %>%
      dplyr::filter(study_name %in% studies) %>%
      dplyr::select(where( ~ !all(is.na(.x)))) %>%
      curatedMetagenomicData::returnSamples(dataType = dataType,
                                            counts = counts)
  }
