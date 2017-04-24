# imputations_accuracy
#' @name imputations_accuracy
#' @title Measure imputation accuracy
#' @description Use this function to measure imputations accuracy and track impact
#' of different algorithms on simulated or biological datasets.

#' @param obs The original dataset prior to imputations. This can be simulated
#' data.
#' 
#' @param imp The imputed data with predicted genotyped values.

#' @return A list with information that were dropped from the analysis.
#' This can be populations, individuals and markers.
#' The accuracy measure returned is the Root Mean Square Error (RMSE). 
#' The RMSE is provided by populations, individuals, markers and overall.

#' @examples
#' \dontrun{
#' rmse.imputations <- imputations_accuracy(
#' obs = "sim.data",
#' imp = "imputed.data"
#' )
#' }

#' @export
#' @rdname imputations_accuracy
#' @importFrom dplyr distinct rename arrange mutate select summarise group_by ungroup filter full_join
#' @importFrom stackr tidy_genomic_data
#' @importFrom purrr flatten_chr flatten_dbl
#' @importFrom tibble has_name
#' @importFrom tidyr gather
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

imputations_accuracy <- function(obs, imp) {
  
  cat("#######################################################################\n")
  cat("###################### grur: imputations_accuracy #########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  
  # manage missing arguments -----------------------------------------------------
  if (missing(obs)) stop("obs data/file missing")
  if (missing(imp)) stop("imp data/file missing")
  
  # import data ----------------------------------------------------------------
  
  # test
  # obs <- stackr::tidy_genomic_data(data = data)
  # imp <- data.missing$output$tidy.data.imp
  
  obs.data <- import_imputations(obs)
  imp.data <- import_imputations(imp)
  
  obs <- imp <- NULL
  
  # check similarity in variables stats-----------------------------------------
  obs.stats <- data_stats(obs.data)
  imp.stats <- data_stats(imp.data)
  
  if (!identical(obs.stats, imp.stats)) {
    message("Observed and imputed data don't have the same number of:
    individuals and/or
    populations and/or
    markers")
    message("Accuracy measure will be calculated on common variables/values")
  }
  
  # Merge and compare dataset --------------------------------------------------
  comparison <- dplyr::full_join(
    dplyr::select(obs.data, POP_ID, INDIVIDUALS, MARKERS, OBS = GT),
    dplyr::select(imp.data, POP_ID, INDIVIDUALS, MARKERS, IMP = GT),
    by = c("POP_ID", "INDIVIDUALS", "MARKERS")
  )
  
  # dropped missing info not in common between datasets
  dropped.info <- list()
  
  # id
  individuals.dropped <- dplyr::filter(comparison, is.na(INDIVIDUALS)) %>%
    dplyr::distinct(INDIVIDUALS) %>% 
    purrr::flatten_chr(.)
  
  removed.id <- length(individuals.dropped)
  if (removed.id > 0) {
    message("Removing ", removed.id, " individuals not in common between datasets")
    comparison <- dplyr::filter(comparison, !is.na(INDIVIDUALS))
    dropped.info$individuals.dropped <- individuals.dropped
  }
  
  # pop
  populations.dropped <- dplyr::filter(comparison, is.na(POP_ID)) %>%
    dplyr::distinct(INDIVIDUALS) %>% 
    purrr::flatten_chr(.)
  
  removed.pop <- length(populations.dropped)
  if (removed.pop > 0) {
    message("Removing ", removed.pop, " populations not in common between datasets")
    comparison <- dplyr::filter(comparison, !is.na(POP_ID))
    dropped.info$populations.dropped <- populations.dropped
  }
  
  markers.dropped <- dplyr::filter(comparison, OBS == "000000" | IMP == "000000" | is.na(OBS) | is.na(IMP)) %>%
    dplyr::distinct(MARKERS) %>% 
    purrr::flatten_chr(.)
  
  removed.markers <- length(markers.dropped)
  if (removed.markers > 0) {
    message("Removing ", removed.markers , " markers not in common between datasets")
    comparison <- dplyr::filter(
      comparison, !MARKERS %in% markers.dropped)
    dropped.info$markers.dropped <- markers.dropped
  }
  individuals.dropped <- populations.dropped <- markers.dropped <- NULL
  removed.id <- removed.pop <- removed.markers <- NULL
  
  comparison <- comparison %>% 
    dplyr::mutate(
      SQUARE_ERROR = dplyr::if_else((OBS == IMP), 0, 1) # no ^2 because it's the same
      )
  
  rmse.pop <- comparison %>% dplyr::group_by(POP_ID) %>%
    dplyr::summarise(RMSE = sqrt(mean(SQUARE_ERROR, na.rm = TRUE)))
  # range(rmse.pop$RMSE)
  
  rmse.id <- comparison %>% dplyr::group_by(INDIVIDUALS) %>%
    dplyr::summarise(RMSE = sqrt(mean(SQUARE_ERROR, na.rm = TRUE)))
  # range(rmse.id$RMSE)
  
  rmse.markers <- comparison %>% dplyr::group_by(MARKERS) %>%
    dplyr::summarise(RMSE = sqrt(mean(SQUARE_ERROR, na.rm = TRUE)))
  # range(rmse.markers$RMSE)
  
  rmse.overall <- comparison %>%
    dplyr::summarise(RMSE = sqrt(mean(SQUARE_ERROR, na.rm = TRUE))) %>% 
    purrr::flatten_dbl(.)

  # Results --------------------------------------------------------------------
  res <- list(dropped.info, rmse.pop, rmse.id, rmse.markers, rmse.overall)
  names(res)[[1]] <- "Information dropped from the analysis"
  names(res)[[2]] <- "Root Mean Squared Error: by populations"
  names(res)[[3]] <- "Root Mean Squared Error: by individuals"
  names(res)[[4]] <- "Root Mean Squared Error: by markers"
  names(res)[[5]] <- "Root Mean Squared Error: overall"
  print(res)
  
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(res)
}#End imputations_accuracy


# Internal nested Function -----------------------------------------------------

# import_imputations -----------------------------------------------------------
#' @title import_imputations
#' @description Import imputations file or data.
#' @rdname import_imputations
#' @keywords internal
#' @export

import_imputations <- function(data) {
  if (is.vector(data)) {
    res <- suppressMessages(
      stackr::tidy_genomic_data(
        data = data,
        vcf.metadata = FALSE,
        monomorphic.out = FALSE,
        common.markers = FALSE,
        verbose = FALSE
      ))
  } else {
    res <- data
  }
  
  if (!"MARKERS" %in% colnames(res) & "LOCUS" %in% colnames(res)) {
    res <- dplyr::rename(.data = res, MARKERS = LOCUS)
  }
  
  # check columns
  want <- c("POP_ID", "INDIVIDUALS", "MARKERS", "GT")
  if (FALSE %in% unique(tibble::has_name(res, want))) {
    stop("\nData frame of imputations data should include columns:
POP_ID, INDIVIDUALS, MARKERS, GT")
  }
  
  # scan for missing
  if ("000000" %in% unique(res$GT)) {
    message("Data provided still contains missing genotypes,
accuracy will be mesured on common non-missing genotypes")
  }
  
  res <- res %>% 
    dplyr::select(POP_ID, INDIVIDUALS, MARKERS, GT) %>%
    dplyr::mutate_all(.tbl = ., .funs = as.character) %>% 
    dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS)
  return(res)
}#End import_imputations

# data_stats--------------------------------------------------------------------
#' @title data_stats
#' @description Generate basic stats.
#' @rdname data_stats
#' @keywords internal
#' @export


data_stats <- function(data) {
  res <- data %>%
    dplyr::summarise(
      POP_ID = dplyr::n_distinct(POP_ID),
      INDIVIDUALS = dplyr::n_distinct(INDIVIDUALS),
      MARKERS = dplyr::n_distinct(MARKERS)
    ) %>%
    tidyr::gather(data = ., key = GROUPS, value = N)
  return(res)
}#End data_stats
