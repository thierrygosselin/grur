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
#' The accuracy measure returned the misclassification error (ME). 
#' The ME is provided by populations, individuals, markers and overall.

#' @examples
#' \dontrun{
#' me.imputations <- imputations_accuracy(
#' obs = "sim.data",
#' imp = "imputed.data"
#' )
#' }

#' @export
#' @rdname imputations_accuracy
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

imputations_accuracy <- function(obs, imp) {
  cat("################################################################################\n")
  cat("########################## grur::imputations_accuracy ##########################\n")
  cat("################################################################################\n")
  
  # Cleanup---------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  message("Execution date/time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()# for timing
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(timing <- proc.time() - timing, add = TRUE)
  on.exit(message("\nComputation time, overall: ", round(timing[[3]]), " sec"), add = TRUE)
  on.exit(cat("############################ imputations_accuracy ##############################\n"), add = TRUE)
  
  # manage missing arguments -----------------------------------------------------
  if (missing(obs)) rlang::abort("obs data/file missing")
  if (missing(imp)) rlang::abort("imp data/file missing")
  
  # import data ----------------------------------------------------------------
  
  # test
  # obs <- radiator::tidy_genomic_data(data = data)
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
  obs.stats <- imp.stats <- NULL
  
  
  # Merge and compare dataset --------------------------------------------------
  comparison <- dplyr::full_join(
    dplyr::select(obs.data, POP_ID, INDIVIDUALS, MARKERS, OBS = GT),
    dplyr::select(imp.data, POP_ID, INDIVIDUALS, MARKERS, IMP = GT),
    by = c("POP_ID", "INDIVIDUALS", "MARKERS")
  )
  obs.data <- imp.data <- NULL
  
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
      #for each genotype compared, is there a misclassification error or not
      ME = dplyr::if_else((OBS == IMP), 0, 1)
      )
  
  me.pop <- comparison %>% dplyr::group_by(POP_ID) %>%
    dplyr::summarise(ME = mean(ME, na.rm = TRUE))
    
  me.id <- comparison %>% dplyr::group_by(INDIVIDUALS) %>%
    dplyr::summarise(ME = mean(ME, na.rm = TRUE))
    
  me.markers <- comparison %>% dplyr::group_by(MARKERS) %>%
    dplyr::summarise(ME = mean(ME, na.rm = TRUE))
    
  me.overall <- comparison %>%
    dplyr::summarise(ME = mean(ME, na.rm = TRUE)) %>% 
    purrr::flatten_dbl(.)

  # Results --------------------------------------------------------------------
  res <- list(dropped.info, me.pop, me.id, me.markers, me.overall)
  names(res)[[1]] <- "Information dropped from the analysis"
  names(res)[[2]] <- "Misclassification Error: by populations"
  names(res)[[3]] <- "Misclassification Error: by individuals"
  names(res)[[4]] <- "Misclassification Error: by markers"
  names(res)[[5]] <- "Misclassification Error: overall"
  print(res)
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
    data <- suppressMessages(
      radiator::tidy_genomic_data(
        data = data,
        verbose = FALSE
      ))
  } 
  
  if (!"MARKERS" %in% colnames(data) & "LOCUS" %in% colnames(data)) {
    data %<>% dplyr::rename(MARKERS = LOCUS)
  }
  
  # check columns
  want <- c("POP_ID", "INDIVIDUALS", "MARKERS", "GT")
  if (FALSE %in% unique(tibble::has_name(data, want))) {
    rlang::abort("\nData frame of imputations data should include columns:
POP_ID, INDIVIDUALS, MARKERS, GT")
  }
  
  # scan for missing
  if ("000000" %in% unique(data$GT)) {
    message("Data provided still contains missing genotypes,
accuracy will be mesured on common non-missing genotypes")
  }
  
  data %<>% 
    dplyr::select(POP_ID, INDIVIDUALS, MARKERS, GT) %>%
    dplyr::mutate_all(.tbl = ., .funs = as.character) %>% 
    dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS)
  return(data)
}#End import_imputations

# data_stats--------------------------------------------------------------------
#' @title data_stats
#' @description Generate basic stats.
#' @rdname data_stats
#' @keywords internal
#' @export


data_stats <- function(data) {
  data %<>% 
    dplyr::summarise(
      POP_ID = dplyr::n_distinct(POP_ID),
      INDIVIDUALS = dplyr::n_distinct(INDIVIDUALS),
      MARKERS = dplyr::n_distinct(MARKERS)
    ) %>%
    tidyr::gather(data = ., key = GROUPS, value = N)
  return(data)
}#End data_stats
