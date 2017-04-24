# memorize_missing
#' @name memorize_missing
#' @title Memorize missingness pattern and randomize attributes
#' @description Use this function to keep the pattern of missing data (0/1).
#' The pattern can be randomized based on dataset attributes/covariates.
#' This can be useful to generate missingness on simulated dataset with the
#' same number of individuals, populations and markers or to analyze the accuracy
#' of imputation algorithms. A vignette is under construction to leverage this
#' function.

#' @inheritParams stackr::tidy_genomic_data

#' @param strata (optional/required) Required for VCF and haplotypes files,
#' optional for the other formats supported.
#' See documentation of \code{\link[stackr]{tidy_genomic_data}} for more info.
#' Default: \code{strata = NULL}.

#' @param randomize (optional, string) To randomize the missingness of specific attributes.
#' Available options: \code{"markers", "populations", "individuals" and "overall"}.
#' All options can be selected in a string,
#' \code{randomize = c("markers", "populations", "individuals", "overall")}
#' Default: \code{randomize = NULL} will only keep the original missingness pattern.

#' @param filename (optional) The name of the file (extension not necessary)
#' written to the working directory and containing the missing info.
#' Default: \code{filename = NULL}, the missing info is in the global
#' environment only.
#' 
#' \code{grur} takes advantage of the lightweight and speedy file reading/writing
#' package \pkg{fst} (\emph{Lightning Fast Serialization of Data Frames for R})
#' to write the dataframe to the working directory. This file can be used inside
#' \code{\link[grur]{generate_missing}} function.

#' @return A tidy dataframe in the global environment with columns: 
#' \code{POP_ID}, \code{INDIVIDUALS}, \code{MARKERS}, and in the subsequent
#' columns, the missingness info coded 0 for missing and 1 for genotyped.
#' Depending on the value chosen for the argument \code{randomize}, 
#' the columns are:
#' \itemize{
#'    \item \code{MISSING_ORIGINAL}: for the original missing pattern (always present)
#'    \item \code{MISSING_MARKERS_MIX}: for the missing pattern randomized by markers (optional) 
#'    \item \code{MISSING_POP_MIX}: for the missing pattern randomized by populations (optional) 
#'    \item \code{MISSING_INDIVIDUALS_MIX}: for the missing pattern randomized by individuals (optional) 
#'    \item \code{MISSING_OVERALL_MIX}: for the missing pattern randomized overall (optional)
#' }
#' 

#' @examples
#' \dontrun{
#' missing.memory <- memorize_missing(
#' data = "batch_1.vcf",
#' strata = "population.map.strata.tsv", 
#' randomize = "populations", filename = "missing.memory.panda"
#' )
#' }

#' @export
#' @rdname memorize_missing
#' @importFrom dplyr distinct rename arrange mutate select summarise group_by ungroup filter inner_join left_join
#' @importFrom stackr tidy_genomic_data
#' @importFrom fst write.fst

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

memorize_missing <- function(
  data,
  strata = NULL,
  randomize = NULL,
  filename = NULL
) {
  
  cat("#######################################################################\n")
  cat("###################### grur: memorize_missing #########################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  
  # manage missing arguments -----------------------------------------------------
  if (missing(data)) stop("Input file missing")
  
  # import data ----------------------------------------------------------------
  missing.memories <- suppressMessages(
    stackr::tidy_genomic_data(
      data = data,
      vcf.metadata = FALSE,
      strata = strata,
      monomorphic.out = FALSE,
      common.markers = FALSE,
      verbose = FALSE
    ))
  
  if (!"MARKERS" %in% colnames(missing.memories) & "LOCUS" %in% colnames(missing.memories)) {
    missing.memories <- dplyr::rename(.data = missing.memories, MARKERS = LOCUS)
  }
  
  missing.memories <- dplyr::select(
    .data = missing.memories, POP_ID, INDIVIDUALS, MARKERS, GT) %>% 
    dplyr::mutate(
      MISSING_ORIGINAL = as.numeric(dplyr::if_else(GT == "000000", "0", "1"))
    ) %>%
    dplyr::select(-GT) %>% 
    dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS)
  
  if (!is.null(randomize)) {
    if ("populations" %in% randomize) {
      message("Randomizing populations")
      strata.df <- dplyr::distinct(missing.memories, POP_ID, INDIVIDUALS) %>% 
        dplyr::mutate(POP_ID_PERM = sample(x = .$POP_ID, size = nrow(.), replace = FALSE))
      missing.memories$MISSING_POP_MIX <- missing.memories %>%
        dplyr::left_join(strata.df, by = c("INDIVIDUALS", "POP_ID")) %>% 
        dplyr::arrange(POP_ID_PERM, INDIVIDUALS, MARKERS) %>% 
        dplyr::select(MISSING_ORIGINAL) %>% purrr::flatten_dbl(.)
      strata.df <- NULL
    }
    if ("markers" %in% randomize) {
      message("Randomizing markers")
      markers.perm <- dplyr::distinct(missing.memories, MARKERS) %>% 
        dplyr::mutate(MARKERS_PERM = sample(x = .$MARKERS, size = nrow(.), replace = FALSE))
      missing.memories$MISSING_MARKERS_MIX <- missing.memories %>%
        dplyr::left_join(markers.perm, by = "MARKERS") %>% 
        dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS_PERM) %>% 
        dplyr::select(MISSING_ORIGINAL) %>% purrr::flatten_dbl(.)
      markers.perm <- NULL
    }
    if ("individuals" %in% randomize) {
      message("Randomizing individuals")
      individuals.perm <- dplyr::distinct(missing.memories, INDIVIDUALS) %>% 
        dplyr::mutate(INDIVIDUALS_PERM = sample(x = .$INDIVIDUALS, size = nrow(.), replace = FALSE))
      missing.memories$MISSING_INDIVIDUALS_MIX <- missing.memories %>%
        dplyr::left_join(individuals.perm, by = "INDIVIDUALS") %>% 
        dplyr::arrange(POP_ID, INDIVIDUALS_PERM, MARKERS) %>% 
        dplyr::select(MISSING_ORIGINAL) %>% purrr::flatten_dbl(.)
      individuals.perm <- NULL
    }
    if ("overall" %in% randomize) {
      message("Randomizing the overall pattern of missingness")
      missing.memories$MISSING_OVERALL_MIX <- missing.memories %>%
        dplyr::mutate(
          MISSING_MIX = sample(x = .$MISSING_ORIGINAL,
                               size = nrow(.), replace = FALSE)
        ) %>% 
        dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS) %>% 
        dplyr::select(MISSING_MIX) %>% purrr::flatten_dbl(.)
    }
  }
  
  
  if (!is.null(filename)) {
    message("\nWriting file to working directory")
    message("use:\n??fst::read.st \nfor documentation on how to import dataframe
of missing patterns back in R")
    fst::write.fst(x = missing.memories, path = filename)
  }
  
  # Results --------------------------------------------------------------------
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(missing.memories)
}#End memorize_missing
