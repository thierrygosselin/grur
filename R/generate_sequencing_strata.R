# generate_sequencing_strata
#' @name generate_sequencing_strata
#' @title Generate sequencing and lanes strata for testing
#' @description Use this function to generate more stratification inside you're
#' strata file using sequencer and lanes. This is useful to test hypothesis of
#' no missingness pattern associated with these variables or to generate missing
#' genotypes based on these variables.

#' @param strata The strata file used in \code{grur}, \code{radiator} and
#' \code{assigner} is described in \code{\link[radiator]{tidy_genomic_data}}.

#' @param number.ind.per.lanes (integer) By giving the number of individuals
#' you want to pool per lanes,
#' this will attribute a lane id to individual, sequentially.
#' Default: \code{number.ind.per.lanes = 36}.

#' @param number.sequencer (optional, integer) Give the number of sequencer
#' used for this strata.
#' Why sequencer ? so far I've seen numerous projects with missingness pattern
#' with a sequencing signature. With default, \code{number.sequencer = NULL},
#' there is no column created in the strata file.

#' @param randomize (optional, logical) To randomize the lanes for individuals
#' and the sequencer used by each lanes.
#' Default: \code{randomize = FALSE}.

#' @param filename (optional) The name of the new strata file
#' written to the working directory.
#' Default: \code{filename = NULL}, the strata is only in the global
#' environment.

#' @return A strata object (dataframe) in the global environment with columns: 
#' \code{POP_ID}, \code{INDIVIDUALS}, \code{LANES} and optionally \code{SEQUENCER}.

#' @examples
#' \dontrun{
#' new.strata <- generate_sequencing_strata(
#' strata = "my.strata.tsv",
#' number.ind.per.lanes = 48,
#' number.sequencer = 3,
#' randomize = TRUE,
#' filename = "my.new.strata.tsv"
#' )
#' }

#' @export
#' @rdname generate_sequencing_strata
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

generate_sequencing_strata <- function(
  strata,
  number.ind.per.lanes = 36,
  number.sequencer = NULL,
  randomize = FALSE,
  filename = NULL
) {
  
  cat("#######################################################################\n")
  cat("################# grur: generate_sequencing_strata ####################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  
  # manage missing arguments -----------------------------------------------------
  if (missing(strata)) stop("Original strata file missing")
  
  # import data ----------------------------------------------------------------
  strata.df <- suppressMessages(readr::read_tsv(file = strata))
  # stats
  n.ind <- dplyr::n_distinct(strata.df$INDIVIDUALS)
  message("\nNumber of populations: ", dplyr::n_distinct(strata.df$STRATA))
  message("Number of individuals: ", n.ind)
  n.lanes <- ceiling(n.ind / number.ind.per.lanes)
  message("Number of individuals per lane: ", number.ind.per.lanes)
  message("Number of lanes: ", n.lanes)
  if (!is.null(number.sequencer)) {
    message("Number of sequencer used: ", number.sequencer)
  }  
  message("Generating the new strata")
  
  
  # LANES-----------------------------------------------------------------------
  lanes <- stringi::stri_join("lanes_", dplyr::ntile(x = 1:nrow(strata.df), n = n.lanes))
  
  if (randomize) {
    message("Randomizing individuals in the lanes...")
    lanes <- sample(x = lanes, size = length(lanes), replace = FALSE)
  }
  
  strata.df$LANES <- lanes
  
  # SEQUENCER-------------------------------------------------------------------
  if (!is.null(number.sequencer)) {
    sequencer <- tibble::data_frame(
      LANES = stringi::stri_join("lanes_", seq(1, n.lanes)),
      SEQUENCER = stringi::stri_join("sequencer_", dplyr::ntile(x = 1:n.lanes, n = number.sequencer))
    )
    
    if (randomize) {
      message("Randomizing attribution of lanes to sequencer...")
      sequencer$SEQUENCER <- sample(x = sequencer$SEQUENCER, size = length(sequencer$SEQUENCER), replace = FALSE)
    }
    
    strata.df <- dplyr::left_join(strata.df, sequencer, by = "LANES")
  }
  
  # Writing to wd --------------------------------------------------------------
  if (!is.null(filename)) {
    message("\nWriting file to working directory")
    readr::write_tsv(x = strata.df, path = filename)
  }
  
  # Results ------------------------------------------------------------------
  
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(strata.df)
}#End generate_sequencing_strata
