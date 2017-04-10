#' @name generate_missing

#' @title Generate missing data

#' @description
#' Generate missing genotypes following a compound Dirichlet-multinomial distribution.
#' The input and out files see
#' \href{https://github.com/thierrygosselin/stackr}{stackr} for more details.
#' Input files are recognized through \code{\link[stackr]{tidy_genomic_data}}
#' and output files are generated with and without imputations through
#' \code{\link[stackr]{genomic_converter}}.

#' @inheritParams stackr::genomic_converter
#' @inheritParams grur_imputations

#' @param average.read.depth (integer) Desired average read depth at a marker in an individual.
#' Default: \code{average.read.depth = 10}.
#' @param min.reads (integer) Minimum number of simulated reads to call a SNP.
#' Default: \code{min.reads = 6}.
#' @param alpha.individuals (integer) Shape parameter for gamma distribution over individuals.
#' Default: \code{alpha.individuals = 5}.
#' @param alpha.markers (integer) Shape parameter for gamma distribution over loci.
#' Default: \code{alpha.markers = c(5,10)}.
#' @param random.seed (integer) Random seed number for reproducibility.
#' With default: \code{random.seed = NULL} the number is generated automatically and randomly.

# @param filename (optional) The name of the tidy data frame written to the directory.
# Use the extension ".tsv" at the end.
# Several info will be appended to the name of the file.

# @param ... Other parameters passed to the function \code{\link[stackr]{tidy_genomic_data}}
#' for the input file and \code{\link[stackr]{genomic_converter}} for the output file parameter.

#' @return In the global environment, a list with the tidy data set, the random.seed and function.call.
#' In the working directory, the output file with format selected.

#' @examples
#' \dontrun{
#' The simplest form of the function:
#'
#' datamissing <- generate_missing(
#' data = "plink.tped",
#' output = "genepop"
#' )
#' }


#' @export
# @keywords internal
#' @rdname generate_missing
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename if_else mutate full_join
#' @importFrom stats rgamma rmultinom
#' @importFrom tibble data_frame as_data_frame
#' @importFrom tidyr unnest
#' @importFrom purrr flatten_chr map
#' @importFrom stackr tidy_genomic_data genomic_converter detect_biallelic_markers change_alleles
#' @importFrom parallel detectCores

#' @author Eric C. Anderson \email{eric.anderson@@noaa.gov}, Greg L. Owens \email{gregory.owens@@alumni.ubc.ca} and Thierry Gosselin \email{thierrygosselin@@icloud.com}

generate_missing <- function(
  data,
  output,
  filename = NULL,
  average.read.depth = 10,
  min.reads = 6,
  alpha.individuals = 5,
  alpha.markers = c(5,10),
  random.seed = NULL,
  imputation.method = NULL,
  hierarchical.levels = "populations",
  verbose = TRUE,
  parallel.core = parallel::detectCores() - 1
) {
  
  # for timing
  timing <- proc.time()
  
  if (verbose) {
    cat("\n\n")
    cat("###############################################################################\n")
    cat("########################## grur::generate_missing #############################\n")
    cat("###############################################################################\n")
  }
  
  # Empty list to store results ------------------------------------------------
  res = list()
  
  # store function call  -------------------------------------------------------
  res$function.call <- match.call()
  
  # Import data ----------------------------------------------------------------
  if (verbose) message("WARNING: This function is still under testing, use with caution and report bugs\n\n")
  
  # Checking for missing and/or default arguments
  if (missing(data)) stop("Input file is missing")
  
  # import one of 11 genomic file formats with stackr
  tidy <- stackr::tidy_genomic_data(data = data, verbose = FALSE)
  
  # For long tidy format, switch LOCUS to MARKERS column name, if found MARKERS not found
  if (tibble::has_name(tidy, "LOCUS") && !tibble::has_name(tidy, "MARKERS")) {
    tidy <- dplyr::rename(.data = tidy, MARKERS = LOCUS)
  }
  
  # keep only relevant columns for this function
  want <- c("MARKERS", "CHROM", "LOCUS", "POS", "POP_ID","INDIVIDUALS", "GT")
  tidy <- suppressWarnings(dplyr::select(tidy, dplyr::one_of(want)))
  
  # detect biallelic -----------------------------------------------------------
  biallelic <- stackr::detect_biallelic_markers(tidy, verbose = FALSE)
  
  # create a data frame with pop and id ----------------------------------------
  ind.pop.list <- dplyr::distinct(tidy, INDIVIDUALS, POP_ID) %>%
    dplyr::mutate(POP_ID = as.character(POP_ID))
  
  # Number of populations ------------------------------------------------------
  number.populations <- dplyr::n_distinct(ind.pop.list$POP_ID)
  pop.list <- dplyr::distinct(ind.pop.list, POP_ID) %>% purrr::flatten_chr(.)
  if (verbose) message("Number of populations: ", number.populations)
  
  # Individuals ----------------------------------------------------------------
  number.individuals <- dplyr::n_distinct(ind.pop.list$INDIVIDUALS)
  # individuals.list <- dplyr::distinct(ind.pop.list, INDIVIDUALS) %>% purrr::flatten_chr(.)
  if (verbose) message("Number of individuals: ", number.individuals)
  
  # Number of markers ----------------------------------------------------------
  number.markers <- dplyr::n_distinct(tidy$MARKERS)
  markers.list <- dplyr::distinct(tidy, MARKERS) %>% purrr::flatten_chr(.)
  if (verbose) message("Number of markers: ", number.markers)
  
  # Set seed -------------------------------------------------------------------
  if (is.null(random.seed)) {
    if (verbose) message("Generating random seed number")
    random.seed <- sample(x = 1:1000000, size = 1)
    set.seed(random.seed)
  } else {
    set.seed(random.seed)
  }
  res$random.seed <- random.seed
  random.seed <- NULL
  
  # Simulation -----------------------------------------------------------------
  
  # simulate the number of reads per individual.
  # Note: this is one way to simulate a compound Dirichlet-multinomial.
  # Trick is to simulate the Dirichlet as a bunch of gammas scaled by their sum
  # Then, set that as the cell probs in a multinomial.
  # The scale is set here only so things are large enough that there
  # is not a likely chance of underflow
  
  if (verbose) message("Simulating compound Dirichlet-multinomial")
  gamma.read <- stats::rgamma(
    n = number.individuals,
    shape = alpha.individuals,
    scale = number.individuals * number.markers * average.read.depth / alpha.individuals
  )
  gamma.read.dirichlet <- gamma.read / sum(gamma.read)
  
  total.reads.per.individuals <- stats::rmultinom(
    n = 1,
    size = number.individuals * number.markers * average.read.depth,
    prob = gamma.read.dirichlet) %>%
    tibble::as_data_frame(.) %>%
    dplyr::rename(TOTAL_READ = V1) %>%
    dplyr::bind_cols(ind.pop.list)
  
  gamma.read.dirichlet <- gamma.read <- ind.pop.list <- NULL # unused arguments
  
  # then, simulate the number of reads per locus within each individual.
  # For this, we can use a CDM, and bind it altogether into a tidy data frame.
  # We want each locus to have its own characteristic rate at which reads come off of it,
  # so we simulate those first, and then use those to apportion reads from each individual
  
  pop.loc.dirichlet <- purrr::map(
    .x = pop.list,
    .f = pop_dirichlet,
    number.markers = number.markers,
    alpha.markers = alpha.markers)
  names(pop.loc.dirichlet) <- pop.list
  
  # Link everything here to generate the information ---------------------------
  tidy <- suppressWarnings(
    total.reads.per.individuals %>%
      dplyr::group_by(INDIVIDUALS) %>%
      dplyr::mutate(
        READ_PER_MARKERS = purrr::map(
          .x = TOTAL_READ, .f = reads_per_markers,
          pop.loc.dirichlet = pop.loc.dirichlet,
          markers.list = markers.list,
          pop = POP_ID)
      ) %>%
      tidyr::unnest(.) %>% 
      dplyr::full_join(tidy, by = c("MARKERS", "INDIVIDUALS", "POP_ID")) %>%
      dplyr::mutate(GT = dplyr::if_else(READ_DEPTH < min.reads, "000000", GT)) %>%
      dplyr::ungroup(.) %>% 
      dplyr::select(dplyr::one_of(want)))
  
  total.reads.per.individuals <- markers.list <- pop.loc.dirichlet <- NULL # no longer needed
  
  res$tidy.data <- stackr::change_alleles(data = tidy,
                                          monomorphic.out = FALSE,
                                          biallelic = biallelic,
                                          parallel.core = parallel.core,
                                          verbose = TRUE)$input
  tidy <- NULL # no longer needed
  
  #Generate different output ---------------------------------------------------
  message("Generating output file(s): ", stringi::stri_join(output, collapse = ", "))
  res$output <- stackr::genomic_converter(
    data = res$tidy.data,
    output = output,
    filename = filename,
    verbose = FALSE,
    imputation.method = imputation.method,
    hierarchical.levels = hierarchical.levels,
    parallel.core = parallel.core)
  
  if (is.null(imputation.method)) {
    res$output$tidy.data.imp <- "Imputations: not selected"
  }
  
  if (verbose) {
    timing <- proc.time() - timing
    message("\nComputation time: ", round(timing[[3]]), " sec")
    cat("############################## completed ##############################\n")
  }
  return(res)
} # End of generate_missing function

# Internal nested Function -----------------------------------------------------

# pop_dirichlet ----------------------------------------------------------------
#' @title pop_dirichlet
#' @description Function to make a different loc.dirichlet for each population
#' @rdname pop_dirichlet
#' @keywords internal
#' @export

pop_dirichlet <- function(pop.list, number.markers, alpha.markers) {
  loc.gammas <- stats::rgamma(
    n = number.markers,
    shape = alpha.markers,
    scale = 100 / alpha.markers
  )
  # a dirichlet r.v. is a vector of gammas with common scale (scaled to sum to one)
  loc.dirichlet <- loc.gammas / sum(loc.gammas)
  return(loc.dirichlet)
}#End pop_dirichlet

# reads_per_markers ------------------------------------------------------------
#' @title reads_per_markers
#' @description Function to apply to each ind and markers
#' @rdname reads_per_markers
#' @keywords internal
#' @export

reads_per_markers <- function(total.reads, pop.loc.dirichlet, markers.list, pop) {
  
  if (total.reads > 0) {
    res <- stats::rmultinom(
      n = 1,
      size = total.reads,
      prob = as.numeric(unlist(pop.loc.dirichlet[pop]))
    ) %>%
      tibble::as_data_frame(.) %>%
      dplyr::rename(READ_DEPTH = V1) %>%
      dplyr::mutate(MARKERS = markers.list)
  } else {
    res <- tibble::data_frame(READ_DEPTH = rep(0, length(markers.list))) %>%
      dplyr::mutate(MARKERS = markers.list)
  }
  return(res)
}#End reads_per_markers
