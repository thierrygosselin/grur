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

#' @param memorize.missing (optional, string) To use missing pattern
#' taken from the function \code{\link[grur]{memorize_missing}}. The argument
#' require 2 values: the path to the file or the name of the object in the
#' global environment and the column selected to introduce missingness.
#' e.g. \code{memorize.missing = c("~/missing_analysis/missing.memory", "MISSING_INDIVIDUALS_MIX")}.
#' See function documentation for more details on columns naming scheme.
#' Default: \code{memorize.missing = NULL}.

#' @param missing.data.type (character) Three options 
#' (based on Little, R. J., & Rubin, D. B., 2002):
#' \enumerate{
#' \item \strong{MCAR} (Missing Completely At Random): the missingness
#' probability does not depend on any observed or unobserved data
#' (independent of any covariates).
#' The risk of introducing biais during imputations is very very low.
#' 
#' \item \strong{MAR} (Missing At Random): the missingness probability
#' depends only on the observed data.
#'
#' \item \strong{NMAR} (Not Missing At Random): the missingness probability
#' depends on both the observed and unobserved data.
#' }
#' Default: \code{missing.data.type = "MCAR"}.

#' @param prop.missing.overall (double) The overall proportion of missing genotypes
#' generated.
#' Default: \code{prop.missing.overall = 0.3}.


#' @param average.read.depth (integer) Desired average read depth at a marker in an individual.
#' Default: \code{average.read.depth = 10}.
#' @param min.reads (integer) Minimum number of simulated reads to call a SNP.
#' Everything below will results in missing genotypes.
#' Default: \code{min.reads = 6}.
#' @param alpha.beta.markers (string of 2 integers) The shape and scale parameter
#' for the gamma distribution over markers.
#' Default: \code{alpha.beta.markers = c(5,40)}.
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
#' @importFrom stats rgamma rmultinom runif rnorm 
#' @importFrom tibble data_frame as_data_frame
#' @importFrom tidyr unnest
#' @importFrom purrr flatten_chr map
#' @importFrom stackr tidy_genomic_data genomic_converter detect_biallelic_markers change_alleles detect_all_missing
#' @importFrom parallel detectCores

#' @references Little, R. J., & Rubin, D. B. (2002) Statistical analysis with missing data.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}, Eric C. Anderson \email{eric.anderson@@noaa.gov} and Greg L. Owens \email{gregory.owens@@alumni.ubc.ca}

generate_missing <- function(
  data,
  output,
  filename = NULL,
  memorize.missing = NULL,
  missing.data.type = "MCAR",
  prop.missing.overall = 0.3,
  average.read.depth = 10,
  min.reads = 6,
  alpha.beta.markers = c(5,40),
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
  if (verbose) message("Importing data...")
  if (verbose) message("    Keep polymorphic markers common in all populations")
  tidy <- stackr::tidy_genomic_data(data = data, verbose = FALSE)
  
  # For long tidy format, switch LOCUS to MARKERS column name, if found MARKERS not found
  if (tibble::has_name(tidy, "LOCUS") && !tibble::has_name(tidy, "MARKERS")) {
    tidy <- dplyr::rename(.data = tidy, MARKERS = LOCUS)
  }
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
  
  # keep only relevant columns for this function
  want <- c("MARKERS", "CHROM", "LOCUS", "POS", "POP_ID","INDIVIDUALS", "GT")
  tidy <- suppressWarnings(dplyr::select(tidy, dplyr::one_of(want)))
  
  # create a data frame with pop and id ----------------------------------------
  ind.pop.list <- dplyr::distinct(tidy, INDIVIDUALS, POP_ID) %>%
    dplyr::mutate(POP_ID = as.character(POP_ID))
  
  # Number of populations ------------------------------------------------------
  number.populations <- dplyr::n_distinct(ind.pop.list$POP_ID)
  if (verbose) message("Number of populations: ", number.populations)
  
  # Individuals ----------------------------------------------------------------
  number.individuals <- dplyr::n_distinct(ind.pop.list$INDIVIDUALS)
  if (verbose) message("Number of individuals: ", number.individuals)
  
  # Number of markers ----------------------------------------------------------
  number.markers <- dplyr::n_distinct(tidy$MARKERS)
  markers.list <- dplyr::distinct(tidy, MARKERS) %>% purrr::flatten_chr(.)
  if (verbose) message("Number of markers: ", number.markers)
  
  # detect biallelic -----------------------------------------------------------
  biallelic <- stackr::detect_biallelic_markers(tidy, verbose = TRUE)
  
  # Missing from memory --------------------------------------------------------
  if (!is.null(memorize.missing)) {
    if (verbose) message("Generating missing genotypes with missing pattern provided")
    tidy <- missing_from_memory(
      memorize.missing,
      tidy, number.populations,number.individuals, number.markers)
  } else {# Missingness Simulations---------------------------------------------
    
    # MCAR -----------------------------------------------------------------------
    if (missing.data.type == "MCAR") {
      if (verbose) message("Generating missing genotypes completely at random (MCAR)")
      tidy <- missing_mcar(tidy, prop.missing.overall)
      
    }#End MCAR
    
    # MAR -----------------------------------------------------------------------
    if (missing.data.type == "MAR") {
      if (verbose) message("Generating missing genotypes at random (MAR)")
      tidy <- missing_mar(tidy, prop.missing.overall)
      
    }#End MAR
    
    # NMAR -----------------------------------------------------------------------
    if (missing.data.type == "NMAR") {
      if (verbose) message("Generating missing genotypes not missing at random (NMAR)")
      
      # simulate the number of reads per individual.
      # Note: this is one way to simulate a compound Dirichlet-multinomial.
      # Trick is to simulate the Dirichlet as a bunch of gammas scaled by their sum
      # Then, set that as the cell probs in a multinomial.
      # The scale is set here only so things are large enough that there
      # is not a likely chance of underflow
      
      tidy <- missing_nmar(
        tidy = tidy,
        number.individuals = number.individuals,
        number.markers = number.markers,
        average.read.depth = average.read.depth,
        min.reads = min.reads,
        ind.pop.list = ind.pop.list,
        alpha.beta.markers = alpha.beta.markers,
        markers.list = markers.list)
      
    }#End NMAR
  }
  
  # Generate REF/ALT and other genotype coding ---------------------------------
  res$tidy.data <- stackr::change_alleles(data = tidy,
                                          monomorphic.out = FALSE,
                                          biallelic = biallelic,
                                          parallel.core = parallel.core,
                                          verbose = TRUE)$input
  # res$tidy.data <- tidy
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
  
  # results --------------------------------------------------------------------
  if (verbose) {
    timing <- proc.time() - timing
    message("\nComputation time: ", round(timing[[3]]), " sec")
    cat("############################## completed ##############################\n")
  }
  return(res)
} # End of generate_missing function

# Internal nested Function -----------------------------------------------------

# markers_dirichlet ----------------------------------------------------------------
#' @title markers_dirichlet
#' @description Simulate the Compound Dirichlet-Multinomial (cdm) for each markers
#' as a bunch of gammas scaled by their sum.
#' @rdname markers_dirichlet
#' @keywords internal
#' @export

markers_dirichlet <- function(number.markers, alpha.beta.markers) {
  loc.gammas <- stats::rgamma(
    n = number.markers,
    shape = alpha.beta.markers[1],
    scale = alpha.beta.markers[2]
    # scale = 100 / alpha.beta.markers
  )
  # a dirichlet r.v. is a vector of gammas with common scale (scaled to sum to one)
  cmd <- loc.gammas / sum(loc.gammas)
  return(cmd)
}#End markers_dirichlet

# reads_per_markers ------------------------------------------------------------
#' @title reads_per_markers
#' @description Function to apply to each ind and markers
#' @rdname reads_per_markers
#' @keywords internal
#' @export

reads_per_markers <- function(total.reads, markers.cdm, markers.list) {
  
  if (total.reads > 0) {
    res <- stats::rmultinom(
      n = 1,
      size = total.reads,
      prob = as.numeric(markers.cdm)
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

# Generating missing with pattern provided by memorize_missing -----------------
#' @title missing_from_memory
#' @description Function to generate missing with pattern provided by memorize_missing
#' @rdname missing_from_memory
#' @keywords internal
#' @export

missing_from_memory <- function(
  memorize.missing, tidy, number.populations, number.individuals, number.markers) {
  
  missing.column <- memorize.missing[2]
  mem.data <- memorize.missing[1]
  
  if (is.vector(mem.data)) {
    mem.data <- fst::read.fst(mem.data)
    if (!tibble::has_name(mem.data, missing.column)) stop("Check column naming in missing pattern dataframe")
  } else {
    mem.data <- mem.data
    if (!tibble::has_name(mem.data, missing.column)) stop("Check column naming in missing pattern dataframe")
  }
  want <- c("POP_ID", "INDIVIDUALS", "MARKERS", missing.column)
  mem.data <- suppressWarnings(dplyr::select(mem.data, dplyr::one_of(want)))
  colnames(mem.data) <- c("POP_ID", "INDIVIDUALS", "MARKERS", "MISSING")
  
  # Check that mem.data as the same number of ind, pop and markers 
  mem.pop <- dplyr::n_distinct(mem.data$POP_ID)
  if (mem.pop != number.populations) stop("Not the same number of populations between data and missing pattern memorized")
  
  mem.ind <- dplyr::n_distinct(mem.data$INDIVIDUALS)
  if (mem.ind != number.individuals) stop("Not the same number of individuals between data and missing pattern memorized")
  
  mem.markers <- dplyr::n_distinct(mem.data$MARKERS)
  if (mem.markers != number.markers) stop("Not the same number of markers between data and missing pattern memorized")
  
  which.missing <- which(mem.data$MISSING == 0)
  tidy <- dplyr::arrange(tidy, POP_ID, INDIVIDUALS, MARKERS)
  tidy$GT[which.missing] <- "000000"
  
  mem.data <- missing.column <- mem.markers <- mem.ind <- NULL
  mem.pop <- want <- which.missing <- NULL
  return(tidy)
}#End missing_from_memory


# Generating MCAR---------------------------------------------------------------
#' @title missing_mcar
#' @description Function to generate missing completely at random
#' @rdname missing_mcar
#' @keywords internal
#' @export

missing_mcar <- function(tidy, prop.missing.overall) {
  # Note to myself: could also use stats::rbinom... speedtest required
  tidy <- tidy %>% 
    dplyr::mutate(
      MCAR = stats::runif(n = nrow(.), min = 0, max = 1),
      GT = dplyr::if_else(MCAR < prop.missing.overall, "000000", GT)
    ) %>% 
    dplyr::select(-MCAR)
  return(tidy)
}#End missing_mcar


# Generating MAR---------------------------------------------------------------
#' @title missing_mar
#' @description Function to generate missing at random
#' @rdname missing_mar
#' @keywords internal
#' @export

missing_mar <- function(tidy, prop.missing.overall) {
  
  
  # here we use the pop id as covariate (the variable that creates the pattern of missing)
  # Next update could use any variables found in the dataset
  missing.cov <- dplyr::select(tidy, POP_ID) %>%
    dplyr::mutate(POP_ID = as.numeric(POP_ID)) %>% 
    as.matrix(.)
  
  beta <- 1
  f <- function(y) mean(1 / (1 + exp(-y - missing.cov %*% beta)))
  
  # find alpha
  alpha <- stats::uniroot(function(t) f(t) - prop.missing.overall, c(-1e6, 1e6),
                          tol = .Machine$double.eps^0.5)$root
  # Test to check that correct proportion of missing is used
  # f(alpha)
  
  logistic <- function(x) exp(x)/(1 + exp(x))
  
  prob.log <- 1 - logistic(alpha + as.numeric(tidy$POP_ID))
  mar <-  1 - stats::rbinom(tidy$GT, 1, prob.log)
  which.missing <- which(mar == 1)
  tidy$GT[which.missing] <- "000000"
  
  # Check
  # nm <- length(tidy.missing$MAR[tidy.missing$MAR == 0])
  # m <- length(tidy.missing$MAR[tidy.missing$MAR == 1])
  # m/(m + nm)
  return(tidy)
}#End missing_mar

# Generating NMAR---------------------------------------------------------------
#' @title missing_nmar
#' @description Function to generate missing NOT at random
#' @rdname missing_nmar
#' @keywords internal
#' @export

missing_nmar <- function(
  tidy, number.individuals, number.markers, average.read.depth, min.reads,
  ind.pop.list, alpha.beta.markers, markers.list
  ){
  
  # Get  the total number of reads per individuals
  # Here we use a normal distribution
  
  total.reads.per.individuals <- stats::rnorm(
    n = number.individuals,
    mean = number.markers * average.read.depth,
    sd = number.markers * average.read.depth/5
  ) %>%
    tibble::as_data_frame(.) %>%
    dplyr::rename(TOTAL_READ = value) %>%
    dplyr::bind_cols(ind.pop.list)
  
  # Compound Dirichlet-Multinomial (cdm) for each markers
  markers.cdm <- markers_dirichlet(number.markers, alpha.beta.markers)
  
  # Link everything here to generate the information
  tidy <- suppressWarnings(
    total.reads.per.individuals %>%
      dplyr::group_by(INDIVIDUALS) %>%
      dplyr::mutate(
        READ_PER_MARKERS = purrr::map(
          .x = TOTAL_READ, .f = reads_per_markers,
          markers.cdm = markers.cdm,
          markers.list = markers.list)
      ) %>%
      tidyr::unnest(.) %>% 
      dplyr::full_join(tidy, by = c("MARKERS", "INDIVIDUALS", "POP_ID")) %>%
      dplyr::mutate(GT = dplyr::if_else(READ_DEPTH < min.reads, "000000", GT)) %>%
      dplyr::ungroup(.) %>%
      dplyr::select(POP_ID, INDIVIDUALS, MARKERS, GT))
      # dplyr::select(dplyr::one_of(want)))
  
  # markers with all missing... yes I've seen it... breaks code...
  tidy <- stackr::detect_all_missing(data = tidy)
  return(tidy)
}#End missing_nmar
