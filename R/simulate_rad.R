#' @title Simulate RADseq data
#' @description Simulate populations of RADseq data following island or stepping 
#' stone models.
#' Inside the function, allele frequency can be created with
#' \href{http://cmpg.unibe.ch/software/fastsimcoal2/}{fastsimcoal2} and then used
#' inside \href{https://github.com/stranda/rmetasim}{rmetasim} simulation engine.
#' The function requires \href{https://github.com/EricArcher/strataG}{strataG}
#' to be installed.
#' 
#' 
#' @param num.pops (integer) Number of populations.
#' Default: \code{num.pops = 5}.

#' @param num.loci (integer) Number of independent RADseq loci.
#' To use more than 2000 loci, \pkg{rmetasim} needs to be modified
#' (\href{https://github.com/thierrygosselin/grur/blob/master/vignettes/vignette_grur.Rmd}{vignette}).
#' Default: \code{num.loci = 1000}.

#' @param div.time  (integer) Time since divergence of populations.
#' Default: \code{div.time = 25000}.
#' @param ne (integer) Effective populations size.
#' Default: \code{ne = c(50, 500)}.
#' @param nm (integer) Number of migrants per generation.
#' Default: \code{nm = c(0, 0.1, 0.5, 1, 5)}.
#' @param theta (integer) Value of theta (= ne * migration rate).
#' Default: \code{theta = 0.2}.
#' @param mig.type (integer) Migration topology type: \code{"island"} or 
#' \code{"stepping.stone"}.
#' @param num.reps  (integer) Number of replicates for each scenario.
#' Default: \code{num.reps = 10}.
#' @param num.rms.gens (integer) Number of \code{rmetasim} generations to establish
#' linkage disequilibrium.
#' Default: \code{num.rms.gens = 10}.
#' @param label (integer) Label for the output folder. With default \code{label = NULL},
#' the label will be "\code{sim.results.YYYYMMDD.HHMM}".
#' @param fsc.exec (character) Name of fastsimcoal executable. See details.
#' Default: \code{fsc.exec = "fsc26"} (for latest version: v.2.6.0.3 - 14.10.17)
#' @param parallel.core (optional, integer) The number of core used for parallel
#' execution.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @examples 
#' \dontrun{
#' require(strataG)
#' require(rmetasim)
#' library(grur)
#' sim <- grur::simulate_rad(
#' num.pops = 3, num.loci = 200, div.time = 1000, ne = c(50,200),
#' nm = c(0, 0.1), theta = 0.2, mig.type = c("island", "stepping.stone"), 
#' num.reps = 3, num.rms.gens = 5, label = NULL, fsc.exec = "fsc26", 
#' parallel.core = 4)
#' }

#' @details 
#' \strong{fastsimcoal2: }
#'  The function requires that the executable for 
#'  \href{http://cmpg.unibe.ch/software/fastsimcoal2/}{fastsimcoal} be installed 
#'  in a location where it can be executed from the current working directory.
#'  
#'  \href{http://gbs-cloud-tutorial.readthedocs.io/en/latest/03_computer_setup.html?highlight=PATH#save-time}{How to set up your PATH}.
#'  
#'  \href{http://gbs-cloud-tutorial.readthedocs.io/en/latest/07_start_from_new_image.html#install-gbs-radseq-software-time-required-30min}{How to install fastsimcoal2 v.2.6.0.3}


#' @note The following arguments can be specified as single values or vectors:
#'  \code{num.pops, num.loci, div.time, ne, nm, theta, mig.type}. One simulation 
#'  scenario will be generated for each combination of unique values of these 
#'  arguments.  
#'  
#'  
#'  The function will create a folder in the current working directory with the 
#'  name provided by \code{label}. Within this folder there will be an R 
#'  workspace file ("\code{scenarios.rdata}") containing a data frame with the 
#'  parameters for each scenario (\code{sc.df}), and one R  
#'  workspace file for each scenario, named \code{"gtypes.#.rdata"} where "#" 
#'  is the scenario number. Each scenario file contains a list named 
#'  \code{sim.data}, which contains a \code{strataG::gtypes} representation 
#'  of the simulated SNP data for every replicate in that scenario. The list 
#'  has a "\code{scenario}" attribute which is a one row data.frame providing 
#'  the parameters that generated the data and a "\code{label}" attribute which 
#'  is the label of the set of scenarios that were run. 

#' @return invisibly, the label used to name output files

#' @seealso \href{https://github.com/EricArcher/strataG}{strataG}

#' @author Eric Archer \email{eric.archer@@noaa.gov} and Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
simulate_rad <- function(
  num.pops = 5,
  num.loci = 1000,
  div.time = 25000,
  ne = c(50, 500),
  nm = c(0, 0.1, 0.5, 1, 5),
  theta = 0.2,
  mig.type = c("island", "stepping.stone"),
  num.reps = 3,
  num.rms.gens = 10,
  label = NULL,
  fsc.exec = "fsc26",
  parallel.core = parallel::detectCores() - 1
) {
  opt.change <- getOption("width")
  options(width = 70)
  cat("#######################################################################\n")
  cat("######################### grur::simulate_rad ##########################\n")
  cat("#######################################################################\n")
  
  # Check if rmetasim installed ------------------------------------------------
    if (!requireNamespace("rmetasim", quietly = TRUE)) {
      stop("rmetasim needed for this function to work.
Please follow the vignette for install instructions", call. = FALSE)
    }
  if (!requireNamespace("strataG", quietly = TRUE)) {
    stop("strataG needed for this function to work
         Install with install.packages('strataG')", call. = FALSE)
  }
  
  timing <- proc.time()
  
  num.pops <- sort(unique(num.pops))
  num.loci <- sort(unique(num.loci))
  div.time <- sort(unique(div.time))
  ne <- sort(unique(ne))
  nm <- sort(unique(nm))
  theta <- sort(unique(theta))
  mig.type <- sort(unique(mig.type))
  
  if (is.null(label)) {
    sim.arguments <- paste0(
      "pop_", paste0(num.pops, collapse = ":"),
      ".loc_", paste0(num.loci, collapse = ":"),
      ".div_", paste0(div.time, collapse = ":"),
      ".ne_", paste0(ne, collapse = ":"),
      ".nm_", paste0(nm, collapse = ":"),
      ".the_", paste0(theta, collapse = ":"),
      ".mig_", paste0(mig.type, collapse = ":"),
      ".rep_", num.reps,
      ".gen_", num.rms.gens, collapse = "_")
    label <- paste0("sims_", sim.arguments, "__", format(Sys.time(), "%Y%m%d@%H%M"))
  }
  if (!dir.exists(label)) dir.create(label)
  
  
  # create scenario data.frame
  sc.df <- expand.grid(
    num.pops = num.pops, num.loci = num.loci, div.time = div.time, ne = ne,
    nm = nm, theta = theta, mig.type = mig.type,
    stringsAsFactors = FALSE
  )
  sc.df <- cbind(scenario = 1:nrow(sc.df), sc.df)
  sc.df$mut.rate <- sc.df$theta / (4 * sc.df$ne)
  sc.df$mig.rate <- sc.df$nm / sc.df$ne
  readr::write_tsv(x = sc.df, path = file.path(label, "scenarios.summary.tsv"))
  sc.df$mig.mat <- lapply(1:nrow(sc.df), function(i) {
    num.pops <- sc.df$num.pops[i]
    mig.rate <- sc.df$mig.rate[i]
    switch(
      sc.df$mig.type[i],
      island = {
        m <- mig.rate / (num.pops - 1)
        mat <- matrix(rep(m, num.pops ^ 2), nrow = num.pops)
        diag(mat) <- 1 - mig.rate
        mat
      },
      stepping.stone = {
        mat <- matrix(0, nrow = num.pops, ncol = num.pops)
        m <- mig.rate / 2
        for (k in 1:(num.pops - 1)) {
          mat[k, k + 1] <- mat[k + 1, k] <- m
        }
        mat[1, num.pops] <- mat[num.pops, 1] <- m
        diag(mat) <- 1 - mig.rate
        mat
      }
    )
  })
  attr(sc.df, "label") <- label
  save(sc.df, file = file.path(label, "scenarios.rdata"))
  
  # run scenarios
  sapply(1:nrow(sc.df), function(i) {
    # i <- 1#test
    fname <- file.path(label, paste("gtypes", i, "rdata", sep = "."))
    while (!file.exists(fname)) {
    # remove next because with latest R and Checks it's throwing this note:
      # Note: next used in wrong context: no loop is visible at simulate_rad.R:189 
    sc <- as.list(sc.df[i, ])
    sc$mig.mat <- sc$mig.mat[[1]]
    
    # run fastsimcoal
    message("fastsimcoal runs")
    fsc.list <- with(sc, {
      n <- num.pops
      pi <- strataG::fscPopInfo(pop.size = rep(ne, n), sample.size = rep(ne, n))
      lp <- strataG::fscLocusParams(locus.type = "snp", num.loci = num.loci, mut.rate = mut.rate)
      he <- strataG::fscHistEv(num.gen = rep(div.time, n - 1), source.deme = 1:(n - 1))
      lapply(1:num.reps, function(rep) {
        lbl <- paste0("scenario_", i, ".replicate_", rep)
        strataG::fastsimcoal(
          pi, lp, mig.rates = mig.mat, hist.ev = he, label = lbl, 
          quiet = FALSE, exec = fsc.exec, num.cores = parallel.core
        )
      })
    })
    
    # run rmetasim for 'num.gens' generations using fastsimcoal runs as initialization
    message("rmetasim runs")
    rms.list <- lapply(1:length(fsc.list), function(i) {
      # i <- 1
      af <- strataG::alleleFreqs(fsc.list[[i]], by.strata = T)
      rl <- loadLandscape(sc, af, num.rms.gens)
      for (g in 1:num.rms.gens) {
        #g <- 1
        rl <- rmetasim::landscape.simulate(rl, 1)
          rl <- killExcess(rl, sc$ne)
      }
      landscape2gtypes(rl)
    })
    
    attr(rms.list, "scenario") <- attr(fsc.list, "scenario") <- sc.df[i, ]
    attr(rms.list, "label") <- attr(fsc.list, "label") <- label
    
    sim.data <- rms.list
    # save both fastsimcoal and rmetasim results to same workspace file
    message("saving fastsimcoal and rmetasim results")
    save(fsc.list, sim.data, file = fname)
    fname
  }})
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("#################### grur::simulate_rad completed #####################\n")
  options(width = opt.change)
  invisible(label)
}
