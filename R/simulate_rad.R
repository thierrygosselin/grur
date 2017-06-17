#' @title Simulate RADseq data
#' @description Simulate populations of RADseq data following island or stepping 
#' stone models
#' @param num.pops number of populations.
#' Default: \code{num.pops = 5}.
#' @param num.loci number of independent RADseq loci.
#' Default: \code{num.loci = 1000}.
#' @param div.time time since divergence of populations.
#' Default: \code{div.time = 25000}.
#' @param ne effective populations size.
#' Default: \code{ne = c(50, 500)}.
#' @param nm number of migrants per generation.
#' Default: \code{nm = c(0, 0.1, 0.5, 1, 5)}.
#' @param theta value of theta (= ne * migration rate).
#' Default: \code{theta = 0.2}.
#' @param mig.type migration topology type: \code{"island"} or 
#' \code{"stepping.stone"}.
#' @param num.reps number of replicates for each scenario.
#' Default: \code{num.reps = 10}.
#' @param num.rms.gens number of \code{rmetasim} generations to establish
#' linkage disequilibrium.
#' Default: \code{num.rms.gens = 5}.
#' @param label label for the output folder. With default \code{label = NULL},
#' the label will be "\code{sim.results.YYYYMMDD.HHMM}".
#' @param fsc.exec name of fastsimcoal executable.
#' Default: \code{fsc.exec = "fsc252"}.
#' @param parallel.core (optional) The number of core used for parallel
#' execution.
#' Default: \code{parallel::detectCores() - 1}.

#' @note The following arguments can be specified as single values or vectors:
#'  \code{num.pops, num.loci, div.time, ne, nm, theta, mig.type}. One simulation 
#'  scenario will be generated for each combination of unique values of these 
#'  arguments.  
#'  
#'  The function requires that the executable for \code{fastsimcoal} be installed 
#'  in a location where it can be executed from the current working directory.  
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
#'   
#' @return invisibly, the label used to name output files
#'
#' @author Eric Archer <eric.archer@@noaa.gov>
#' 
#' @importFrom strataG fscPopInfo fscLocusParams fscHistEv fastsimcoal alleleFreqs df2gtypes
#' @importFrom rmetasim landscape.simulate landscape.ploidy landscape.democol landscape.locusvec landscape.new.empty landscape.new.intparam landscape.new.floatparam landscape.new.switchparam landscape.new.local.demo landscape.new.epoch landscape.new.locus landscape.mig.matrix
#' @export
#' 
simulate_rad <- function(
  num.pops = 5,
  num.loci = 1000,
  div.time = 25000,
  ne = c(50, 500),
  nm = c(0, 0.1, 0.5, 1, 5),
  theta = 0.2,
  mig.type = c("island", "stepping.stone"),
  num.reps = 10,
  num.rms.gens = 5,
  label = NULL,
  fsc.exec = "fsc252",
  parallel.core = NULL
) {
  if (is.null(label)) {
    label <- paste0("sim.results.", format(Sys.time(), "%Y%m%d.%H%M"))
  }
  if (!dir.exists(label)) dir.create(label)
  
  num.pops <- sort(unique(num.pops))
  num.loci <- sort(unique(num.loci))
  div.time <- sort(unique(div.time))
  ne <- sort(unique(ne))
  nm <- sort(unique(nm))
  theta <- sort(unique(theta))
  mig.type <- sort(unique(mig.type))
  
  # create scenario data.frame
  sc.df <- expand.grid(
    num.pops = num.pops, num.loci = num.loci, div.time = div.time, ne = ne,
    ne = ne, nm = nm, theta = theta, mig.type = mig.type,
    stringsAsFactors = FALSE
  )
  sc.df <- cbind(scenario = 1:nrow(sc.df), sc.df)
  sc.df$mut.rate <- sc.df$theta / (4 * sc.df$ne)
  sc.df$mig.rate <- sc.df$nm / sc.df$ne
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
    fname <- file.path(label, paste("gtypes", i, "rdata", sep = "."))
    if (file.exists(fname)) next
    
    sc <- as.list(sc.df[i, ])
    sc$mig.mat <- sc$mig.mat[[1]]
    
    # run fastsimcoal
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
    rms.list <- lapply(1:length(fsc.list), function(i) {
      af <- strataG::alleleFreqs(fsc.list[[i]], by.strata = T)
      rl <- loadLandscape(sc, af, num.rms.gens)
      for (g in 1:num.rms.gens) {
        rl <- rmetasim::landscape.simulate(rl, 1)
        rl <- killExcess(rl, sc$ne)
      }
      landscape2gtypes(rl)
    })
    
    attr(rms.list, "scenario") <- attr(fsc.list, "scenario") <- sc.df[i, ]
    attr(rms.list, "label") <- attr(fsc.list, "label") <- label
    
    sim.data <- rms.list
    # save both fastsimcoal and rmetasim results to same workspace file
    save(fsc.list, sim.data, file = fname)
    fname
  })
  
  invisible(label)
}
