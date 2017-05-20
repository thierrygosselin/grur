#' @title Simulate SNPs
#' @description Simulate populations of SNP data following island or stepping 
#'   stone models
#'   
#' @param num.pops number of populations
#' @param num.loci number of independent SNP loci
#' @param div.time time since divergence of populations
#' @param Ne effective populations size
#' @param Nm number of migrants per generation
#' @param theta value of theta (= Ne * migration rate)
#' @param mig.type migration topology type: \code{"island"} or 
#'   \code{"stepping.stone"}
#' @param num.reps number of replicates for each scenario
#' @param num.rms.gens number of \code{rmetasim} generations to establish 
#'   linkage disequilibrium
#' @param label label for the output folder. If \code{NULL}, the label will be 
#'   "\code{sim.results.YYYYMMDD.HHMM}".
#' @param fsc.exec name of fastsimcoal executable
#' @param num.cores number of cores to use
#'   
#' @note The following arguments can be specified as single values or vectors:
#'  \code{num.pops, num.loci, div.time, Ne, Nm, theta, mig.type}. One simulation 
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
simulateSNPs <- function(
  num.pops = 5,
  num.loci = 1000,
  div.time = 25000,
  Ne = c(50, 500),
  Nm = c(0, 0.1, 0.5, 1, 5),
  theta = 0.2,
  mig.type = c("island", "stepping.stone"),
  num.reps = 10,
  num.rms.gens = 5,
  label = NULL,
  fsc.exec = "fsc252",
  num.cores = NULL
) {
  if(is.null(label)) {
    label <- paste0("sim.results.", format(Sys.time(), "%Y%m%d.%H%M"))
  }
  if(!dir.exists(label)) dir.create(label)
  
  num.pops <- sort(unique(num.pops))
  num.loci <- sort(unique(num.loci))
  div.time <- sort(unique(div.time))
  Ne <- sort(unique(Ne))
  Nm <- sort(unique(Nm))
  theta <- sort(unique(theta))
  mig.type <- sort(unique(mig.type))
  
  # create scenario data.frame
  sc.df <- expand.grid(
    num.pops = num.pops, num.loci = num.loci, div.time = div.time, Ne = Ne,
    Ne = Ne, Nm = Nm, theta = theta, mig.type = mig.type,
    stringsAsFactors = FALSE
  )
  sc.df <- cbind(scenario = 1:nrow(sc.df), sc.df)
  sc.df$mut.rate <- sc.df$theta / (4 * sc.df$Ne)
  sc.df$mig.rate <- sc.df$Nm / sc.df$Ne
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
        for(k in 1:(num.pops - 1)) {
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
    if(file.exists(fname)) next
    
    sc <- as.list(sc.df[i, ])
    sc$mig.mat <- sc$mig.mat[[1]]
    
    # run fastsimcoal
    fsc.list <- with(sc, {
      n <- num.pops
      pi <- fscPopInfo(pop.size = rep(Ne, n), sample.size = rep(Ne, n))
      lp <- fscLocusParams(locus.type = "snp", num.loci = num.loci, mut.rate = mut.rate)
      he <- fscHistEv(num.gen = rep(div.time, n - 1), source.deme = 1:(n - 1))
      lapply(1:num.reps, function(rep) {
        lbl <- paste0("scenario_", i, ".replicate_", rep)
        fastsimcoal(
          pi, lp, mig.rates = mig.mat, hist.ev = he, label = lbl, 
          quiet = FALSE, exec = fsc.exec, num.cores = num.cores
        )
      })
    })
    
    # run rmetasim for 'num.gens' generations using fastsimcoal runs as initialization
    rms.list <- lapply(1:length(fsc.list), function(i) {
      af <- alleleFreqs(fsc.list[[i]], by.strata = T)
      rl <- loadLandscape(sc, af, num.rms.gens)
      for(g in 1:num.rms.gens) {
        rl <- landscape.simulate(rl, 1)
        rl <- killExcess(rl, sc$Ne)
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
