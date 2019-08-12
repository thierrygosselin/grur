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
#'  parameters for each scenario (\code{scenarios}), and one R  
#'  workspace file for each scenario, named \code{"gtypes.#.rdata"} where "#" 
#'  is the scenario number. Each scenario file contains a list named 
#'  \code{sim.data}, which contains a \code{strataG::gtypes} representation 
#'  of the simulated SNP data for every replicate in that scenario. The list 
#'  has a "\code{scenario}" attribute which is a one row data.frame providing 
#'  the parameters that generated the data and a "\code{label}" attribute which 
#'  is the label of the set of scenarios that were run. 

#' @return invisibly, the label used to name output files

#' @seealso \href{https://github.com/EricArcher/strataG}{strataG}

#' @author Eric Archer \email{eric.archer@@noaa.gov} and 
#'   Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
simulate_rad <- function(
  num.pops = 5,
  num.loci = 1000,
  div.time = 25000,
  ne = 50, #c(50, 500),
  nm = c(0, 5), #c(0, 0.1, 0.5, 1, 5),
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
    stop(
      "rmetasim needed for this function to work. ", 
      "Please follow the vignette for install instructions", 
      call. = FALSE
    )
  }
  if (!requireNamespace("strataG", quietly = TRUE)) {
    stop(
      "strataG needed for this function to work. ",
      "Install with install.packages('strataG')", 
      call. = FALSE
    )
  }
  
  start.time <- Sys.time()
  
  num.pops <- sort(unique(num.pops))
  num.loci <- sort(unique(num.loci))
  div.time <- sort(unique(div.time))
  ne <- sort(unique(ne))
  nm <- sort(unique(nm))
  theta <- sort(unique(theta))
  mig.type <- sort(unique(mig.type))
  
  if (is.null(label)) {
    sim.arguments <- paste0(
      "_pop.", paste0(num.pops, collapse = "."),
      "_loc.", paste0(num.loci, collapse = "."),
      "_div.", paste0(div.time, collapse = "."),
      "_ne.", paste0(ne, collapse = "."),
      "_nm.", paste0(nm, collapse = "."),
      "_the.", paste0(theta, collapse = "."),
      "_mig.", paste0(mig.type, collapse = "."),
      "_rep.", num.reps,
      "_gen.", num.rms.gens
    )
    label <- paste0(
      "sims", sim.arguments, "_", format(Sys.time(), "%Y%m%d@%H%M")
    )
  }
  if(dir.exists(label)) unlink(label, recursive = TRUE, force = TRUE)
  dir.create(label)
  
  # create scenario data.frame
  scenarios <- expand.grid(
    num.pops = num.pops, num.loci = num.loci, div.time = div.time, ne = ne,
    nm = nm, theta = theta, mig.type = mig.type,
    stringsAsFactors = FALSE
  )
  scenarios <- cbind(scenario = 1:nrow(scenarios), scenarios)
  scenarios$mut.rate <- scenarios$theta / (4 * scenarios$ne)
  scenarios$mig.rate <- scenarios$nm / scenarios$ne
  attr(scenarios, "label") <- label
  save(scenarios, file = file.path(label, paste0(label, "_scenarios.rdata")))
  
  # run scenarios
  for(sc.i in 1:nrow(scenarios)) {    
    fname <- file.path(label, paste0(label, "_genotypes_", sc.i, ".rdata"))
    if(file.exists(fname)) next
    
    cat(format(Sys.time()), "---- Scenario", sc.i, "----\n")
    sc <- as.list(scenarios[sc.i, ])
    sc$mig.mat <- make_mig_mat(sc$mig.rate, sc$num.pops, sc$mig.type)
    p <- run_fsc_sim(sc = sc, num.rep = num.reps)
    sim.list <- lapply(1:num.reps, function(sim.i) {
      fsc <- strataG::fscReadArp(p, sim = c(1, sim.i), drop.mono = TRUE) 
      rms <- fsc %>% 
        calc_freqs() %>% 
        run_rmetasim(sc = sc, num.gens = num.rms.gens) %>% 
        rmetasim::landscape.make.genind() %>% 
        strataG::genind2gtypes() %>% 
        strataG::as.data.frame()
      list(fsc = fsc, rms = rms)
    })
    strataG::fscCleanup(p$label)
    
    message("saving fastsimcoal and rmetasim results")
    save(sim.list, sc, label, file = fname)
    fname
  }
  
  timing <- difftime(Sys.time(), start.time)
  message("\nComputation time: ", format(round(timing, 2)))
  cat("#################### grur::simulate_rad completed #####################\n")
  options(width = opt.change)
  
  invisible(label)
}
