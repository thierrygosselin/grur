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
#' @param use.wd (optional, logical) Use the working directory for output files? 
#' Default: \code{use.wd = FALSE} means that files are written to temporary folder.

#' @examples 
#' \dontrun{
#' library(grur)
#' sim <- grur::simulate_rad(
#'   num.pops = 3, 
#'   num.loci = 200, 
#'   div.time = 1000, 
#'   ne = c(50,200),
#'   nm = c(0, 0.1), 
#'   theta = 0.2, 
#'   mig.type = c("island", "stepping.stone"), 
#'   num.reps = 3,
#'   num.rms.gens = 5,
#'   label = NULL, 
#'   fsc.exec = "fsc26", 
#'   parallel.core = 4
#'  )
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

#' @return invisibly, a list containing the parameters used to run the 
#'   simulation and the filenames of the genotypes for each scenario
#'   and replicate

#' @seealso \href{https://github.com/EricArcher/strataG}{strataG}

#' @author Eric Archer \email{eric.archer@@noaa.gov} and 
#'   Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
simulate_rad <- function(
  num.pops = 5,
  num.loci = 1000,
  div.time = 25000,
  ne = 50,
  nm = c(0, 5), 
  theta = 0.2,
  mig.type = c("island", "stepping.stone"),
  num.reps = 10,
  num.rms.gens = 10,
  label = NULL,
  fsc.exec = "fsc26",
  parallel.core = parallel::detectCores() - 1,
  use.wd = FALSE
) {
  opt.change <- getOption("width")
  options(width = 70)
  cat("#######################################################################\n")
  cat("######################### grur::simulate_rad ##########################\n")
  cat("#######################################################################\n")
  
  # # Check if rmetasim installed ------------------------------------------------
  # if (!requireNamespace("rmetasim", quietly = TRUE)) {
  #   stop(
  #     "rmetasim needed for this function to work. ", 
  #     "Please follow the vignette for install instructions", 
  #     call. = FALSE
  #   )
  # }
  # if (!requireNamespace("strataG", quietly = TRUE)) {
  #   stop(
  #     "strataG needed for this function to work. ",
  #     "Install with devtools::install_github('ericarcher/strataG')", 
  #     call. = FALSE
  #   )
  # }
  
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
      "_pop.", paste0(num.pops, collapse = ","),
      "_loc.", paste0(num.loci, collapse = ","),
      "_div.", paste0(div.time, collapse = ","),
      "_ne.", paste0(ne, collapse = ","),
      "_nm.", paste0(nm, collapse = ","),
      "_the.", paste0(theta, collapse = ","),
      "_mig.", paste0(mig.type, collapse = ","),
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
    nm = nm, theta = theta, mig.type = mig.type, num.rms.gens = num.rms.gens,
    stringsAsFactors = FALSE
  )
  scenarios <- cbind(scenario = 1:nrow(scenarios), scenarios)
  scenarios$mut.rate <- scenarios$theta / (4 * scenarios$ne)
  scenarios$mig.rate <- scenarios$nm / scenarios$ne
  
  # create parameter list
  params <- list(
    label = label,
    scenarios = scenarios,
    num.rep = num.reps,
    fsc.exec = fsc.exec,
    num.cores = parallel.core,
    use.wd = use.wd
  )
  params$rep.df <- expand.grid( 
    rep = 1:num.reps, 
    sc = 1:nrow(scenarios),
    stringsAsFactors = FALSE
  )[, c("sc", "rep")]
  
  # preload rmetasim landscapes
  params$Rland <- lapply(1:nrow(params$scenarios), .setupScRland, params = params)

  # run replicates
  sc.rep.vec <- 1:nrow(params$rep.df)
  fnames <- if(params$num.cores == 1) {  
    tryCatch(lapply(sc.rep.vec, .runWithLabel, params = params))
  } else {
    grur_future(.x = sc.rep.vec, .f = .runScRep, flat.future = "walk", parallel.core = params$num.cores, params = params)
  }
  
  timing <- difftime(Sys.time(), start.time)
  message("\nComputation time: ", format(round(timing, 2)))
  cat("#################### grur::simulate_rad completed #####################\n")
  options(width = opt.change)
  
  params <- list(
    label = params$label,
    scenarios = params$scenarios,
    num.rep = params$num.rep,
    replicate.fnames = stats::setNames(
      cbind(params$rep.df, unlist(fnames), stringsAsFactors = FALSE),
      c("scenario", "replicate", "fname")
    )
  )
  save(
    params, 
    file = file.path(params$label, paste0(params$label, "_params.rdata"))
  )
  invisible(params)
}


#' @noRd
#' 
.runWithLabel <- function(rep.i, params) {
  cat(paste0(
    format(Sys.time()), 
    " ---- Scenario ", 
    params$scenarios$scenario[params$rep.df$sc[rep.i]], 
    ", Replicate ", 
    params$rep.df$rep[rep.i],
    " (", 
    round(100 * rep.i / nrow(params$rep.df)), 
    "%) ----\n"
  ))
  .runScRep(rep.i, params)
}


#' @noRd
#' 
.runScRep <- carrier::crate(function(rep.i, params) {
  `%>%` <- magrittr::`%>%`
  sc.num <- params$rep.df$sc[rep.i]
  rep.num <- params$rep.df$rep[rep.i]
  
  tryCatch(
    {
      p <- grur::.runFscSim(rep.i, params)
      gen.data <- strataG::fscReadArp(p)
      strataG::fscCleanup(p$label, p$folder)
      sc <- params$scenarios[sc.num, ]
      if (sc$num.rms.gen > 0) {
        cat(format(Sys.time()), "running rmetasim...\n")
        gen.data <- gen.data %>%
          grur::.calcFreqs() %>%
          grur::.runRmetasim(Rland = params$Rland[[sc.num]], sc = sc) %>%
          strataG::landscape2df()
      }
      
      fname <- file.path(
        params$label,
        paste0(params$label, "_scenario.", sc$scenario, "_replicate.", rep.num, ".csv")
      )
      utils::write.csv(gen.data, file = fname, row.names = FALSE)
      fname
    },
    error = function(e) {
      stop(
        format(Sys.time()),
        " Scenario ", params$scenarios$scenario[sc.num],
        ", Replicate ", rep.num, 
        ": ", e
      )
    }
  )
})

#' @noRd
#' @export
#' @keywords internal
.runFscSim <- function(rep.i, params) {
  sc.i <- params$rep.df$sc[rep.i]
  sc <- params$scenarios[sc.i, ]
  
  deme.list <- lapply(1:sc$num.pops, function(i) {
    strataG::fscDeme(deme.size = sc$ne, sample.size = sc$ne)
  })
  deme.list$ploidy <- 2
  
  strataG::fscWrite(
    demes = do.call(strataG::fscSettingsDemes, deme.list),
    migration = if(sc$num.pops > 1) {
      mig.mat <- .makeMigMat(sc$mig.rate, sc$num.pops, sc$mig.type) 
      strataG::fscSettingsMigration(mig.mat)
    } else NULL,
    events = .makeEventSettings(sc$div.time, sc$num.pops),
    genetics = strataG::fscSettingsGenetics(
      strataG::fscBlock_snp(1, sc$mut.rate), 
      num.chrom = sc$num.loci
    ),
    label = paste0(
      params$label, ".sc_", sc$scenario, ".rep_", params$rep.df$rep[rep.i]
    ),
    use.wd = params$use.wd
  ) %>% 
    strataG::fscRun(num.cores = 1, exec = params$fsc.exec)
}

#' @noRd
#' 
.setupScRland <- function(sc.num, params) {  
  if(params$scenarios$num.rms.gen[sc.num] == 0) return(NA)
  sc <- params$scenarios[sc.num, ]
  localS <- matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2)
  localR <- matrix(c(0, 0, 1, 0), nrow = 2, ncol = 2)
  localM <- matrix(c(0, 0, 0, 1), nrow = 2, ncol = 2)
  S <- M <- matrix(0, nrow = sc$num.pops * 2, ncol = sc$num.pops * 2)
  diag(S) <- diag(M) <- 1
  R <- if(sc$num.pops == 1) NULL else {
    rmetasim::landscape.mig.matrix(
      h = sc$num.pops, s = 2, mig.model = "custom", 
      R.custom = .makeMigMat(sc$mig.rate, sc$num.pops, sc$mig.type)
    )$R
  }
  
  rmetasim::landscape.new.empty() %>% 
    rmetasim::landscape.new.intparam(
      h = sc$num.pops, s = 2, cg = 0, ce = 0, totgen = sc$num.rms.gens + 1
    ) %>% 
    rmetasim::landscape.new.switchparam() %>% 
    rmetasim::landscape.new.floatparam() %>% 
    rmetasim::landscape.new.local.demo(localS, localR, localM) %>% 
    rmetasim::landscape.new.epoch(R = R, carry = rep(sc$ne, sc$num.pops)) 
}

#' @noRd
#' @export
#' @keywords internal
.runRmetasim <- function(freqs, Rland, sc) {
  for (i in 1:length(freqs$global)) {
    Rland <- rmetasim::landscape.new.locus(
      Rland, type = 2, ploidy = 2, mutationrate = 0,
      transmission = 0, numalleles = 2, allelesize = 1,
      frequencies = freqs$global[[i]], states = names(freqs$global[[i]])
    )
  }
  
  Rland %>% 
    rmetasim::landscape.new.individuals(rep(c(sc$ne, 0), sc$num.pops)) %>% 
    rmetasim::landscape.setpopfreq(freqs$pop) %>% 
    rmetasim::landscape.simulate(numit = sc$num.rms.gens)
}

#' @noRd
#' @export
#' @keywords internal
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' 
.calcFreqs <- function(snps) {
  if (ncol(snps) < 3) stop("no loci present")
  mac.df <- snps %>% 
    strataG::df2gtypes(ploidy = 2) %>% 
    strataG::as.data.frame(coded = T) %>% 
    dplyr::select(-.data$id) %>% 
    grur::rad_long(
      x = ., 
      cols = .data$stratum,
      names_to = "locus",
      values_to = "mac"
    ) %>% 
    dplyr::mutate(
      stratum = as.numeric(factor(.data$stratum)),
      locus = as.numeric(factor(.data$locus))
    )
  
  list(
    global = lapply(split(mac.df, mac.df$locus), function(loc.df) {
      .alleleProp(loc.df$mac)
    }),
    pop = lapply(split(mac.df, mac.df$stratum), function(st.df) {
      lapply(split(st.df, st.df$locus), function(loc.df) .alleleProp(loc.df$mac))
    })
  )
}

#' @noRd
#' 
.makeMigMat <- function(mig.rate, num.pops, 
                       type = c("island", "stepping.stone")) {
  if(is.na(mig.rate)) return(NULL)
  type <- match.arg(type)
  mig.mat <- switch(
    type,      
    island = {
      m <- mig.rate / (num.pops - 1)
      matrix(rep(m, num.pops ^ 2), nrow = num.pops)
    },
    stepping.stone = {
      mat <- matrix(0, nrow = num.pops, ncol = num.pops)
      m <- mig.rate / 2
      for (k in 1:(num.pops - 1)) {
        mat[k, k + 1] <- mat[k + 1, k] <- m
      }
      mat[1, num.pops] <- mat[num.pops, 1] <- m
      mat
    }
  )
  diag(mig.mat) <- 1 - mig.rate
  mig.mat
}

#' @noRd
#' 
.alleleProp <- function(mac) {
  maf <- mean(mac == 0) + (mean(mac == 1) / 2)
  c('1' = maf, '2' = 1 - maf)
}

#' @noRd
#' 
.makeEventSettings <- function(dvgnc.time, num.pops) {
  if(num.pops == 1) return(NULL)
  pop.pairs <- t(utils::combn(num.pops, 2) - 1)
  pop.pairs <- pop.pairs[pop.pairs[, 1] == 0, , drop = FALSE]
  do.call(
    strataG::fscSettingsEvents, 
    lapply(
      1:nrow(pop.pairs),
      function(i) {
        strataG::fscEvent(dvgnc.time, pop.pairs[i, 2], pop.pairs[i, 1])
      }
    )
  )
}
