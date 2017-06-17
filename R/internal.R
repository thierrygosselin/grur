#' @keywords internal
#' 
landscape.new.ind.genos <- function(rland, PopulationSizes, AlleleFreqs = NULL) {
  if (!is.null(AlleleFreqs))
    if (length(AlleleFreqs) != length(rland$loci)) {
      stop("AlleleFreqs must be NULL or a list the same length as rland$loci")
    }
  if (length(PopulationSizes) != rland$intparam$stages * rland$intparam$habitats) {
    stop("PopulationSizes must be a vector of length equal to the product of #stages and #habitats")
  }
  
  loccols <- length(rmetasim::landscape.locusvec(rland))
  rland$individuals <- matrix(
    nrow = sum(PopulationSizes),
    ncol = rmetasim::landscape.democol() + length(rmetasim::landscape.locusvec(rland))
  )
  
  rland$intparam$nextid <- nrow(rland$individuals) + 1
  
  rland$individuals[, 1] <- inverse.rle(
    list(lengths = PopulationSizes, values = 1:length(PopulationSizes))
  ) #set the demographic stages
  rland$individuals[, 2] <- 0 #unused
  rland$individuals[, 3] <- rland$intparam$currentgen #time of birth
  rland$individuals[, 4] <- 1:nrow(rland$individuals) # ID
  rland$individuals[, 5] <- 0 #mother's id
  rland$individuals[, 6] <- 0 #father's id
  
  pops <- rland$individuals[, 1] %/% rland$intparam$stages + 1
  
  for (i in 1:length(AlleleFreqs)) {
    ainds <- 1:dim(AlleleFreqs[[i]])[1]
    frqs <- rowSums(AlleleFreqs[[i]][, "freq", , drop = FALSE])
    props <- frqs / sum(frqs)
    
    #make sure that the loci have the correct allele specification
    rland$loci[[i]]$alleles <- vector("list", length(props))
    for (j in ainds) {
      rland$loci[[i]]$alleles[[j]]$aindex <- j
      rland$loci[[i]]$alleles[[j]]$birth <- rland$intparam$currentgen
      rland$loci[[i]]$alleles[[j]]$prop <- props[j]
      rland$loci[[i]]$alleles[[j]]$state <- dimnames(AlleleFreqs[[i]])[[1]][j]
    }
    
    ###now specify genotypes based on within pop freqs
    cols <- rmetasim::landscape.democol() + which(rmetasim::landscape.locusvec(rland) == i)
    for (j in 1:dim(AlleleFreqs[[i]])[3]) {
      indx <- which(pops == j)
      num.als <- length(indx) * length(cols)
      rland$individuals[indx, cols] <- sample(ainds, num.als, T, AlleleFreqs[[i]][, "prop", j])
    }
  }
  
  rland
}

#' @keywords internal
#' 
loadLandscape <- function(sc, AlleleFreqs, num.gens) {
  rl <- rmetasim::landscape.new.intparam(
    rmetasim::landscape.new.empty(), h = sc$num.pops, 
    s = 2, cg = 0, ce = 0, totgen = num.gens + 1
  )
  
  rl <- rmetasim::landscape.new.floatparam(rmetasim::landscape.new.switchparam(rl))
  
  localS <- matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2)
  localR <- matrix(c(0, 0, 1.2, 0), nrow = 2, ncol = 2)
  localM <- matrix(c(0, 0, 0, 1.2), nrow = 2, ncol = 2)
  rl <- rmetasim::landscape.new.local.demo(rl, localS, localR, localM)
  
  S <- M <- matrix(0, nrow = sc$num.pops * 2, ncol = sc$num.pops * 2)
  diag(S) <- diag(M) <- 1
  R <- rmetasim::landscape.mig.matrix(
    h = nrow(sc$mig.mat), s = 2, mig.model = "custom", R.custom = sc$mig.mat
  )$R
  rl <- rmetasim::landscape.new.epoch(rl, R = R, carry = rep(sc$Ne, sc$num.pops))
  
  #just make loci that have the correct type and ploidy
  for (i in 1:length(AlleleFreqs)) {
    rl <- rmetasim::landscape.new.locus(
      rl, type = 2, ploidy = 2, mutationrate = 0, #sc$mut.rate,
      transmission = 0, numalleles = 2, states = NULL
    )
  }
  
  landscape.new.ind.genos(rl, rep(c(sc$Ne, 0), sc$num.pops), AlleleFreqs)
}

#' @keywords internal
#' 
landscape2gtypes <- function(Rland) {
  pl <- rmetasim::landscape.ploidy(Rland)
  strata <- Rland$individuals[, 1] %/% Rland$intparam$stages + 1
  gen.data <- Rland$individuals[, -(1:rmetasim::landscape.democol())]
  rownames(gen.data) <- Rland$individuals[, 4]
  loc.names <- paste0("Locus", 1:length(pl))
  colnames(gen.data) <- paste(rep(loc.names, each = pl[1]), 1:pl[1], sep = ".")
  gen.data <- cbind(strata = strata, gen.data)
  strataG::df2gtypes(gen.data, ploidy = pl[1], id.col = NULL, strata.col = 1, loc.col = 2)
}

#' @keywords internal
#' 
killExcess <- function(rl, n) {
  to.kill <- tapply(1:nrow(rl$individuals), rl$individuals[, 1], function(i) {
    if (length(i) > n) {
      sample(i, length(i) - n)
    } else NULL
  })
  to.kill <- unlist(unname(to.kill))
  if (length(to.kill) > 0) rl$individuals <- rl$individuals[-to.kill, ]
  rl
}
