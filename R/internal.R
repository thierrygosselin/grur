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

#' @title ind_genotyped_helper
#' @description Help individual's genotyped threshold
#' @rdname ind_genotyped_helper
#' @keywords internal
ind_genotyped_helper <- function(x) {
  # x <- res$missing.genotypes.ind
  # Set the breaks for the figure
  max.ind <- dplyr::n_distinct(x$INDIVIDUALS)
  
  
  threshold.helper.overall <- x %>%
    dplyr::summarise(
      `0` = length(PERCENT[PERCENT >= 0]),
      `10` = length(PERCENT[PERCENT >= 10]),
      `20` = length(PERCENT[PERCENT >= 20]),
      `30` = length(PERCENT[PERCENT >= 30]),
      `40` = length(PERCENT[PERCENT >= 40]),
      `50` = length(PERCENT[PERCENT >= 50]),
      `60` = length(PERCENT[PERCENT >= 60]),
      `70` = length(PERCENT[PERCENT >= 70]),
      `80` = length(PERCENT[PERCENT >= 80]),
      `90` = length(PERCENT[PERCENT >= 90]),
      `100` = length(PERCENT[PERCENT == 100])
    ) %>%
    tidyr::gather(data = ., key = GENOTYPED_THRESHOLD, value = NUMBER_INDIVIDUALS) %>%
    dplyr::mutate(POP_ID = rep("OVERALL", n()))
  
  threshold.helper.pop <- x %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::summarise(
      `0` = length(PERCENT[PERCENT >= 0]),
      `10` = length(PERCENT[PERCENT >= 10]),
      `20` = length(PERCENT[PERCENT >= 20]),
      `30` = length(PERCENT[PERCENT >= 30]),
      `40` = length(PERCENT[PERCENT >= 40]),
      `50` = length(PERCENT[PERCENT >= 50]),
      `60` = length(PERCENT[PERCENT >= 60]),
      `70` = length(PERCENT[PERCENT >= 70]),
      `80` = length(PERCENT[PERCENT >= 80]),
      `90` = length(PERCENT[PERCENT >= 90]),
      `100` = length(PERCENT[PERCENT == 100])
    ) %>%
    tidyr::gather(data = ., key = GENOTYPED_THRESHOLD, value = NUMBER_INDIVIDUALS, -POP_ID)
  
  mean.pop <- threshold.helper.pop %>%
    dplyr::group_by(GENOTYPED_THRESHOLD) %>%
    dplyr::summarise(
      NUMBER_INDIVIDUALS = round(mean(NUMBER_INDIVIDUALS), 0)
    ) %>%
    dplyr::mutate(POP_ID = rep("MEAN_POP", n()))
  
  threshold.helper <- suppressWarnings(
    dplyr::bind_rows(threshold.helper.pop, mean.pop, threshold.helper.overall) %>%
      dplyr::mutate(
        GENOTYPED_THRESHOLD = as.numeric(GENOTYPED_THRESHOLD),
        POP_ID = factor(POP_ID, levels = c(levels(x$POP_ID), "MEAN_POP", "OVERALL"), ordered = TRUE)
      ))
  threshold.helper.pop <- mean.pop <- threshold.helper.overall <- x <- NULL
  
  
  
  #Function to replace plyr::round_any
  rounder <- function(x, accuracy, f = round) {
    f(x / accuracy) * accuracy
  }
  
  if (max.ind >= 1000) {
    y.breaks.by <- rounder(max.ind/10, 100, ceiling)
    y.breaks.max <- rounder(max.ind, 1000, ceiling)
    y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
  } else {
    y.breaks.by <- rounder(max.ind/10, 10, ceiling)
    y.breaks.max <- rounder(max.ind, 100, ceiling)
    y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
  }
  
  axis.title.element.text.fig <- ggplot2::element_text(
    size = 12, family = "Helvetica", face = "bold")
  axis.text.element.text.fig <- ggplot2::element_text(
    size = 10, family = "Helvetica")
  
  plot.ind.geno.threshold <- ggplot2::ggplot(
    threshold.helper,
    ggplot2::aes(x = GENOTYPED_THRESHOLD, y = NUMBER_INDIVIDUALS)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
    ggplot2::scale_x_continuous(name = "Individual's missing genotyped threshold (percent)") +
    ggplot2::scale_y_continuous(name = "Individuals\n(blacklisted number)", breaks = y.breaks, limits = c(0, y.breaks.max)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = axis.title.element.text.fig,
      axis.title.y = axis.title.element.text.fig,
      axis.text.x = axis.text.element.text.fig,
      axis.text.y = axis.text.element.text.fig
    ) +
    ggplot2::facet_grid(~POP_ID)
  # plot.ind.geno.threshold
  return(plot.ind.geno.threshold)
}#End ind_genotyped_helper

#' @title blacklists_id_generator
#' @description Generate blacklist of ids
#' @rdname blacklists_id_generator
#' @keywords internal
blacklists_id_generator <- function(x, y, path.folder) {
  blacklist <- list()
  blacklist.id.missing.geno <- y %>%
    dplyr::filter(PERCENT >= x) %>%
    dplyr::ungroup(.) %>%
    dplyr::select(INDIVIDUALS)
  if (length(blacklist.id.missing.geno$INDIVIDUALS) > 0) {
    blacklist.name <- stringi::stri_join("blacklist.id.missing.", x)
    readr::write_tsv(
      blacklist.id.missing.geno,
      stringi::stri_join(path.folder, "/", as.name(blacklist.name), ".tsv"))
    blacklist[[blacklist.name]] <- blacklist.id.missing.geno
  } else {
    blacklist <- NULL
  }
  return(blacklist)
}#End blacklists_id_generator


#' @title whitelists_markers_generator
#' @description Generate whitelists of markers
#' @rdname whitelists_markers_generator
#' @keywords internal
whitelists_markers_generator <- function(x, y, path.folder) {
  whitelist <- list()
  tidy.col <- colnames(y)
  markers.meta <- purrr::keep(
    .x = tidy.col,
    .p = tidy.col %in% c("MARKERS", "CHROM", "LOCUS", "POS"))
  
  whitelist.missing.geno <- dplyr::ungroup(y) %>%
    dplyr::filter(MISSING_GENOTYPE_PROP <= x) %>%
    dplyr::select(dplyr::one_of(markers.meta)) %>% 
    dplyr::distinct(MARKERS, .keep_all = TRUE)
  
  if (length(whitelist.missing.geno$MARKERS) > 0) {
    whitelist.name <- stringi::stri_join("whitelist.markers.missing.max.", x)
    readr::write_tsv(
      whitelist.missing.geno,
      stringi::stri_join(path.folder, "/", as.name(whitelist.name), ".tsv"))
    whitelist[[whitelist.name]] <- whitelist.missing.geno
  } else {
    whitelist <- NULL
  }
  return(whitelist)
}#End whitelists_markers_generator


#' @title markers_genotyped_helper
#' @description Help individual's genotyped threshold
#' @rdname markers_genotyped_helper
#' @keywords internal
markers_genotyped_helper <- function(x, y) {
  # x <- res$missing.genotypes.markers.pop
  # Set the breaks for the figure
  max.markers <- dplyr::n_distinct(x$MARKERS)
  
  threshold.helper.overall <- y %>% 
    dplyr::ungroup(.) %>% 
    dplyr::summarise(
      `0` = length(PERCENT[PERCENT >= 0]),
      `10` = length(PERCENT[PERCENT >= 10]),
      `20` = length(PERCENT[PERCENT >= 20]),
      `30` = length(PERCENT[PERCENT >= 30]),
      `40` = length(PERCENT[PERCENT >= 40]),
      `50` = length(PERCENT[PERCENT >= 50]),
      `60` = length(PERCENT[PERCENT >= 60]),
      `70` = length(PERCENT[PERCENT >= 70]),
      `80` = length(PERCENT[PERCENT >= 80]),
      `90` = length(PERCENT[PERCENT >= 90]),
      `100` = length(PERCENT[PERCENT == 100])
    ) %>%
    tidyr::gather(data = ., key = GENOTYPED_THRESHOLD, value = NUMBER_MARKERS) %>%
    dplyr::mutate(POP_ID = rep("OVERALL", n()))
  
  threshold.helper.pop <- x %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::summarise(
      `0` = length(PERCENT[PERCENT >= 0]),
      `10` = length(PERCENT[PERCENT >= 10]),
      `20` = length(PERCENT[PERCENT >= 20]),
      `30` = length(PERCENT[PERCENT >= 30]),
      `40` = length(PERCENT[PERCENT >= 40]),
      `50` = length(PERCENT[PERCENT >= 50]),
      `60` = length(PERCENT[PERCENT >= 60]),
      `70` = length(PERCENT[PERCENT >= 70]),
      `80` = length(PERCENT[PERCENT >= 80]),
      `90` = length(PERCENT[PERCENT >= 90]),
      `100` = length(PERCENT[PERCENT == 100])
    ) %>%
    tidyr::gather(data = ., key = GENOTYPED_THRESHOLD, value = NUMBER_MARKERS, -POP_ID)
  
  mean.pop <- threshold.helper.pop %>%
    dplyr::group_by(GENOTYPED_THRESHOLD) %>%
    dplyr::summarise(
      NUMBER_MARKERS = round(mean(NUMBER_MARKERS), 0)
    ) %>%
    dplyr::mutate(POP_ID = rep("MEAN_POP", n()))
  
  threshold.helper <- suppressWarnings(
    dplyr::bind_rows(threshold.helper.pop, mean.pop, threshold.helper.overall) %>%
      dplyr::mutate(
        GENOTYPED_THRESHOLD = as.numeric(GENOTYPED_THRESHOLD),
        POP_ID = factor(POP_ID, levels = c(levels(x$POP_ID), "MEAN_POP", "OVERALL"), ordered = TRUE)
      ))
  threshold.helper.pop <- mean.pop <- threshold.helper.overall <- x <- y <- NULL
  
  
  
  #Function to replace plyr::round_any
  rounder <- function(x, accuracy, f = round) {
    f(x / accuracy) * accuracy
  }
  
  if (max.markers >= 1000) {
    y.breaks.by <- rounder(max.markers / 10, 100, ceiling)
    y.breaks.max <- rounder(max.markers, 1000, ceiling)
    y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
  } else {
    y.breaks.by <- rounder(max.markers / 10, 10, ceiling)
    y.breaks.max <- rounder(max.markers, 100, ceiling)
    y.breaks <- seq(0, y.breaks.max, by = y.breaks.by)
  }
  
  
  plot.markers.geno.threshold <- ggplot2::ggplot(
    threshold.helper,
    ggplot2::aes(x = GENOTYPED_THRESHOLD, y = NUMBER_MARKERS)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
    ggplot2::scale_x_continuous(name = "Marker's missing genotyped threshold (percent)") +
    ggplot2::scale_y_continuous(name = "Markers\n(whitelisted number)", breaks = y.breaks, limits = c(0, y.breaks.max)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica"),#, angle = 90, hjust = 1, vjust = 0.5),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    ) +
    ggplot2::facet_grid(~POP_ID)
  # plot.markers.geno.threshold
  return(plot.markers.geno.threshold)
}#End markers_genotyped_helper


#' @title generate_pcoa_plot
#' @description Generate the PCoA plots
#' @rdname generate_pcoa_plot
#' @keywords internal
#' @importFrom ggpubr ggarrange
generate_pcoa_plot <- function(
  strata.select,
  pc.to.do,
  vectors,
  variance.component,
  path.folder,
  write.plot
) {
  pcoa.plots <- list()
  
  pcoa_plot <- function(
    pc.to.do,
    vectors,
    variance.component) {
    pcoa.plots <- list()
    
    # pc.to.do <- pc.to.do[[1]]
    # vectors <- res$vectors
    # strata.select <- "POP_ID"
    
    pcx <- pc.to.do[1]
    pcy <- pc.to.do[2]
    element.text.fig <- ggplot2::element_text(
      size = 12, family = "Helvetica", face = "bold")
    
    # size = MISSING_GENOTYPE_PERCENT,
    # ggplot2::scale_size_area(name = "Individual's missing\ngenotypes proportion", max_size = 5) +
      
    ibm.plot <- ggplot2::ggplot(
      vectors,
      ggplot2::aes_string(
        x = stringi::stri_join("Axis.", pcx),
        y = stringi::stri_join("Axis.", pcy), size = vectors$MISSING_GENOTYPE_PERCENT),
      environment = environment()) +
      ggplot2::geom_point(ggplot2::aes_string(colour = strata.select), alpha = 0.5) +
      ggplot2::labs(x = stringi::stri_join("PCo", pcx, " [", variance.component[pcx,2], "]")) +
      ggplot2::labs(y = stringi::stri_join("PCo", pcy, " [", variance.component[pcy,2], "]")) +
      ggplot2::scale_size_area(name = "Individual's\nmissing genotypes\n(proportion)", max_size = 5) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title.x = element.text.fig,
        axis.title.y = element.text.fig,
        legend.title = element.text.fig,
        legend.text = element.text.fig
      )
    ibm_plot_name <- stringi::stri_join(
      "ibm.plot.pco", pcx, ".pco", pcy, ".strata.", strata.select)
    pcoa.plots[[ibm_plot_name]] <- ibm.plot
    return(pcoa.plots)
  }#End pcoa_plot
  
  pcoa.plots.strata <- purrr::map(
    .x = pc.to.do, .f = pcoa_plot,
    vectors = vectors,
    variance.component = variance.component) %>%
    purrr::flatten(.)
  
  # pcoa.plots.strata <- arrange_plots_legend(
  #   pcoa.plots.strata, ncol = 2, nrow = 3,
  #   position = "right")
  # 
  # ggplot2::labs(title = stringi::stri_join("Principal Coordinates Analysis (PCoA)\n Identity by Missing (IBM) with strata = ", i)) +
  
  pcoa.plots.strata <- ggpubr::ggarrange(
    pcoa.plots.strata[[1]],
    pcoa.plots.strata[[2]],
    pcoa.plots.strata[[3]],
    pcoa.plots.strata[[4]],
    pcoa.plots.strata[[5]],
    pcoa.plots.strata[[6]],
    ncol = 2, nrow = 3, legend = "right", common.legend = TRUE
  )
  
  if (write.plot) {
    ggplot2::ggsave(
      filename = stringi::stri_join(path.folder, "/ibm.plots.strata.", strata.select, ".pdf"),
      plot = pcoa.plots.strata,
      width = 20, height = 15,
      dpi = 600, units = "cm",
      useDingbats = FALSE)
  }
  
  ibm_strata_name <- stringi::stri_join("ibm.strata.", strata.select)
  pcoa.plots[[ibm_strata_name]] <- pcoa.plots.strata
  
  return(pcoa.plots)
}#End generate_pcoa_plot


# tested alternative to ggpubr or cowplot package to reduce installed packages...
# the problem is my limited understanding of grid and grid extra 
# note to myself: you managed to write the combined plots but you failed in
# keeping the combined plots in an object that can be return in the result list

# arrange_plots_legend <- function(plots, ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
#   # plots <- list(...)
#   # plots <- unlist(plots)
#   position <- match.arg(position)
#   
#   g <- ggplot2::ggplotGrob(plots[[1]] + ggplot2::theme(legend.position = position))$grobs
#   legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
#   lheight <- sum(legend$height)
#   lwidth <- sum(legend$width)
#   
#   # gridExtra::grid.arrange(
#   #   do.call(gridExtra::arrangeGrob, lapply(plots, function(x)
#   #     x + ggplot2::theme(legend.position="none"))),
#   #   legend,
#   #   ncol = 1,
#   #   heights = grid::unit.c(grid::unit(1, "npc") - lheight, lheight))
#   gl <- lapply(plots, function(x) x + ggplot2::theme(legend.position="none"))
#   gl <- c(gl, ncol = ncol, nrow = nrow)
#   combined <- switch(
#     position,
#     "bottom" = gridExtra::arrangeGrob(
#       do.call(gridExtra::arrangeGrob, gl),
#       legend,
#       ncol = 1,
#       heights = grid::unit.c(grid::unit(1, "npc") - lheight, lheight)),
#     "right" = gridExtra::arrangeGrob(
#       do.call(gridExtra::arrangeGrob, gl),
#       legend,
#       ncol = 2,
#       widths = grid::unit.c(grid::unit(1, "npc") - lwidth, lwidth)))
#   
#   grid::grid.newpage()
#   grid::grid.draw(combined)
#   return(combined)
# }


#' @title pct_missing_by_total
#' @description Generates plot missing by total
#' @rdname pct_missing_by_total
#' @importFrom rlang .data UQ ":="
#' @importFrom dplyr group_by left_join summarise mutate filter ungroup select arrange
#' @importFrom ggplot2 ggplot aes_string geom_segment geom_point geom_abline geom_hline aes scale_x_log10 facet_wrap labs
#' @importFrom purrr map
#' @importFrom tidyr nest unnest
#' @importFrom broom tidy
#' @importFrom stats lm
#' @keywords internal

pct_missing_by_total <- function(strata.select, data, ci = 0.95, path.folder, write.plot = TRUE) {
  res <- list() # to store results
  
  data <- dplyr::rename(
    .data = data,
    STRATA_SELECT = .data[[rlang::UQ(strata.select)]])
  
  # the expected % missing by random (1 / number of factor levels)
  ran.pct.miss <- 1 / length(unique(data$STRATA_SELECT))
  
  # convert GT_BIN to missing or not
  data$is.missing <- data$GT_MISSING_BINARY == 0
  
  
  # summarize missingness by locus and factor column
  miss.smry <- dplyr::group_by(.data = data, MARKERS, STRATA_SELECT) %>% 
    dplyr::summarise(num.missing.col = sum(is.missing)) %>%
    dplyr::left_join(
      data %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::summarise(num.missing.total = sum(is.missing)),
      by = "MARKERS"
    ) %>%
    dplyr::filter(num.missing.total > 0) %>%
    dplyr::mutate(
      pct.missing = num.missing.col / num.missing.total
    ) %>%
    dplyr::ungroup(.)
  
  # summarize overall percent missing by factor column
  pct.miss.col <- miss.smry %>%
    dplyr::group_by(STRATA_SELECT) %>%
    dplyr::summarise(pct.missing = sum(num.missing.col) / sum(num.missing.total)) %>%
    dplyr::ungroup(.) %>% 
    dplyr::arrange(dplyr::desc(pct.missing)) %>% 
    dplyr::mutate(STRATA_SELECT = as.character(STRATA_SELECT))
  
  # Keep string to reorder factor based on marginal % missingness
  level.miss <- pct.miss.col$STRATA_SELECT
  pct.miss.col <- dplyr::mutate(
    .data = pct.miss.col,
    STRATA_SELECT = factor(STRATA_SELECT, levels = level.miss, ordered = TRUE))

  # summarize % missing for each total number missing in each factor level (label)
  miss.smry <- miss.smry %>%
    dplyr::group_by(STRATA_SELECT, num.missing.total) %>%
    dplyr::summarize(
      n = length(pct.missing),
      mean.miss = mean(pct.missing),
      lci = stats::quantile(pct.missing, (1 - ci) / 2),
      uci = stats::quantile(pct.missing, (1 + ci) / 2)
    ) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(
      STRATA_SELECT = factor(STRATA_SELECT, levels = level.miss, ordered = TRUE))
  
  # Tidy result of lm
  lm.res <- miss.smry %>% 
    split(x = ., f = .$STRATA_SELECT) %>% 
    purrr::map_df(~ broom::tidy(stats::lm(
      mean.miss ~ log10(num.missing.total), weights = n, data = .)), .id = "STRATA_SELECT")
  
  lm.res.name <- stringi::stri_join(
    "pct.missing.lm.", "strata.", strata.select)
  res[[lm.res.name]] <- lm.res
  
  
  miss.lm.coefs <- lm.res %>%
    dplyr::select(-c(std.error, statistic, p.value)) %>% 
    tidyr::spread(data = ., key = term, value = estimate) %>% 
    dplyr::rename(a = `(Intercept)`, b = `log10(num.missing.total)`) %>% 
    dplyr::mutate(
      STRATA_SELECT = factor(STRATA_SELECT, levels = level.miss, ordered = TRUE)) %>% 
    dplyr::arrange(STRATA_SELECT)
  
  axis.title.element <- ggplot2::element_text(
    size = 12, family = "Helvetica", face = "bold")
  axis.text.element <- ggplot2::element_text(
    size = 10, family = "Helvetica")
  fig <- ggplot2::ggplot(
    miss.smry, ggplot2::aes_string(x = "num.missing.total")) +
    ggplot2::geom_segment(
      ggplot2::aes_string(xend = "num.missing.total", y = "lci", yend = "uci"),
      color = "gray50") +
    ggplot2::geom_point(ggplot2::aes_string(y = "mean.miss", size = "n")) +
    ggplot2::scale_size(breaks = c(0, 10, 100, 1000, 10000)) +
    ggplot2::geom_abline(data = miss.lm.coefs,
                         ggplot2::aes(intercept = a, slope = b),
                         color = "yellow") +
    ggplot2::geom_hline(ggplot2::aes_string(yintercept = "ran.pct.miss")) +
    ggplot2::geom_hline(data = pct.miss.col,
                        ggplot2::aes_string(yintercept = "pct.missing"),
                        color = "red") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = stringi::stri_join("Strata: ", strata.select),
      x = expression(paste("Total number of missing genotypes (", log[10], ")")),
      y = "Missing genotypes (mean percentage)",
      size = "Markers missing (number)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = axis.title.element,
      axis.text.x = axis.text.element,
      axis.title.y = axis.title.element,
      axis.text.y = axis.text.element,
      legend.title = axis.title.element,
      legend.text = axis.text.element
    ) +
    ggplot2::facet_wrap(~ STRATA_SELECT)
  # fig
  
  fig.name <- stringi::stri_join(
    "pct.missing.plot", ".strata.", strata.select)
  
  if (write.plot) {
    ggplot2::ggsave(
      filename = stringi::stri_join(path.folder, "/", fig.name, ".pdf"),
      plot = fig,
      width = 20, height = 15,
      dpi = 600, units = "cm",
      useDingbats = FALSE)
  }
  
  res[[fig.name]] <- fig
  
  return(res)
}#End pct_missing_by_total
