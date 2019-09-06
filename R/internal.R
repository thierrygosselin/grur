#' @title ind_genotyped_helper
#' @description Help individual's genotyped threshold
#' @rdname ind_genotyped_helper
#' @export
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
    ggplot2::scale_x_continuous(name = "Individual's missing genotyped threshold (percent)", breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
    ggplot2::scale_y_continuous(name = "Individuals\n(blacklisted number)", breaks = y.breaks, limits = c(0, y.breaks.max)) +
    ggplot2::theme(
      axis.title.x = axis.title.element.text.fig,
      axis.title.y = axis.title.element.text.fig,
      axis.text.x = axis.text.element.text.fig,
      axis.text.y = axis.text.element.text.fig
    ) +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(~POP_ID)
  # plot.ind.geno.threshold
  return(plot.ind.geno.threshold)
}#End ind_genotyped_helper

#' @title blacklists_id_generator
#' @description Generate blacklist of ids
#' @rdname blacklists_id_generator
#' @export
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
#' @export
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
  
  n.whitelisted.markers <- nrow(whitelist.missing.geno)
  n.markers <- nrow(y)
  if (n.whitelisted.markers > 0 && n.whitelisted.markers < n.markers) {
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
#' @export
#' @keywords internal
markers_genotyped_helper <- function(x, y) {
  # x <- res$missing.genotypes.markers.pop
  # Set the breaks for the figure
  max.markers <- dplyr::n_distinct(x$MARKERS)
  
  threshold.helper.overall <- y %>% 
    dplyr::ungroup(.) %>% 
    dplyr::summarise(
      `0` = length(PERCENT[PERCENT == 0]),
      `10` = length(PERCENT[PERCENT <= 10]),
      `20` = length(PERCENT[PERCENT <= 20]),
      `30` = length(PERCENT[PERCENT <= 30]),
      `40` = length(PERCENT[PERCENT <= 40]),
      `50` = length(PERCENT[PERCENT <= 50]),
      `60` = length(PERCENT[PERCENT <= 60]),
      `70` = length(PERCENT[PERCENT <= 70]),
      `80` = length(PERCENT[PERCENT <= 80]),
      `90` = length(PERCENT[PERCENT <= 90]),
      `100` = length(PERCENT[PERCENT <= 100])
    ) %>%
    tidyr::gather(data = ., key = GENOTYPED_THRESHOLD, value = NUMBER_MARKERS) %>%
    dplyr::mutate(POP_ID = rep("OVERALL", n()))
  
  threshold.helper.pop <- x %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::summarise(
      `0` = length(PERCENT[PERCENT == 0]),
      `10` = length(PERCENT[PERCENT <= 10]),
      `20` = length(PERCENT[PERCENT <= 20]),
      `30` = length(PERCENT[PERCENT <= 30]),
      `40` = length(PERCENT[PERCENT <= 40]),
      `50` = length(PERCENT[PERCENT <= 50]),
      `60` = length(PERCENT[PERCENT <= 60]),
      `70` = length(PERCENT[PERCENT <= 70]),
      `80` = length(PERCENT[PERCENT <= 80]),
      `90` = length(PERCENT[PERCENT <= 90]),
      `100` = length(PERCENT[PERCENT <= 100])
    ) %>%
    tidyr::gather(data = ., key = GENOTYPED_THRESHOLD, value = NUMBER_MARKERS, -POP_ID)
  
  mean.pop <- threshold.helper.pop %>%
    dplyr::group_by(GENOTYPED_THRESHOLD) %>%
    dplyr::summarise(NUMBER_MARKERS = round(mean(NUMBER_MARKERS), 0)) %>%
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
    ggplot2::scale_x_continuous(name = "Marker's missing genotyped threshold (percent)", breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
    ggplot2::scale_y_continuous(name = "Markers\n(whitelisted number)", breaks = y.breaks, limits = c(0, y.breaks.max)) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica"),#, angle = 90, hjust = 1, vjust = 0.5),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    ) +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(~POP_ID)
  # plot.markers.geno.threshold
  return(plot.markers.geno.threshold)
}#End markers_genotyped_helper


#' @title generate_pcoa_plot
#' @description Generate the PCoA plots
#' @rdname generate_pcoa_plot
#' @keywords internal
#' @export
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
    
    ibm.plot <- ggplot2::ggplot(
      vectors,
      ggplot2::aes_string(
        x = stringi::stri_join("Axis.", pcx),
        y = stringi::stri_join("Axis.", pcy), size = vectors$MISSING_GENOTYPE_PERCENT),
      environment = environment()) +
      ggplot2::geom_point(ggplot2::aes_string(colour = strata.select), alpha = 0.5) +
      ggplot2::labs(x = stringi::stri_join("PCo", pcx, " [", variance.component[pcx,2], "]")) +
      ggplot2::labs(y = stringi::stri_join("PCo", pcy, " [", variance.component[pcy,2], "]")) +
      ggplot2::scale_size_area(name = "Individual's\nmissing genotypes\n(percent)", max_size = 4) +
      ggplot2::theme(
        axis.title.x = element.text.fig,
        axis.title.y = element.text.fig,
        legend.title = element.text.fig,
        legend.text = element.text.fig
      ) +
      ggplot2::theme_bw()
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
  
# message("ggarrange problem")
# print(pcoa.plots.strata[[1]])
# message("worked")
  pcoa.plots.strata <- ggpubr::ggarrange(
    pcoa.plots.strata[[1]],
    pcoa.plots.strata[[2]],
    pcoa.plots.strata[[3]],
    pcoa.plots.strata[[4]],
    pcoa.plots.strata[[5]],
    pcoa.plots.strata[[6]],
    ncol = 2, nrow = 3, legend = "right", common.legend = TRUE
  )
  
  # message("write plot problem")
  
  if (write.plot) {
    plot.name <- stringi::stri_join("ibm.plots.strata.", strata.select, ".pdf")
    ggplot2::ggsave(
      filename = file.path(path.folder, plot.name),
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
#' @keywords internal
#' @export

pct_missing_by_total <- function(
  strata.select, data, ci = 0.95, path.folder, write.plot = TRUE) {
  
  # rename strata column and convert GT_BIN to missing or not
  data %<>%  
    # dplyr::rename(STRATA_SELECT = data[[!!(strata.select)]]) %>% 
    dplyr::rename(STRATA_SELECT = !! strata.select) %>% 
    dplyr::mutate(is.missing = GT_MISSING_BINARY == 0)
  
  # count number of strata
  n.strata <- dplyr::n_distinct(data$STRATA_SELECT)
  
  # summarize missingness by locus and factor column
  miss.smry <- data %>% 
    dplyr::group_by(MARKERS, STRATA_SELECT) %>% 
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
    dplyr::summarise(
      pct.missing = sum(num.missing.col) / sum(num.missing.total)
    ) %>%
    dplyr::ungroup(.) %>% 
    dplyr::arrange(dplyr::desc(pct.missing)) %>% 
    dplyr::mutate(STRATA_SELECT = as.character(STRATA_SELECT))
  
  # Keep string to reorder factor based on marginal % missingness
  level.miss <- pct.miss.col$STRATA_SELECT
  pct.miss.col <- pct.miss.col %>% 
    dplyr::mutate(
      STRATA_SELECT = factor(STRATA_SELECT, levels = level.miss, ordered = TRUE)
    )

  # summarize % missing for each total number missing in each factor level (label)
  miss.smry %<>% 
    dplyr::group_by(STRATA_SELECT, num.missing.total) %>%
    dplyr::summarize(
      n = length(pct.missing),
      mean.miss = mean(pct.missing),
      lci = stats::quantile(pct.missing, (1 - ci) / 2),
      uci = stats::quantile(pct.missing, (1 + ci) / 2)
    ) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(
      STRATA_SELECT = factor(STRATA_SELECT, levels = level.miss, ordered = TRUE)
    )
  
  # to store results
  res <- list() 
  
  # Tidy result of lm
  lm.res.name <- stringi::stri_join("pct.missing.lm.", "strata.", strata.select)
  res[[lm.res.name]] <- miss.smry %>% 
    split(x = ., f = .$STRATA_SELECT) %>% 
    purrr::map_df(
      ~ broom::tidy(
        stats::lm(
          mean.miss ~ log10(num.missing.total), 
          weights = n, 
          data = .
        )
      ), 
      .id = "STRATA_SELECT"
    )
  
  miss.lm.coefs <- res[[lm.res.name]] %>%
    dplyr::select(-c(std.error, statistic, p.value)) %>% 
    tidyr::spread(data = ., key = term, value = estimate) %>% 
    dplyr::rename(a = `(Intercept)`, b = `log10(num.missing.total)`) %>% 
    dplyr::mutate(
      STRATA_SELECT = factor(STRATA_SELECT, levels = level.miss, ordered = TRUE)
    ) %>% 
    dplyr::arrange(STRATA_SELECT)
  
  axis.title.element <- ggplot2::element_text(
    size = 12, family = "Helvetica", face = "bold"
  )
  axis.text.element <- ggplot2::element_text(size = 10, family = "Helvetica")
  
  # generate figure
  fig.name <- stringi::stri_join("pct.missing.plot", ".strata.", strata.select)
  res[[fig.name]] <- ggplot2::ggplot(
    miss.smry, ggplot2::aes_string(x = "num.missing.total")) +
    ggplot2::geom_segment(
      ggplot2::aes_string(xend = "num.missing.total", y = "lci", yend = "uci"),
      color = "gray50"
    ) +
    ggplot2::geom_point(ggplot2::aes_string(y = "mean.miss", size = "n")) +
    ggplot2::scale_size(breaks = c(0, 10, 100, 1000, 10000)) +
    ggplot2::geom_abline(
      data = miss.lm.coefs,
      ggplot2::aes(intercept = a, slope = b),
      color = "yellow"
    ) +
    # the expected % missing by random (1 / number of factor levels)
    ggplot2::geom_hline(yintercept = 1 / dplyr::n_distinct(data$STRATA_SELECT)) +
    ggplot2::geom_hline(
      data = pct.miss.col,
      ggplot2::aes_string(yintercept = "pct.missing"),
      color = "red"
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = stringi::stri_join("Strata: ", strata.select),
      x = expression(
        paste("Total number of missing genotypes (", log[10], ")")
      ),
      y = "Missing genotypes (mean percentage)",
      size = "Markers missing (number)"
    ) +
    ggplot2::theme(
      axis.title.x = axis.title.element,
      axis.text.x = axis.text.element,
      axis.title.y = axis.title.element,
      axis.text.y = axis.text.element,
      legend.title = axis.title.element,
      legend.text = axis.text.element
    ) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~ STRATA_SELECT)
  # fig
  
  if (write.plot) {
    ggplot2::ggsave(
      filename = stringi::stri_join(path.folder, "/", fig.name, ".pdf"),
      plot = res[[fig.name]],
      width = n.strata * 3, 
      height = n.strata * 3,
      dpi = 600, 
      units = "cm",
      useDingbats = FALSE, 
      limitsize = FALSE
    )
  }
  
  return(res)
}#End pct_missing_by_total



#' @title missing_rda
#' @description The function uses Permutation test Anova and
#' Redundancy Analysis to highlight potential patterns of missingness in the data
#' based on the strata provided
#' @rdname missing_rda
#' @keywords internal
#' @export
# IMPORT FROM ------------------------------------------------------------------
#' @importFrom dplyr distinct rename arrange mutate select summarise group_by ungroup filter inner_join left_join
#' @importFrom ggpubr ggarrange
#' @importFrom ape pcoa
#' @importFrom radiator tidy_genomic_data change_pop_names ibdg_fh detect_all_missing write_rad
#' @importFrom cowplot plot_grid align_plots
#' @importFrom Matrix Matrix
#' @importFrom vegan anova.cca rda
#' @importFrom ape pcoa
#' @importFrom adespatial dist.ldc
#' @importFrom pbmcapply pbmclapply
#' @importFrom data.table as.data.table dcast.data.table melt.data.table

missing_rda <- function(
  data, 
  strata, 
  permutations = 1000, 
  parallel.core = parallel::detectCores() - 1
) {
  res <- list()
  if (is.vector(data)) {
    res$data.pcoa <- radiator::read_rad(data)
  } else {
    res$data.pcoa <- data
  }
  data <- NULL
  res$data.pcoa %<>% dplyr::ungroup(.)
  
  # PCoA on Jaccard distance (binary)
  res$data.pcoa <- ape::pcoa(
    D = sqrt(adespatial::dist.ldc(
      Y = res$data.pcoa,
      method = "jaccard",
      binary = TRUE, silent = TRUE)), 
    rn = strata$INDIVIDUALS)
  # D is Euclidean because the function outputs D[jk] = sqrt(1-S[jk])
  
  # Anova on the pcoa data
  rda_anova <- function(
    strata.select, 
    data.pcoa, 
    strata, 
    permutations = 100,
    parallel.core = parallel::detectCores() - 1) {
    
    # strata.select <- "POP_ID" #test
    # data.pcoa <- res$data.pcoa #test
    
    res.rda <- list()#to store results
    
    # Generate the function formula --------------------------------------------
    # based on increasing levels for available strata
    grur_count <- function(x) length(unique(x))
    strata.formula <- dplyr::select(strata, dplyr::one_of(strata.select)) %>%
      dplyr::summarise_all(.tbl = ., .funs = grur_count) %>% 
      tidyr::gather(data = ., key = TERMS, value = LEVELS) %>%
      dplyr::arrange(LEVELS) %>% 
      dplyr::select(TERMS) %>%
      purrr::flatten_chr(.)
    
    formula.grur <- stats::as.formula(
      paste("data.pcoa$vectors ~ ", paste(strata.formula, collapse= "+")))
    
    # Check how many strata.select are used and messages -----------------------
    # if (length(strata.select) > 1) {
    message("Redundancy Analysis using strata: ",
            stringi::stri_join(strata.select, collapse = ", "))
    rda_strata_name <- stringi::stri_join("rda.strata", strata.select, collapse = "_", sep = "_")
    anova_strata_name <- stringi::stri_join("anova.strata.", strata.select, collapse = "_", sep = "_")
    # } else {
    #   rda_strata_name <- stringi::stri_join("rda.strata.", strata.select)
    #   anova_strata_name <- stringi::stri_join("anova.strata.", strata.select)
    #   strata.select <- rlang::sym(strata.select)
    #   message("Redundancy Analysis using strata: ", strata.select)
    # }
    message("RDA model formula: ", format(formula.grur))
    
    #RDA -----------------------------------------------------------------------
    
    data.rda <- vegan::rda(formula.grur, strata)
    # data.rda <- vegan::rda(rlang::`f_rhs<-`(data.pcoa$vectors ~ ., strata.select), strata)
    # data.rda <- vegan::rda(rlang::`f_rhs<-`(data.pcoa$vectors ~ ., stats::reformulate(termlabels = strata.select)), strata)
    # data.rda <- vegan::rda(stats::reformulate(termlabels = strata.select, response = data.pcoa$vectors), data = strata)
    # data.rda <- vegan::rda(res$data.pcoa$vectors ~ POP_ID + POP_TYPE + ECOTYPE, strata)
    
    # ANOVA overall test
    message("Permutation test for Redundancy Analysis using strata: ", stringi::stri_join(strata.select, collapse = ", "))
    # data.anova <- vegan::anova.cca(object = data.rda, permutations = permutations, parallel = parallel.core)
    # data.anova <- suppressWarnings(broom::tidy(vegan::anova.cca(object = data.rda, permutations = permutations, parallel = parallel.core)))
    
    data.anova <- suppressWarnings(broom::tidy(vegan::anova.cca(object = data.rda, by = "terms", model = "direct", permutations = permutations, parallel = parallel.core)))
    p.value.message <- dplyr::select(data.anova, STRATA = term, VARIANCE = Variance, P_VALUE = p.value) %>%
      dplyr::filter(STRATA %in% strata.formula)
    
    message("\nHypothesis based on the strata provided")
    message("    Null Hypothesis (H0): No pattern of missingness in the data between strata")
    message("    Alternative Hypothesis (H1): Presence of pattern(s) of missingness in the data between strata\n")
    # message("    p-value: ", data.anova$p.value[1], "\n")
    print(p.value.message)
    message("    note: low p-value -> reject the null hypothesis\n")
    #return results
    res.rda[[rda_strata_name]] <- data.rda
    res.rda[[anova_strata_name]] <- data.anova
    return(res.rda)
  }
  
  strata.select <- purrr::keep(.x = colnames(strata),
                               .p = !colnames(strata) %in% "INDIVIDUALS")
  
  if (length(strata.select) > 1) {
    message("\nSeparate analysis of strata\n")
    res$rda_separate_strata <- purrr::map(
      .x = strata.select, 
      .f = rda_anova,
      data.pcoa = res$data.pcoa,
      strata = strata,
      permutations = permutations,
      parallel.core = parallel.core) %>%
      purrr::flatten(.)
    
    message("\nCombined strata analysis\n")
    
    res$rda_combined_strata <- rda_anova(
      strata.select = strata.select,
      data.pcoa = res$data.pcoa,
      strata = strata,
      permutations = permutations,
      parallel.core = parallel.core)
  } else {
    res$rda_combined_strata <- rda_anova(
      strata.select = strata.select,
      data.pcoa = res$data.pcoa,
      strata = strata,
      permutations = permutations,
      parallel.core = parallel.core)
  }
  
  return(res)
}#End missing_rda


#' @title make_mig_mat
#' @description Make migration matrix for simulation
#' @rdname make_mig_mat
#' @keywords internal
#' @export
make_mig_mat <- function(mig.rate, num.pops, 
                         type = c("island", "stepping.stone")) {
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


#' @title run_fsc_sim
#' @description Run fastsimcoal
#' @rdname run_fsc_sim
#' @keywords internal
#' @export
run_fsc_sim <- function(sc, num.rep, num.cores, exec) {
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
  
  demes <- do.call(
    strataG::fscSettingsDemes,
    c(
      lapply(1:sc$num.pops, function(i) {
        strataG::fscDeme(deme.size = sc$ne, sample.size = sc$ne)
      }),
      ploidy = 2
    )
  )
  
  strataG::fscWrite(
    demes = demes,
    migration = if(sc$num.pops > 1) {
      strataG::fscSettingsMigration(sc$mig.mat)
    } else NULL,
    events = .makeEventSettings(sc$div.time, sc$num.pops),
    genetics = strataG::fscSettingsGenetics(
      strataG::fscBlock_snp(1, 1e-4), num.chrom = 1000
    ),
    label = "grur.fsc.sim"
  ) %>% 
    strataG::fscRun(
      num.sims = num.rep, num.cores = num.cores, quiet = TRUE, exec = exec
    )
}


#' @title calc_freqs
#' @description Calculate by-population and by-locus allele frequencies
#' @rdname calc_freqs
#' @keywords internal
#' @export
calc_freqs <- function(snps) {
  .alleleProp <- function(mac) {
    maf <- mean(mac == 0) + (mean(mac == 1) / 2)
    c('1' = maf, '2' = 1 - maf)
  }
  
  mac.df <- snps %>% 
    strataG::df2gtypes(ploidy = 2) %>% 
    strataG::as.data.frame(coded = T) %>% 
    dplyr::select(-id) %>% 
    tidyr::gather("locus", "mac", -stratum) %>% 
    dplyr::mutate(
      stratum = as.numeric(factor(stratum)),
      locus = as.numeric(factor(locus))
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


#' @title run_rmetasim
#' @description Run Rmetasim
#' @rdname run_rmetasim
#' @keywords internal
#' @export
run_rmetasim <- function(freqs, sc, num.gens) {
  cat(format(Sys.time()), "running rmetasim\n")
  
  localS <- matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2)
  localR <- matrix(c(0, 0, 1.2, 0), nrow = 2, ncol = 2)
  localM <- matrix(c(0, 0, 0, 1.2), nrow = 2, ncol = 2)
  S <- M <- matrix(0, nrow = sc$num.pops * 2, ncol = sc$num.pops * 2)
  diag(S) <- diag(M) <- 1
  R <- rmetasim::landscape.mig.matrix(
    h = nrow(sc$mig.mat), 
    s = 2, 
    mig.model = "custom", 
    R.custom = sc$mig.mat
  )$R
  
  Rland <- rmetasim::landscape.new.empty() %>% 
    rmetasim::landscape.new.intparam(
      h = sc$num.pops, s = 2, cg = 0, ce = 0, totgen = num.gens + 1
    ) %>% 
    rmetasim::landscape.new.switchparam() %>% 
    rmetasim::landscape.new.floatparam() %>% 
    rmetasim::landscape.new.local.demo(localS, localR, localM) %>% 
    rmetasim::landscape.new.epoch(R = R, carry = rep(sc$ne, sc$num.pops)) 
  
  for(i in 1:length(freqs$global)) {
    Rland <- rmetasim::landscape.new.locus(
      Rland, type = 2, ploidy = 2, mutationrate = 0,
      transmission = 0, numalleles = 2, allelesize = 1,
      frequencies = freqs$global[[i]], states = names(freqs$global[[i]])
    )
  }
  
  Rland <- Rland %>% 
    rmetasim::landscape.new.individuals(rep(c(sc$ne, 0), sc$num.pops)) %>% 
    rmetasim::landscape.setpopfreq(freqs$pop) 
  
  for(gen.i in 1:num.gens) {
    Rland <- rmetasim::landscape.simulate(Rland, numit = 1)
    to.kill <- tapply(
      1:nrow(Rland$individuals), 
      Rland$individuals[, 1], 
      function(i) if (length(i) > sc$ne) sample(i, length(i) - sc$ne) else NULL
    )
    to.kill <- unlist(unname(to.kill))
    if (length(to.kill) > 0) Rland$individuals <- Rland$individuals[-to.kill, ]
  }
  
  Rland
}
