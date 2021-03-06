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
    tidyr::pivot_longer(
      data = .,
      cols = tidyselect::everything(),
      names_to = "GENOTYPED_THRESHOLD",
      values_to = "NUMBER_INDIVIDUALS"
    ) %>%
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
    tidyr::pivot_longer(
      data = .,
      cols = -POP_ID,
      names_to = "GENOTYPED_THRESHOLD",
      values_to = "NUMBER_INDIVIDUALS"
    )
  
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
      file.path(path.folder,  paste0(as.name(blacklist.name), ".tsv")))
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
    tidyr::pivot_longer(
      data = .,
      cols = tidyselect::everything(),
      names_to = "GENOTYPED_THRESHOLD",
      values_to = "NUMBER_MARKERS"
    ) %>% 
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
    tidyr::pivot_longer(
      data = .,
      cols = -POP_ID,
      names_to = "GENOTYPED_THRESHOLD",
      values_to = "NUMBER_MARKERS"
    )
  
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
    ncol = 2, 
    nrow = 3, 
    legend = "right", 
    common.legend = TRUE
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
    grur::rad_wide(x = ., names_from = term, values_from = estimate, tidy = TRUE) %>% 
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
      tidyr::pivot_longer(
        data = .,
        cols = tidyselect::everything(),
        names_to = "TERMS",
        values_to = "LEVELS"
      ) %>%
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



# parallel_core_opt ------------------------------------------------------------
#' @title parallel_core_opt
#' @description Optimization of parallel core argument for radiator
#' @keywords internal
#' @export
parallel_core_opt <- function(parallel.core = NULL, max.core = NULL) {
  # strategy:
  # minimum of 1 core and a maximum of all the core available -2
  # even number of core
  # test
  # parallel.core <- 1
  # parallel.core <- 2
  # parallel.core <- 3
  # parallel.core <- 11
  # parallel.core <- 12
  # parallel.core <- 16
  # max.core <- 5
  # max.core <- 50
  # max.core <- NULL
  
  # Add-ons options
  # to control the max and min number to use...
  
  if (is.null(parallel.core)) {
    parallel.core <- parallel::detectCores() - 2
  } else {
    parallel.core <- floor(parallel.core / 2) * 2
    parallel.core <- max(1, min(parallel.core, parallel::detectCores() - 2))
  }
  
  if (is.null(max.core)) {
    parallel.core.opt <- parallel.core
  } else {
    parallel.core.opt <- min(parallel.core, floor(max.core / 2) * 2)
  }
  return(parallel.core.opt)
}#End parallel_core_opt

# using future and future.apply -------------------------------------------------
#' @name grur_future
#' @title grur parallel function
#' @description Updating grur to use future
# @inheritParams future::plan
# @inheritParams future::availableCores
#' @inheritParams future.apply::future_apply
#' @rdname grur_future
#' @export
#' @keywords internal
grur_future <- function(
  .x,
  .f,
  flat.future = c("int", "chr", "dfr", "dfc", "walk", "drop"),
  split.vec = FALSE,
  split.with = NULL,
  split.chunks = 4L,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  opt.change <- getOption("width")
  options(width = 70)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(if (parallel.core > 1L) future::plan(strategy = "sequential"), add = TRUE)
  
  # argument for flattening the results
  flat.future <- match.arg(
    arg = flat.future,
    choices = c("int", "chr", "dfr", "dfc", "walk", "drop"),
    several.ok = FALSE
  )
  
  # splitting into chunks-------------------------------------------------------
  if (split.vec && is.null(split.with)) {
    # d: data, data length, data size
    # sv: split vector
    d <- .x
    df <- FALSE
    if (any(class(d) %in% c("tbl_df","tbl","data.frame"))) {
      d <- nrow(d)
      df <- TRUE
    }
    if (length(d) > 1L) d <- length(d)
    stopifnot(is.integer(d))
    sv <- as.integer(floor((split.chunks * (seq_len(d) - 1) / d) + 1))
    # sv <- as.integer(floor((parallel.core * cpu.rounds * (seq_len(d) - 1) / d) + 1))
    stopifnot(length(sv) == d)
    
    # split
    if (df) {
      .x$SPLIT_VEC <- sv
      .x %<>% dplyr::group_split(.tbl = ., "SPLIT_VEC", .keep = FALSE)
    } else {
      .x %<>% split(x = ., f = sv)
    }
  }
  if (!is.null(split.with)) {
    # check
    if (length(split.with) != 1 || !is.character(split.with)) {
      rlang::abort(message = "Contact author: problem with parallel computation")
    }
    .data <- NULL
    stopifnot(rlang::has_name(.x, split.with))
    if (split.vec) {
      sv <- dplyr::distinct(.x, .data[[split.with]])
      d <- nrow(sv)
      sv$SPLIT_VEC <- as.integer(floor((split.chunks * (seq_len(d) - 1) / d) + 1))
      .x %<>%
        dplyr::left_join(sv, by = split.with) %>%
        dplyr::group_split(.tbl = ., "SPLIT_VEC", .keep = FALSE)
    } else {
      .x %<>% dplyr::group_split(.tbl = ., .data[[split.with]], .keep = TRUE)
    }
  }
  
  
  
  if (parallel.core == 1L) {
    future::plan(strategy = "sequential")
  } else {
    parallel.core <- parallel_core_opt(parallel.core = parallel.core)
    lx <- length(.x)
    if (lx < parallel.core) {
      future::plan(strategy = "multisession", workers = lx)
    } else {
      future::plan(strategy = "multisession", workers = parallel.core)
    }
  }
  
  # .x <- future.apply::future_apply(X = .x, FUN = .f, ...)
  # capture dots
  # d <- rlang::dots_list(..., .ignore_empty = "all", .preserve_empty = TRUE, .homonyms = "first")
  # if (bind.rows) .x %<>% dplyr::bind_rows(.)
  
  
  
  # Run the function in parallel and account for dots-dots-dots argument
  rad_map <- switch(flat.future,
                    int = {furrr::future_map_int},
                    chr = {furrr::future_map_chr},
                    dfr = {furrr::future_map_dfr},
                    dfc = {furrr::future_map_dfc},
                    walk = {furrr::future_walk},
                    drop = {furrr::future_map}
  )
  
  opts <- furrr::furrr_options(globals = FALSE)
  if (length(list(...)) == 0) {
    .x %<>% rad_map(.x = ., .f = .f, .options = opts)
  } else {
    .x %<>% rad_map(.x = ., .f = .f, ..., .options = opts)
  }
  return(.x)
}#End grur_future


# PIVOT-GATHER-CAST ------------------------------------------------------------
# rationale for doing this is that i'm tired of using tidyverse or data.table semantics
# tidyr changed from gather/spread to pivot_ functions but their are still very slow compared
# to 1. the original gather/spread and data.table equivalent...

#' @title rad_long
#' @description Gather, melt and pivot_longer
#' @rdname rad_long
#' @keywords internal
#' @export

rad_long <- function(
  x,
  cols = NULL,
  measure_vars = NULL,
  names_to = NULL,
  values_to = NULL,
  variable_factor = TRUE,
  keep_rownames = FALSE,
  tidy = FALSE
){
  
  
  # tidyr
  if (tidy) {
    x %>%
      tidyr::pivot_longer(
        data = .,
        cols = -cols,
        names_to = names_to,
        values_to = values_to
      )
  } else {# data.table
    x %>%
      data.table::as.data.table(., keep.rownames = keep_rownames) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = cols,
        measure.vars = measure_vars,
        variable.name = names_to,
        value.name = values_to,
        variable.factor = variable_factor
      ) %>%
      tibble::as_tibble(.)
  }
}#rad_long

#' @title rad_wide
#' @description Spread, dcast and pivot_wider
#' @rdname rad_wide
#' @keywords internal
#' @export
rad_wide <- function(
  x ,
  formula = NULL,
  names_from = NULL,
  values_from = NULL,
  values_fill = NULL,
  sep = "_",
  tidy = FALSE
  
){
  # tidyr
  if (tidy) {
    x %<>%
      tidyr::pivot_wider(
        data = .,
        names_from = names_from,
        values_from = values_from,
        values_fill = values_fill
      )
  } else {# data.table
    x  %>%
      data.table::as.data.table(.) %>%
      data.table::dcast.data.table(
        data = .,
        formula =  formula,
        value.var = values_from,
        sep = sep,
        fill = values_fill
      ) %>%
      tibble::as_tibble(.)
  }
}#rad_wide
