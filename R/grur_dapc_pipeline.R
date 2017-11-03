#' @title grur DAPC and likelihood-based genetic clustering pipeline based on adegenet
#' @description Analyse data with adegenet DAPC and snapclust.
#' @param x Genind object
#' @param strata (character, optional) Default: \code{strata = "POP_ID"}.
#' @param pop.levels (character string, optional) Default: \code{pop.levels = NULL}.
#' @param n.pca.max (integer, optional) Default: \code{n.pca.max = 300}.
#' @param n.rep (integer, optional) Default: \code{n.rep = 100}.
#' @param training.set (double, optional) Default: \code{training.set = 0.9}.
#' @param pc (integer string, optional) Default: \code{pc = 1:4}.
#' @param plot.filename (character, optional) Default: \code{plot.filename = NULL}.
#' @param parallel.core (integer, optional) The number of core used for parallel
#' execution during the pipeline.
#' Default: \code{parallel::detectCores() - 1}.
#' @return coming
#' @export
#' @rdname dapc_clustering
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' 
dapc_clustering <- function(
  x,
  strata = "POP_ID",
  pop.levels = NULL,
  n.pca.max = 300,
  n.rep = 100,
  training.set = 0.9,
  pc = 1:4,
  plot.filename = NULL,
  parallel.core = parallel::detectCores() - 1) {
  
  opt.change <- getOption("width")
  options(width = 70)
  cat("\n\n#######################################################################\n")
  cat("################## grur::dapc_clustering pipeline #####################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  if(!is.null(plot.filename)) message("Analyzing: ", plot.filename)
  
  if (missing(x)) stop("genind/genlight object required")
  if (is.null(pop.levels)) pop.levels <- levels(x@strata[,"POP_ID"])
  pop.num <- unique(stringi::stri_detect_regex(str = pop.levels, pattern = "^[0-9]+$"))
  if (length(pop.num) == 1 && pop.num) pop.levels <- as.character(sort(as.numeric(pop.levels)))
  
  # Color palette with grey:
  # dapc_color_palette_grey <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # Color palette with black:
  # dapc_color_palette_black <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  
  # Admixture figures ----------------------------------------------------------
  # if NA present, use alpha score
  if (anyNA(x@tab)) {
    message("DAPC calculated using alpha-score")
    dapc.best.optim.a.score <- adegenet::optim.a.score(adegenet::dapc(x, n.da = length(levels(x@strata[,strata])), n.pca = round((length(adegenet::indNames(x))/3) - 1, 0)), pop = x@strata[,strata], plot = FALSE)$best
    dapc.data <- adegenet::dapc(x, n.da = length(levels(x@strata[,strata])), n.pca = dapc.best.optim.a.score, pop = x@strata[,strata])
    
    dapc.assignment <- data.frame(
      POP_ID = x@strata[,"POP_ID"],
      INDIVIDUALS = adegenet::indNames(x),
      GROUP_PRIOR = dapc.data$grp,
      GROUP_POST = dapc.data$assign, dapc.data$posterior
    ) %>%
      tidyr::gather(key = ANCESTRY, value = "PERC", -c(POP_ID, INDIVIDUALS, GROUP_PRIOR, GROUP_POST)) %>%
      dplyr::mutate(POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE)) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS)
    
    n.pop <- dplyr::n_distinct(dapc.assignment$POP_ID)
    
    # figure
    dapc.plot <- ggplot2::ggplot(dapc.assignment, ggplot2::aes(x = INDIVIDUALS, y = PERC, fill = ANCESTRY)) +
      ggplot2::geom_bar(stat = "identity",  position = "fill") +
      ggplot2::scale_y_continuous(expand  =  c(0,  0)) + #in order for the line not to expand beyond the graph!
      ggplot2::labs(y = "Ancestry") +
      ggplot2::labs(x = "Individuals") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major.x  =  ggplot2::element_blank(),
        panel.grid.minor.x  =  ggplot2::element_blank(),
        panel.grid.major.y  =  ggplot2::element_line(colour = "grey60",  linetype = "dashed"),
        axis.title.x = ggplot2::element_text(size = 10,  family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = 10,  family = "Helvetica", face = "bold"),
        axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica"),
        axis.ticks.x = ggplot2::element_blank(),
        strip.text.y = ggplot2::element_text(angle = 0),
        panel.spacing.y = ggplot2::unit(0.5, "lines"),
        panel.spacing.x = ggplot2::unit(0.1, "lines"),
        legend.position = "right"
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(
        label.position = "right",
        title.position = "top",
        title.hjust = 0.5,
        label.hjust = 0.5)
      ) +
      ggplot2::facet_grid(~POP_ID, scales = "free", space = "free_x")
    
    res = list(dapc.best.optim.a.score = dapc.best.optim.a.score,
               dapc.data = dapc.data,
               dapc.assignment = dapc.assignment,
               dapc.plot = dapc.plot)
  } else {# No NA use XVAL
    message("DAPC calculated using cross-validation")
    xval <- adegenet::xvalDapc(
      x = x@tab,
      grp = factor(as.character(x@strata[,strata])),
      # grp = adegenet::pop(x),
      n.pca.max = n.pca.max,
      n.rep = n.rep ,
      training.set = training.set,
      result = "groupMean",
      center = TRUE,
      scale = FALSE,
      xval.plot = TRUE,
      parallel = "multicore",
      ncpus = parallel.core
    )
    
    dapc.data <- xval$DAPC
    
    # for admixture type figure:
    dapc.assignment <- data.frame(
      POP_ID = x@strata[,"POP_ID"],
      INDIVIDUALS = adegenet::indNames(x),
      GROUP_PRIOR = dapc.data$grp,
      GROUP_POST = dapc.data$assign, dapc.data$posterior
    ) %>%
      tidyr::gather(key = ANCESTRY, value = "PERC", -c(POP_ID, INDIVIDUALS, GROUP_PRIOR, GROUP_POST)) %>%
      dplyr::mutate(
        POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE),
        ANCESTRY = stringi::stri_replace_all_fixed(ANCESTRY, "X", "K", vectorize_all = FALSE),
        ANCESTRY = factor(ANCESTRY, ordered = TRUE)
      ) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS)
    
    n.pop <- dplyr::n_distinct(dapc.assignment$POP_ID)
    
    
    # figure
    dapc.plot <- ggplot2::ggplot(dapc.assignment, ggplot2::aes(x = INDIVIDUALS, y = PERC, fill = ANCESTRY)) +
      ggplot2::geom_bar(stat = "identity",  position = "fill") +
      ggplot2::scale_y_continuous(expand  =  c(0,  0)) + #in order for the line not to expand beyond the graph!
      # ggplot2::scale_fill_manual(values = dapc_color_palette_black) +
      # ggplot2::scale_fill_brewer(type = "div", palette = palette.dapc) +
      ggplot2::labs(
        title = "adegenet DAPC Cross-Validation",
        x = "Individuals",
        y = "Ancestry") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 12, family = "Helvetica", face = "bold" ),
        panel.grid.major.x  =  ggplot2::element_blank(),
        panel.grid.minor.x  =  ggplot2::element_blank(),
        panel.grid.major.y  =  ggplot2::element_line(colour = "grey60",  linetype = "dashed"),
        axis.title.x = ggplot2::element_text(size = 10,  family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = 10,  family = "Helvetica", face = "bold"),
        axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica"),
        axis.ticks.x = ggplot2::element_blank(),
        strip.text.y = ggplot2::element_text(angle = 0),
        panel.spacing.y = ggplot2::unit(0.5, "lines"),
        panel.spacing.x = ggplot2::unit(0.1, "lines"),
        legend.position = "right"
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(
        label.position = "right",
        title.position = "top",
        title.hjust = 0.5,
        label.hjust = 0.5)
      ) +
      ggplot2::facet_grid(~POP_ID, scales = "free", space = "free_x")
    # dapc.plot
    res = list(xval = xval, dapc.assignment = dapc.assignment, dapc.plot = dapc.plot)
  }#xval
  
  if (!is.null(plot.filename)) {
    write_dapc_plot(x = res, plot.name = plot.filename, width = n.pop * 3, height = (n.pop * 3) / 2)
  }
  
  # PC figures ----------------------------------------------------------------- 
  # pc <- 1:4
  pc.to.do <- utils::combn(pc, 2, simplify = FALSE)
  
  # for PC type figure
  pop.coord <- tibble::as_data_frame(dapc.data$grp.coord) %>%
    dplyr::mutate(
      INDIVIDUALS = as.character(row.names(tibble::as_data_frame(dapc.data$grp.coord))),
      POP_ID = INDIVIDUALS,
      GROUP = INDIVIDUALS,
      ASSIGN = INDIVIDUALS
    )
  
  pc.data <- tibble::data_frame(
    INDIVIDUALS = as.character(row.names(tibble::as_data_frame(dapc.data$ind.coord))),
    POP_ID = dapc.data$grp,
    GROUP = dapc.data$grp,
    ASSIGN = dapc.data$assign) %>% 
    dplyr::bind_cols(tibble::as_data_frame(dapc.data$ind.coord)) %>% 
    dplyr::mutate(
      POP_ID = factor(x = POP_ID, levels = pop.levels),
      GROUP = factor(x = GROUP, levels = pop.levels),
      ASSIGN = factor(x = ASSIGN, levels = pop.levels)
    )
  
  # variance
  variance.component <-  tibble::data_frame(EIGENVALUES = dapc.data$eig) %>%
    dplyr::mutate(VARIANCE_PROP = round(EIGENVALUES/sum(EIGENVALUES), 2))
  
  res$pc.plots <- generate_dapc_pc_plot(pc.to.do = pc.to.do, pc.data = pc.data, variance.component = variance.component, plot.filename = plot.filename)
  
  
  # Densities of discriminant function -----------------------------------------
  density.prep <- pc.data %>% 
    tidyr::gather(data = ., key = "PC", value = "VALUE", -c(INDIVIDUALS, POP_ID, ASSIGN, GROUP))
  
  element.text.fig <- ggplot2::element_text(
    size = 12, family = "Helvetica", face = "bold")
  facet_names <- ggplot2::as_labeller(c(`LD1` = "PC1", `LD2` = "PC2", `LD3` = "PC3", `LD4` = "PC4"))
  
  res$ridges.plots <- ggplot2::ggplot(density.prep, ggplot2::aes(x = VALUE, y = POP_ID, na.rm = TRUE)) +
    ggridges::geom_density_ridges(ggplot2::aes(fill = POP_ID), alpha = 0.4, rel_min_height = 0.01, scale = 3) +
    # ggplot2::geom_density(ggplot2::aes(fill = POP_ID), alpha = 0.4) +
    ggplot2::labs(x = "Discriminant functions") +
    # ggplot2::labs(y = "Density") +
    # ggplot2::expand_limits(y = 0) +
    # ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(colour = NULL))) +
    # ggplot2::theme_minimal(base_size = 14) + 
    ggplot2::scale_x_continuous(expand = c(0.01, 0)) +
    ggplot2::scale_y_discrete(expand = c(0.01, 0)) +
    ggplot2::theme_bw() +
    # ggridges::theme_ridges(grid = TRUE) +
    ggplot2::theme(
      legend.position="none",
      axis.title.x = element.text.fig,
      axis.title.y = element.text.fig,
      axis.text.y = ggplot2::element_text(vjust = 0),
      # legend.title = element.text.fig,
      # legend.text = element.text.fig,
      strip.text.y = element.text.fig,
      strip.text.x = element.text.fig) +
    # ggplot2::facet_wrap(~PC, nrow = 1, ncol = 4, scales = "free")
    ggplot2::facet_grid(~PC, scales = "free", labeller = ggplot2::labeller(PC = facet_names))
  
  if (!is.null(plot.filename)) {
    ridges.plot.filename <- stringi::stri_join(plot.filename, "_dapc.ridges.plots.pdf")
    ggplot2::ggsave(
      filename = ridges.plot.filename, plot = res$ridges.plots,
      width = n.pop * 3, height = (n.pop * 3) / 2,
      dpi = 600, units = "cm", device = "pdf", limitsize = FALSE,
      useDingbats = FALSE)
  }
  
  
  # Clustering with snapclust --------------------------------------------------
  message("Likelihood-based genetic clustering with adegenet snapclust")
  n.k <- length(unique(x@strata[,strata]))
  k.string <- x@strata[,strata]
  k.string <- factor(k.string, levels = pop.levels, ordered = FALSE)
  clust <- adegenet::snapclust(x = x, k = n.k, pop.ini = k.string)
  # clust <- adegenet::snapclust(x = x, k = n.k, pop.ini = "ward")
  # clust <- adegenet::snapclust(x = x, k = n.k, pop.ini = "kmeans")
  # names(clust)
  # adegenet::compoplot(clust)
  clust.assignment <- data.frame(x@strata, clust$proba) %>% 
    tidyr::gather(key = ANCESTRY, value = "PERC", -c(POP_ID, INDIVIDUALS)) %>%
    dplyr::mutate(
      POP_ID = factor(POP_ID, levels = pop.levels, ordered = TRUE),
      ANCESTRY = stringi::stri_replace_all_fixed(ANCESTRY, "X", "K", vectorize_all = FALSE),
      ANCESTRY = factor(ANCESTRY, ordered = TRUE)
    ) %>%
    dplyr::arrange(POP_ID, INDIVIDUALS, ANCESTRY)
  
  n.pop <- dplyr::n_distinct(clust.assignment$POP_ID)
  
  # figure
  message("    Generating snapclust plot")
  res$clust.plot <- ggplot2::ggplot(clust.assignment, ggplot2::aes(x = INDIVIDUALS, y = PERC, fill = ANCESTRY)) +
    ggplot2::geom_bar(stat = "identity",  position = "fill") +
    ggplot2::scale_y_continuous(expand  =  c(0,  0)) + #in order for the line not to expand beyond the graph!
    ggplot2::labs(y = "Ancestry",
                  x = "Individuals",
                  title = "Likelihood-based genetic clustering with adegenet snapclust") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, family = "Helvetica", face = "bold" ),
      panel.grid.major.x  =  ggplot2::element_blank(),
      panel.grid.minor.x  =  ggplot2::element_blank(),
      panel.grid.major.y  =  ggplot2::element_line(colour = "grey60",  linetype = "dashed"),
      axis.title.x = ggplot2::element_text(size = 10,  family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 10,  family = "Helvetica", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8, family = "Helvetica"),
      axis.ticks.x = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(angle = 0),
      panel.spacing.y = ggplot2::unit(0.5, "lines"),
      panel.spacing.x = ggplot2::unit(0.1, "lines"),
      legend.position = "right"
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      label.position = "right",
      title.position = "top",
      title.hjust = 0.5,
      label.hjust = 0.5)
    ) +
    ggplot2::facet_grid(~POP_ID, scales = "free", space = "free_x")
  # clust.plot
  if (!is.null(plot.filename)) {
    clust.plot.filename <- stringi::stri_join(plot.filename, "_snapclust.plot.pdf")
    ggplot2::ggsave(
      filename = clust.plot.filename, plot = res$clust.plot,
      width = n.pop * 3, height = (n.pop * 3) / 2,
      dpi = 600, units = "cm", device = "pdf", limitsize = FALSE,
      useDingbats = FALSE)
  }
  
  # Results --------------------------------------------------------------------
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("############### grur::dapc_clustering pipeline completed ##############\n")
  options(width = opt.change)
  return(res)
}#End xval_dacp_fig


#' @title write_dapc_plot
#' @description Write the admixture dapc plots
#' @rdname write_dapc_plot
#' @keywords internal
#' @export
write_dapc_plot <- function(x, plot.name, width = 25, height = 10, dpi = 600) {
  dapc.plot.name <- stringi::stri_join(plot.name, "_dapc.admix.plot.pdf")
  message("    DAPC admixture plot: ", plot.name)
  print(x$dapc.plot)
  ggplot2::ggsave(
    filename = dapc.plot.name, plot = x$dapc.plot,
    width = width, height = height,
    dpi = dpi, units = "cm", device = "pdf", limitsize = FALSE,
    useDingbats = FALSE)
}#End write_dapc_plot


#' @title generate_dapc_pc_plot
#' @description Generate the PC plots
#' @rdname generate_dapc_pc_plot
#' @keywords internal
#' @export
#' @importFrom ggpubr ggarrange
generate_dapc_pc_plot <- function(
  pc.to.do,
  pc.data,
  variance.component,
  plot.filename = NULL
) {
  message("    Generating DAPC pc plots")
  pc_plot <- function(
    pc.to.do,
    pc.data,
    variance.component) {
    pc.plots <- list()
    
    # pc.to.do <- pc.to.do[[1]]
    pcx <- pc.to.do[1]
    pcy <- pc.to.do[2]
    element.text.fig <- ggplot2::element_text(
      size = 12, family = "Helvetica", face = "bold")
    
    plot <- ggplot2::ggplot(
      pc.data,
      ggplot2::aes_string(
        x = stringi::stri_join("LD", pcx),
        y = stringi::stri_join("LD", pcy)),
      environment = environment()) +
      ggplot2::geom_point(ggplot2::aes(colour = POP_ID), alpha = 0.5) +
      ggplot2::labs(x = stringi::stri_join("PC", pcx, " [", variance.component[pcx,2], "]")) +
      ggplot2::labs(y = stringi::stri_join("PC", pcy, " [", variance.component[pcy,2], "]")) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title.x = element.text.fig,
        axis.title.y = element.text.fig,
        legend.title = element.text.fig,
        legend.text = element.text.fig
      )
    
    plot_name <- stringi::stri_join(
      "plot.pc", pcx, ".pc", pcy)
    pc.plots[[plot_name]] <- plot
    return(pc.plots)
  }#End pc_plot
  
  pc.plots <- purrr::map(
    .x = pc.to.do, .f = pc_plot,
    pc.data = pc.data,
    variance.component = variance.component) %>%
    purrr::flatten(.)
  
  pc.plots <- ggpubr::ggarrange(
    pc.plots[[1]],
    pc.plots[[2]],
    pc.plots[[3]],
    pc.plots[[4]],
    pc.plots[[5]],
    pc.plots[[6]],
    ncol = 2, nrow = 3, legend = "right", common.legend = TRUE
  )
  
  if (!is.null(plot.filename)) {
    ggplot2::ggsave(
      filename = stringi::stri_join(plot.filename, "_dapc.pc.plots.pdf"),
      plot = pc.plots,
      width = 20, height = 15,
      dpi = 600, units = "cm",
      useDingbats = FALSE, device = "pdf", limitsize = FALSE)
  }
  return(pc.plots)
}#End generate_dapc_pc_plot

