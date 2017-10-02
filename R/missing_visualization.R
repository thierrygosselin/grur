#' @name missing_visualization
#' @title Visualize missing genotypes in genomic data set
#' @description Use this function to visualize pattern of missing data.
#'  \itemize{
#'    \item \strong{Input file:} various file formats are supported
#'    (see \code{data} argument below).
#'    \item \strong{Filters:} genotypes, markers, individuals and populations can be
#'   filtered and/or selected in several ways using blacklist,
#'   whitelist and other arguments.
#'    \item \strong{IBM-PCoA:} conduct identity-by-missingness analyses using
#'    Principal Coordinates Analysis, PCoA (also called Multidimensional Scaling, MDS).
#'    \item \strong{FH measure vs missingness}: missingness at the individual level
#'    is contrasted against FH, a new measure of IBDg
#'    (Keller et al., 2011; Kardos et al., 2015; Hedrick & Garcia-Dorado, 2016)
#'    FH is based on the excess in the observed number of homozygous
#'    genotypes within an individual relative to the mean number of homozygous
#'    genotypes expected under random mating.
#'    IBDg is a proxy measure of the realized proportion of the genome
#'    that is identical by descent.
#'    Within this function, we're using a modified version of the measure
#'    described in (Keller et al., 2011; Kardos et al., 2015).
#'    The new measure is population-wise and tailored for RADseq data
#'    (see \code{\link[radiator]{ibdg_fh}} for details).
#'    \item \strong{Figures and Tables:} figures and summary tables of
#'    missing information at the marker, individual and population level are
#'    generated.
#'    \item \strong{Whitelist:} create whitelists of markers based on
#'    desired thresholds of missing genotypes.
#'    \item \strong{Blacklists:} create blacklists of individuals based on
#'    desired thresholds of missing genotypes.
#'    \item \strong{Tidy data:} if the filename argument is used, the
#'    function also output the data in a tidy format.
#' }

#' @inheritParams radiator::tidy_genomic_data
#' @param strata (optional/required) Required for VCF and haplotypes files,
#' optional for the other formats supported.
#'
#' The strata file is a tab delimited file with a minimum of 2 columns headers:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' If a \code{strata} file is specified with all file formats that don't
#' require it, the strata argument will have precedence on the population
#' groupings used internally in those file formats.
#' The \code{STRATA} column can be any hierarchical grouping.
#' To create a strata file see \code{\link[radiator]{individuals2strata}}.
#' If you have already run
#' \href{http://catchenlab.life.illinois.edu/stacks/}{stacks} on your data,
#' the strata file is similar to a stacks \emph{population map file}, make sure you
#' have the required column names (\code{INDIVIDUALS} and \code{STRATA}).
#' For \code{missing_visualization} function, use additional columns in the strata
#' file to store metadata that you want to look for pattern of missingness.
#' e.g. lanes, chips, sequencers, etc.
#' Default: \code{strata = NULL}.

#' @param strata.select (optional, character) Use this argument to select the column
#' from the strata file to generate the PCoA-IBM plot. More than 1 column you
#' want to visualize, use a string of character
#' e.g. \code{strata.select = c("POP_ID", "LANES", "SEQUENCER", "WATERSHED")} to test
#' 4 grouping columns inside the \code{strata} file.
#' Default: \code{strata.select = "POP_ID"}

#' @param distance.method (character) The distance measure to be used.
#' This must be one of "euclidean", "maximum", "manhattan", "canberra",
#' "binary" or "minkowski". The function uses \code{\link[stats]{dist}}.
#' Default: \code{distance.method = "euclidean"}.

#' @param ind.missing.geno.threshold (string) Percentage of missing genotype
#' allowed per individuals (to create the blacklists).
#' Default:\code{ind.missing.geno.threshold = c(10, 20, 30, 40, 50, 60, 70, 80, 90)}.

#' @param filename (optional) Name of the tidy data set,
#' written to the directory created by the function.

#' @param write.plot (optional, logical) When \code{write.plot = TRUE}, the function
#' will write to the directory created by the function the plots, except the heatmap
#' that take longer to generate. For this, do it manually following example below.
#' Default: \code{write.plot = TRUE}.

#' @return A list is created with several objects: the tidy data,
#' the principal coordinates
#' with eigenvalues of the PCoA, the identity-by-missingness plot, several
#' summary tables and plots of missing information
#' per individuals, populations and markers. Blacklisted ids are also included.
#' Whitelists of markers with different level of missingness are also generated
#' automatically.
#' A heatmap showing the missing values in black and genotypes in grey provide a
#' general overview of the missing data. The heatmap is usually long to generate,
#' and thus, it's just included as an object in the list and not written in the folder.

#' @examples
#' \dontrun{
#' Using a  VCF file, the simplest for of the function:
#' ibm.koala <- missing_visualization(
#' data = "batch_1.vcf",
#' strata = "population.map.strata.tsv"
#' )
#' # To see what's inside the list
#' names(ibm.koala)
#' # To view the heatmap:
#' ibm.koala$heatmap
#' # To view the IBM analysis plot:
#' ibm.koala$ibm_plot
#' }

#' @references Legendre, P. and Legendre, L. (1998) Numerical Ecology,
#' 2nd English edition. Amsterdam: Elsevier Science BV.
#' @references Keller MC, Visscher PM, Goddard ME (2011)
#' Quantification of inbreeding due to distant ancestors and its detection
#'  using dense single nucleotide polymorphism data. Genetics, 189, 237–249.
#' @references Kardos M, Luikart G, Allendorf FW (2015)
#' Measuring individual inbreeding in the age of genomics: marker-based
#' measures are better than pedigrees. Heredity, 115, 63–72.
#' @references Hedrick PW, Garcia-Dorado A. (2016)
#' Understanding Inbreeding Depression, Purging, and Genetic Rescue.
#' Trends in Ecology and Evolution. 2016;31: 940-952.
#' doi:10.1016/j.tree.2016.09.005

#' @export
#' @rdname missing_visualization
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid geom_histogram aes_string scale_fill_manual theme_bw stat_smooth geom_boxplot ggsave
#' @importFrom dplyr distinct rename arrange mutate select summarise group_by ungroup filter inner_join left_join
#' @importFrom stringi stri_join stri_replace_all_fixed stri_replace_all_regex
#' @importFrom utils count.fields
#' @importFrom readr read_tsv write_tsv
#' @importFrom data.table fread melt.data.table as.data.table
#' @importFrom ape pcoa
#' @importFrom stats dist lm
#' @importFrom tibble data_frame
#' @importFrom radiator tidy_genomic_data change_pop_names ibdg_fh detect_all_missing
#' @importFrom cowplot plot_grid align_plots

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

missing_visualization <- function(
  data,
  vcf.metadata = FALSE,
  strata = NULL,
  strata.select = "POP_ID",
  distance.method = "euclidean",
  ind.missing.geno.threshold = c(10, 20, 30, 40, 50, 60, 70, 80, 90),
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  blacklist.id = NULL,
  blacklist.genotype = NULL,
  whitelist.markers = NULL,
  monomorphic.out = TRUE,
  max.marker = NULL,
  snp.ld = NULL,
  common.markers = FALSE,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  write.plot = TRUE
) {
  opt.change <- getOption("width")
  options(width = 70)
  cat("#######################################################################\n")
  cat("#################### grur::missing_visualization ######################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  res <- list() #empty list to store results
  
  # manage missing arguments -----------------------------------------------------
  if (missing(data)) stop("Input file missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  # Date and time --------------------------------------------------------------
  file.date <- stringi::stri_replace_all_fixed(
    Sys.time(),
    pattern = " EDT", replacement = "") %>%
    stringi::stri_replace_all_fixed(
      str = .,
      pattern = c("-", " ", ":"), replacement = c("", "@", ""),
      vectorize_all = FALSE) %>%
    stringi::stri_sub(str = ., from = 1, to = 13)
  
  path.folder.message <- stringi::stri_join("missing_visualization_", file.date, sep = "")
  path.folder <- stringi::stri_join(getwd(),"/", "missing_visualization_", file.date, sep = "")
  dir.create(file.path(path.folder))
  
  message("Folder created: \n", path.folder.message)
  file.date <- NULL #unused object
  # import data ----------------------------------------------------------------
  message("\nImporting data")
  if (!is.null(filename)) filename <- stringi::stri_join(path.folder, "/", filename)
  res$tidy.data <- radiator::tidy_genomic_data(
    data = data,
    vcf.metadata = FALSE,
    blacklist.id = blacklist.id,
    blacklist.genotype = blacklist.genotype,
    whitelist.markers = whitelist.markers,
    monomorphic.out = monomorphic.out,
    max.marker = max.marker,
    snp.ld = snp.ld,
    common.markers = common.markers,
    strata = strata,
    pop.select = pop.select,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    filename = filename,
    parallel.core = parallel.core,
    verbose = FALSE
  )
  
  # strata.df --------------------------------------------------------
  strata.df <- dplyr::ungroup(res$tidy.data) %>%
    dplyr::select(-dplyr::one_of(c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT",
                                   "GT_VCF", "GT_VCF_NUC", "GT", "GT_BIN"))) %>% 
    dplyr::distinct(INDIVIDUALS, .keep_all = TRUE)
  
  # New column with GT_MISSING_BINARY O for missing 1 for not missing...
  res$tidy.data <- dplyr::mutate(
    .data = res$tidy.data,
    GT_MISSING_BINARY = dplyr::if_else(GT == "000000", "0", "1"),
    GT_MISSING_BINARY = as.numeric(GT_MISSING_BINARY)
  )
  
  # some statistics  -----------------------------------------------------------
  message("\nInformations:")
  strata.stats <- strata.df %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::tally(.) %>%
    dplyr::mutate(STRATA = stringi::stri_join(POP_ID, n, sep = " = "))
  
  duplicate.id <- nrow(strata.df) - length(unique(strata.df$INDIVIDUALS))
  
  n.chrom = dplyr::n_distinct(res$tidy.data$CHROM)
  n.locus = dplyr::n_distinct(res$tidy.data$LOCUS)
  n.snp = dplyr::n_distinct(res$tidy.data$MARKERS)
  
  # output the proportion of missing genotypes
  na.before <- dplyr::summarise(
    .data = res$tidy.data,
    MISSING = round(length(GT_MISSING_BINARY[GT_MISSING_BINARY == 0])/length(GT_MISSING_BINARY), 6)) %>%
    purrr::flatten_dbl(.) %>% format(., scientific = FALSE)
  
  marker.problem <- radiator::detect_all_missing(data = res$tidy.data)
  if (marker.problem$marker.problem) {
    res$tidy.data <- marker.problem$data
    message("Marker problem: some marker(s) are missing all genotypes")
    message("    removing markers, see blacklist for details: blacklist.markers.all.missing.tsv")
    res$blacklist.markers.all.missing <- marker.problem$blacklist.markers.all.missing
  }
  marker.problem <-  NULL
  
  n.pop <- dplyr::n_distinct(strata.df$POP_ID)
  n.ind <- dplyr::n_distinct(strata.df$INDIVIDUALS)
  message("Number of populations: ", n.pop)
  message("Number of individuals: ", n.ind)
  message("Number of ind/pop:\n", stringi::stri_join(strata.stats$STRATA, collapse = "\n"))
  message("\nNumber of duplicate id: ", duplicate.id)
  message("Number of chrom/scaffolds: ", n.chrom)
  message("Number of locus: ", n.locus)
  message("Number of SNPs: ", n.snp)
  
  message("\nProportion of missing genotypes: ", na.before)

  # Identity-by-missingness (IBM) analysis -------------------------------------
  # MultiDimensional Scaling analysis (MDS) - Principal Coordinates Analysis (PCoA)
  message("\n\nIdentity-by-missingness (IBM) analysis using\n    Principal Coordinate Analysis (PCoA)...")
  
  input.pcoa <- res$tidy.data %>%
    dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT_MISSING_BINARY) %>%
    dplyr::group_by(POP_ID, INDIVIDUALS) %>%
    tidyr::spread(data = ., key = MARKERS, value = GT_MISSING_BINARY)
  
  # we need rownames for this
  suppressWarnings(rownames(input.pcoa) <- input.pcoa$INDIVIDUALS)
  input.pcoa <-  dplyr::ungroup(input.pcoa) %>% dplyr::select(-POP_ID, -INDIVIDUALS)
  
  # euclidean distances between the rows
  # distance.method <- "euclidean"
  # d <- stats::dist(x = input.pcoa, method = distance.method)
  
  # alternative tested
  # d <- vegan::vegdist(x = input.pcoa, method = distance.method) # longer than stats::dist
  
  # for metric PCoA/MDS
  # ibm <- ape::pcoa(D = d) 
  ibm <- ape::pcoa(
    D = stats::dist(
    x = input.pcoa, method = distance.method)) #Legendre's ape
  input.pcoa <- NULL
  
  # ibm$correction # test
  # ibm$note # test
  # ibm$values # test
  # ibm$vectors # test
  # ibm$trace # test
  
  # Should broken_stick values be reported?
  
  # variance
  variance.component <-  tibble::data_frame(EIGENVALUES = ibm$values$Eigenvalues) %>%
    dplyr::mutate(
      VARIANCE_PROP = round(EIGENVALUES/sum(EIGENVALUES), 2)
    )
  
  
  # alternative tested giving the same results:
  # ibm <- stats::cmdscale(d, eig = TRUE, k = 2)
  
  # for non-metric PCoA/MDS
  # ibm <- MASS::isoMDS(d, k=2) # k is the number of dim
  
  # alternative: sammon
  # ibm <- MASS::sammon(d, k =2)
  
  # all gives the same results...
  
  # prep the data for figure:
  # for MASS::isoMDS and stats::cmdscale
  # res$vectors <- tibble::data_frame(INDIVIDUALS = rownames(ibm$points), V1 = ibm$points[,1], V2 = ibm$points[,2]) %>%
  #   dplyr::inner_join(strata.df, by = "INDIVIDUALS")
  # res$vectors <- tibble::data_frame(INDIVIDUALS = rownames(ibm$vectors), V1 = ibm$vectors[,1], V2 = ibm$vectors[,2]) %>%
  # dplyr::inner_join(strata.df, by = "INDIVIDUALS")
  
  # with vegan and ape
  
  # integrate missingness per individuals
  res$missing.genotypes.ind <- res$tidy.data %>%
    dplyr::select(MARKERS, INDIVIDUALS, POP_ID, GT_MISSING_BINARY) %>%
    dplyr::group_by(INDIVIDUALS, POP_ID) %>%
    dplyr::summarise(
      MISSING_GENOTYPE = length(GT_MISSING_BINARY[GT_MISSING_BINARY == 0]),
      MARKER_NUMBER = length(MARKERS),
      MISSING_GENOTYPE_PROP = MISSING_GENOTYPE/MARKER_NUMBER,
      PERCENT = round((MISSING_GENOTYPE_PROP)*100, 2)
    ) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(POP_ID, INDIVIDUALS)
  
  strata.missing <- dplyr::inner_join(
    strata.df,
    dplyr::select(res$missing.genotypes.ind, INDIVIDUALS, MISSING_GENOTYPE_PERCENT = PERCENT),
    by = "INDIVIDUALS")
    
    
  res$vectors <- dplyr::inner_join(
    strata.missing,
    tibble::rownames_to_column(df = data.frame(ibm$vectors), var = "INDIVIDUALS")
    , by = "INDIVIDUALS"
  )
  
  ibm <- NULL
  
  # adjust pop_id
  res$vectors <- radiator::change_pop_names(data = res$vectors, pop.levels = pop.labels, pop.labels = pop.labels)
  
  
  
  
  message("Generating Identity by missingness plot")
  pc.to.do <- utils::combn(1:4, 2, simplify = FALSE)
  
  # if (length(strata.select) > 1) {
    res$ibm.plots <- purrr::map(
      .x = strata.select, 
      .f = generate_pcoa_plot,
      pc.to.do = pc.to.do,
      vectors = res$vectors,
      variance.component = variance.component,
      path.folder = path.folder, write.plot = write.plot) %>%
      purrr::flatten(.)
  # } else {
  #   res$ibm.plots <- generate_pcoa_plot(
  #     strata.select = strata.select,
  #     pc.to.do = pc.to.do,
  #     vectors = res$vectors,
  #     variance.component = variance.component,
  #     path.folder = path.folder, write.plot = write.plot)
  # }
    variance.component <- pc.to.do <- NULL
  
  # Eric's code here -----------------------------------------------------------
    # Note to Eric: I think it's fits very here ...
    # This is a refinement of the previous summaries where we want to see if the
    # percent missing for each level of a given factor
    # (population, plate, region, etc) is related to the number missing at a
    # locus. I’m starting to see the issue of missingness as being decomposable
    # into influences from external factors and from the loci themselves.
    
    # For instance, if a given population tends to have missing data
    # (say due to poor quality samples), then when there is missing data at a
    # locus, it should be attributable to a population regardless of how many
    # samples are missing and the frequency should be about the same.
    # If there is something about a locus, it should be missing at the same
    # rate across all populations, and there should be more missing than normal.
    # This function looks only at loci where there are missing data and calculates
    # the median percent missing for all loci with a given total number of
    # samples missing (e.g., all loci where 1 sample is missing, all loci
    # where 2 samples are missing, etc.). It calculates this for each level of
    # a given factor (i.e., population). This median is plotted as a function of
    # the total number missing. Overlayed on that plot is the percent expected
    # at random (black line) and the overall percent missing in that population
    # (red line). Note that these two values are percents of missing genotypes
    # (# of genotypes missing in population / total number of missing genotypes).
    message("Analysing percentage missing ...") # I'll leave better message to you
    res$pct.missing.total <- purrr::map(
      .x = strata.select,
      .f = pct_missing_by_total,
      data = res$tidy.data,
      ci = 0.95,
      path.folder = path.folder,
      write.plot = write.plot) %>%
      purrr::flatten(.)
    
    # Heatmap --------------------------------------------------------------------
  res$heatmap <- res$tidy.data %>%
    dplyr::mutate(
      GT_MISSING_BINARY = as.character(GT_MISSING_BINARY),
      Missingness = stringi::stri_replace_all_regex(
        GT_MISSING_BINARY,
        pattern = c("^0$", "^1$"),
        replacement = c("missing", "genotyped"),
        vectorize_all = FALSE)) %>%
    ggplot2::ggplot(data = .,(ggplot2::aes(y = MARKERS, x = as.character(INDIVIDUALS)))) +
    ggplot2::geom_tile(ggplot2::aes(fill = Missingness)) +
    ggplot2::scale_fill_manual(values = c("grey", "black")) +
    ggplot2::labs(y = "Markers") +
    ggplot2::labs(x = "Individuals") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    ) +
    ggplot2::facet_grid(~POP_ID, scales = "free", space = "free_x")
  # res$heatmap
  
  # Missing summary ------------------------------------------------------------
  message("Generating missing information summary tables and plot")
  
  # Individuals-----------------------------------------------------------------
  message("Missingness per individuals")
  res$missing.genotypes.ind <- res$tidy.data %>%
    dplyr::select(MARKERS, INDIVIDUALS, POP_ID, GT_MISSING_BINARY) %>%
    dplyr::group_by(INDIVIDUALS, POP_ID) %>%
    dplyr::summarise(
      MISSING_GENOTYPE = length(GT_MISSING_BINARY[GT_MISSING_BINARY == 0]),
      MARKER_NUMBER = length(MARKERS),
      MISSING_GENOTYPE_PROP = MISSING_GENOTYPE/MARKER_NUMBER,
      PERCENT = round((MISSING_GENOTYPE_PROP)*100, 2)
    ) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(POP_ID, INDIVIDUALS)
  
  # Figures
  axis.title.element.text.fig <- ggplot2::element_text(
    size = 12, family = "Helvetica", face = "bold")
  axis.text.element.text.fig <- ggplot2::element_text(
    size = 10, family = "Helvetica")

  # manhattan and violin plots
  missing.genotypes.ind.plots <- suppressMessages(ggplot2::ggplot(
    data = res$missing.genotypes.ind,
    ggplot2::aes(x = POP_ID, y = MISSING_GENOTYPE_PROP, colour = POP_ID)) +
    ggplot2::geom_jitter(alpha = 0.5) +
    ggplot2::geom_violin(trim = TRUE, fill = NA) +
    ggplot2::geom_boxplot(width = 0.1, fill = NA, outlier.colour = NA, outlier.fill = NA) +
    ggplot2::labs(y = "Individual's missing genotypes (proportion)") +
    ggplot2::labs(x = "Populations") +
    ggplot2::labs(colour = "Populations") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.title.x = axis.title.element.text.fig,
      axis.text.x = axis.text.element.text.fig,
      axis.title.y = axis.title.element.text.fig,
      axis.text.y = axis.text.element.text.fig
    ) +
    ggplot2::coord_flip())
  # missing.genotypes.ind.plots
  
  # histogram
  missing.genotypes.ind.histo <- suppressMessages(ggplot2::ggplot(
    data = res$missing.genotypes.ind,
    ggplot2::aes(x = MISSING_GENOTYPE_PROP)) +
    ggplot2::geom_histogram() +
    ggplot2::labs(x = "Individual's missing genotypes (proportion)") +
    ggplot2::labs(y = "Individuals (number)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = axis.title.element.text.fig,
      axis.title.y = axis.title.element.text.fig,
      axis.text.x = axis.text.element.text.fig,
      axis.text.y = axis.text.element.text.fig,
    ))
  # missing.genotypes.ind.histo
  
  # helper plot for individual's genotyped threshold
  ind.genotyped.helper.plot <- ind_genotyped_helper(res$missing.genotypes.ind)
  
  
  # merge plots
  plots <- cowplot::align_plots(
    missing.genotypes.ind.plots, ind.genotyped.helper.plot, align = "v", axis = "l")
  top.row.plot <- suppressMessages(cowplot::plot_grid(
    plots[[1]], missing.genotypes.ind.histo, labels = c("A", "B"), align = "h"))
  
  res$missing.genotypes.ind.plots <- suppressMessages(cowplot::plot_grid(
    top.row.plot, plots[[2]],
    ncol = 1, nrow = 2, labels = c("", "C"), rel_heights = c(1.5, 1)))
  plots <- top.row.plot <- ind.genotyped.helper.plot <- NULL
  
  if (write.plot) {
    ggplot2::ggsave(
      filename = stringi::stri_join(path.folder, "/missing.genotypes.ind.plots.pdf"),
      plot = res$missing.genotypes.ind.plots,
      width = n.pop * 4, height = n.pop * 2.5,
      dpi = 600, units = "cm", useDingbats = FALSE)
  }
  # res$missing.genotypes.ind.plots
  
  
  # Blacklist individuals-------------------------------------------------------
  # ind.missing.geno.threshold <- c(10, 20,30,50,60,70)
  blacklists <- purrr::map(
    .x = ind.missing.geno.threshold,
    .f = blacklists_id_generator,
    y = res$missing.genotypes.ind,
    path.folder = path.folder) %>% purrr::flatten(.)
  
  
  message("Blacklist(s) of individuals generated: ", length(blacklists))
  blacklists.stats <- purrr::map_df(.x = blacklists, .f = nrow) %>% 
    tidyr::gather(data = ., key = "BLACKLIST", value = "n") %>%
    dplyr::transmute(BLACKLIST = stringi::stri_join(BLACKLIST, n, sep = " = "))
  message("    Number of individual(s) blacklisted per blacklist generated:\n", stringi::stri_join("    ", blacklists.stats$BLACKLIST, collapse = "\n"))
  res <- c(res, blacklists)
  blacklists.stats <- blacklists <- NULL
  
  # for (i in ind.missing.geno.threshold) {
  #   # i <- 30
  #   blacklist.id.missing.geno <- res$missing.genotypes.ind %>%
  #     dplyr::filter(PERCENT >= i) %>%
  #     dplyr::ungroup(.) %>%
  #     dplyr::select(INDIVIDUALS)
  #   if (length(blacklist.id.missing.geno$INDIVIDUALS) > 0) {
  #     blacklist_name <- stringi::stri_join("blacklist.id.missing.", i)
  #     res[[blacklist_name]] <- blacklist.id.missing.geno
  #     readr::write_tsv(
  #       blacklist.id.missing.geno,
  #       stringi::stri_join(path.folder, "/", blacklist_name, ".tsv"))
  #   }
  # }
  
  
  
  # FH -------------------------------------------------------------------------
  message("Calulation of FH: a measure of IBDg")
  fh <- radiator::ibdg_fh(data = res$tidy.data, verbose = FALSE)
  # fh <- ibdg_fh(data = res$tidy.data, verbose = FALSE)
  
  res$missing.genotypes.ind.fh <- suppressWarnings(
    dplyr::full_join(
      res$missing.genotypes.ind,
      fh$fh
      # dplyr::select(.data = fh, INDIVIDUALS, FH)
      , by = c("INDIVIDUALS", "POP_ID")
    )
  )
  
  res$missing.genotypes.ind.fh.plots <- ggplot2::ggplot(
    res$missing.genotypes.ind.fh, ggplot2::aes(y = FH, x = MISSING_GENOTYPE_PROP)) +
    ggplot2::geom_point() +
    ggplot2::stat_smooth(method = stats::lm, level = 0.99) +
    # labs(title = "Correlation between missingness and inbreeding coefficient") +
    ggplot2::labs(y = "Individual IBDg (FH)") +
    ggplot2::labs(x = "Missing genotype (proportion)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    )
  
  # merge plots
  plots <- cowplot::align_plots(
    res$missing.genotypes.ind.fh.plots, fh$fh.manhattan.box.plot, align = "v", axis = "l")
  bottom.row.plot <- suppressMessages(cowplot::plot_grid(
    plots[[2]], fh$fh.distribution.plot, labels = c("B", "C"), align = "h"))
  
  res$missing.genotypes.ind.fh.plots <- cowplot::plot_grid(
    plots[[1]], bottom.row.plot,
    ncol = 1, nrow = 2, labels = c("A", ""), rel_heights = c(1.5, 1))
  # res$missing.genotypes.ind.fh.plots
  fh <- plots <- bottom.row.plot <- NULL
  
  if (write.plot) {
    ggplot2::ggsave(
      filename = stringi::stri_join(path.folder, "/missing.genotypes.ind.fh.plots.pdf"),
      plot = res$missing.genotypes.ind.fh.plots,
      width = n.pop * 4, height = n.pop * 2.5,
      dpi = 600, units = "cm", useDingbats = FALSE)
  }
  
  # Populations-----------------------------------------------------------------
  message("Missingness per populations")
  res$missing.genotypes.pop <- res$missing.genotypes.ind %>%
    dplyr::select(INDIVIDUALS, POP_ID, MISSING_GENOTYPE_PROP, PERCENT) %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::summarise(
      MISSING_GENOTYPE_PROP = mean(MISSING_GENOTYPE_PROP, na.rm = TRUE),
      PERCENT = round(MISSING_GENOTYPE_PROP, 2)
    )
  
  # Markers---------------------------------------------------------------------
  message("Missingness per markers")
  
  res$missing.genotypes.markers.overall <- res$tidy.data %>%
    dplyr::select(MARKERS, INDIVIDUALS, GT_MISSING_BINARY) %>%
    dplyr::group_by(MARKERS) %>%
    dplyr::summarise(
      MISSING_GENOTYPE = length(GT_MISSING_BINARY[GT_MISSING_BINARY == 0]),
      INDIVIDUALS_NUMBER = length(INDIVIDUALS),
      MISSING_GENOTYPE_PROP = MISSING_GENOTYPE / INDIVIDUALS_NUMBER,
      PERCENT = round(MISSING_GENOTYPE_PROP * 100, 2)
    ) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(MARKERS)
  
  res$missing.genotypes.markers.pop <- res$tidy.data %>%
    dplyr::select(MARKERS, INDIVIDUALS, POP_ID, GT_MISSING_BINARY) %>%
    dplyr::group_by(MARKERS, POP_ID) %>%
    dplyr::summarise(
      MISSING_GENOTYPE = length(GT_MISSING_BINARY[GT_MISSING_BINARY == 0]),
      INDIVIDUALS_NUMBER = length(INDIVIDUALS),
      MISSING_GENOTYPE_PROP = MISSING_GENOTYPE / INDIVIDUALS_NUMBER,
      PERCENT = round(MISSING_GENOTYPE_PROP * 100, 2)
    ) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(POP_ID, MARKERS)
  
  markers.missing.geno.threshold <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 ,0.8, 0.9)
  whitelists <- purrr::map(
    .x = markers.missing.geno.threshold,
    .f = whitelists_markers_generator,
    y = res$missing.genotypes.markers.overall,
    path.folder = path.folder) %>% purrr::flatten(.)
  
  
  message("Whitelist(s) of markers generated: ", length(whitelists))
  whitelists.stats <- purrr::map_df(.x = whitelists, .f = nrow) %>% 
    tidyr::gather(data = ., key = "WHITELIST", value = "n") %>%
    dplyr::transmute(WHITELIST = stringi::stri_join(WHITELIST, n, sep = " = "))
  message("    Number of markers whitelisted per whitelist generated:\n", stringi::stri_join("    ", whitelists.stats$WHITELIST, collapse = "\n"))
  res <- c(res, whitelists)
  whitelists.stats <- whitelists <- NULL
  
  # Figure markers
  
  # violin plots
  missing.genotypes.markers.plots <- suppressMessages(ggplot2::ggplot(
    data = res$missing.genotypes.markers.pop,
    ggplot2::aes(x = POP_ID, y = MISSING_GENOTYPE_PROP, colour = POP_ID)) +
    # ggplot2::geom_jitter(alpha = 0.5) +
    ggplot2::geom_violin(trim = TRUE, fill = NA) +
    ggplot2::geom_boxplot(width = 0.1, fill = NA, outlier.colour = NA, outlier.fill = NA) +
    ggplot2::labs(y = "Marker's missing genotypes (proportion)") +
    ggplot2::labs(x = "Populations") +
    ggplot2::labs(colour = "Populations") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, family = "Helvetica"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10, family = "Helvetica")
    ) +
    ggplot2::coord_flip())
  # missing.genotypes.markers.plots
  
  markers.histo <- suppressMessages(ggplot2::ggplot(
    data = res$missing.genotypes.markers.overall, ggplot2::aes(x = MISSING_GENOTYPE_PROP)) +
    ggplot2::geom_histogram() +
    ggplot2::labs(x = "Marker's missing genotypes (proportion)") +
    ggplot2::labs(y = "Markers (number)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, family = "Helvetica"),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    ))
  # markers.histo
  
  # helper plot for markers's genotyped threshold
  markers.helper.plot <- markers_genotyped_helper(
    x = res$missing.genotypes.markers.pop, y = res$missing.genotypes.markers.overall)
  
  # merge plots
  plots <- cowplot::align_plots(
    missing.genotypes.markers.plots, markers.helper.plot, align = "v", axis = "l")
  top.row.plot <- suppressMessages(cowplot::plot_grid(
    plots[[1]], markers.histo, labels = c("A", "B"), align = "h"))
  
  res$missing.genotypes.markers.plots <- suppressMessages(cowplot::plot_grid(
    top.row.plot, plots[[2]],
    ncol = 1, nrow = 2, labels = c("", "C"), rel_heights = c(1.5, 1)))
  # res$missing.genotypes.markers.plots
  plots <- top.row.plot <- markers.helper.plot <- markers.histo <- NULL
  
  if (write.plot) {
    ggplot2::ggsave(
      filename = stringi::stri_join(path.folder, "/missing.genotypes.markers.plots.pdf"),
      plot = res$missing.genotypes.markers.plots,
      width = n.pop * 4, height = n.pop * 2.5,
      dpi = 600, units = "cm", useDingbats = FALSE)
  }
  
  # res$missing.genotypes.markers.plot
  
  # # Missingness per markers and per populations
  # message("Missingness per markers and populations")
  #
  # missing.genotypes.markers.pop <- dplyr::ungroup(res$tidy.data) %>%
  #   dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT_MISSING_BINARY) %>%
  #   dplyr::group_by(MARKERS, POP_ID, GT_MISSING_BINARY) %>%
  #   dplyr::tally(.) %>%
  #   dplyr::ungroup(.) %>%
  #   tidyr::complete(
  #     data = .,
  #     GT_MISSING_BINARY,
  #     nesting(MARKERS, POP_ID),
  #     fill = list(n = 0)
  #   ) %>%
  #   dplyr::group_by(MARKERS, POP_ID) %>%
  #   dplyr::summarise(MISSING_GENOTYPE_PROP = n[GT_MISSING_BINARY == 0] / sum(n))
  
  # Results --------------------------------------------------------------------
  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  cat("################ grur::missing_visualization completed ################\n")
  options(width = opt.change)
  return(res)
}
