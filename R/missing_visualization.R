#' @name missing_visualization
#' @title Visualize missing genotypes in genomic data set
#' @description Use this function to visualize pattern of missing data.
#'  \itemize{
#'    \item \strong{Input file:} various file formats are supported
#'    (see \code{data} argument below).
#'    \item \strong{IBM-PCoA:} conduct identity-by-missingness analyses using
#'    Principal Coordinates Analysis, PCoA (also called Multidimensional Scaling, MDS).
#'    \item \strong{RDA}: Redundancy Analysis using the strata provided to test the 
#'    null hypothesis of no pattern of missingness between strata.
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
#'    \item \strong{Whitelists:} create whitelists of markers based on
#'    desired thresholds of missing genotypes.
#'    \item \strong{Blacklists:} create blacklists of individuals based on
#'    desired thresholds of missing genotypes.
#'    \item \strong{Tidy data:} if the filename argument is used, the
#'    function also output the data in the directory in a tidy format (see details).
#' }

#' @inheritParams radiator::tidy_genomic_data

#' @param strata (optional)
#' The strata file is a tab delimited file with a minimum of 2 columns headers:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' If a \code{strata} file is specified the strata argument will have 
#' precedence on the population groupings (\code{POP_ID}) used internally.
#' The \code{STRATA} column can be any hierarchical grouping.
#' For \code{missing_visualization} function, use additional columns in the strata
#' file to store metadata that you want to look for pattern of missingness.
#' e.g. lanes, chips, sequencers, etc.
#' Note that you need different values inside the \code{STRATA} for the function
#' to work.
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

#' @param ... (optional) Advance mode that allows to pass further arguments
#' for fine-tuning the function. Also used for legacy arguments (see advance mode or
#' special sections below).


#' @return A list is created with several objects:
#' the principal coordinates
#' with eigenvalues of the PCoA, the identity-by-missingness plot, several
#' summary tables and plots of missing information
#' per individuals, populations and markers. Blacklisted ids are also included.
#' Whitelists of markers with different level of missingness are also generated
#' automatically.
#' A heatmap showing the missing values in black and genotypes in grey provide a
#' general overview of the missing data. The heatmap is usually long to generate,
#' and thus, it's just included as an object in the list and not written in the folder.

#' @details
#' \strong{filename}
#' 
#' The function uses \code{\link[fst]{write.fst}},
#' to write the tidy data frame in the directory.
#' The file extension appended to
#' the \code{filename} provided is \code{.rad}.
#' The file is written with the
#' \href{https://github.com/fstpackage/fst}{Lightning Fast Serialization of Data Frames for R} package.
#' To read the tidy data file back in R use \code{\link[fst]{read.fst}}.

#' @examples
#' \dontrun{
#' #Using a  VCF file, the simplest for of the function:
#' ibm.koala <- missing_visualization(
#' data = "batch_1.vcf",
#' strata = "population.map.strata.tsv"
#' )
#' 
#' # To see what's inside the list
#' names(ibm.koala)
#' 
#' # To view the heatmap:
#' ibm.koala$heatmap
#' 
#' # To save the heatmap
#' # move to the appropriate directory
#' ggplot2::ggsave(
#' filename = "heatmap.missing.pdf",
#' plot = ibm.koala$heatmap,
#' width = 15, height = 20,
#' dpi = 600, units = "cm", useDingbats = FALSE)
#' 
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
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com} and Eric Archer \email{eric.archer@@noaa.gov}

missing_visualization <- function(
  data,
  strata = NULL,
  strata.select = "POP_ID",
  distance.method = "euclidean",
  ind.missing.geno.threshold = c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90),
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  write.plot = TRUE,
  ...
) {
  # testing
  # path.folder = NULL
  
  verbose <- TRUE
  cat("################################################################################\n")
  cat("######################## grur::missing_visualization ###########################\n")
  cat("################################################################################\n")
  # Cleanup---------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date/time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()# for timing
  res <- list() # store imputed data and output...
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(timing <- proc.time() - timing, add = TRUE)
  on.exit(if (verbose) message("\nComputation time, overall: ", round(timing[[3]]), " sec"), add = TRUE)
  on.exit(if (verbose) cat("############################ missing_visualization #############################\n"), add = TRUE)
  
  
  # Function call and dotslist -------------------------------------------------
  rad.dots <- radiator::radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = "path.folder",
    deprecated = c(
      "pop.levels", "pop.labels", "pop.select", "blacklist.id", "blacklist.genotype",
      "common.markers", "monomorphic.out", "snp.ld"
    ),
    verbose = TRUE
  )
  
  # For testing
  
  # vcf.metadata = FALSE
  # strata = NULL
  # strata.select = "POP_ID"
  # distance.method = "euclidean"
  # ind.missing.geno.threshold = c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90)
  # filename = NULL
  # parallel.core = parallel::detectCores() - 1
  # write.plot = TRUE
  # path.folder = NULL
  # 
  
  # manage missing arguments -----------------------------------------------------
  if (missing(data)) rlang::abort("Input file missing")
  # Folders---------------------------------------------------------------------
  path.folder <- radiator::generate_folder(
    f = path.folder,
    rad.folder = "missing_visualization",
    internal = FALSE,
    file.date = file.date,
    prefix_int = FALSE,
    verbose = verbose)
  if (!is.null(filename)) filename <- file.path(path.folder, filename)
  
  # write the dots file
  radiator::write_rad(
    data = rad.dots,
    path = path.folder,
    filename = stringi::stri_join("grur_missing_visualization_data_args_", file.date, ".tsv"),
    tsv = TRUE,
    internal = FALSE,
    verbose = verbose
  )
  
  # import data ----------------------------------------------------------------
  message("\nImporting data")
  data.type <- radiator::detect_genomic_format(data)
  if (data.type %in% c("SeqVarGDSClass", "gds.file", "vcf.file")) {
    # GDS & VCF
    if (data.type %in% c("SeqVarGDSClass", "gds.file", "vcf.file")) {
      want <- unique(c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID", "GT_BIN", strata.select))
      
      if (!"SeqVarTools" %in% utils::installed.packages()[,"Package"]) {
        rlang::abort('Please install SeqVarTools for this option:\n
                   install.packages("BiocManager")
                   BiocManager::install("SeqVarTools")')
      }
      
      if (data.type == "vcf.file") {
        data <- radiator::read_vcf(
          data = data, 
          strata = strata, 
          vcf.stats = FALSE,
          parallel.core = parallel.core, 
          verbose = FALSE, 
          internal = TRUE, 
          path.folder = path.folder)
        data.type <- "SeqVarGDSClass"
        # tidy.data <- gds2tidy(gds = gds, parallel.core = parallel.core, 
      }
      if (data.type == "gds.file") {
        data <- radiator::read_rad(data, verbose = verbose)
        data.type <- "SeqVarGDSClass"
      }
      
      individuals <- radiator::generate_id_stats(
        gds = data, 
        coverage = FALSE,
        path.folder = path.folder, 
        file.date = file.date, 
        parallel.core = parallel.core, 
        verbose = TRUE)
      id.na.stats <- individuals$stats
      individuals <- individuals$info
      strata.df <- dplyr::select(individuals, -c(MISSING_PROP, HETEROZYGOSITY)) %>% 
        dplyr::rename(POP_ID = STRATA)
      
      tidy.data <- suppressWarnings(
        radiator::gds2tidy(
          gds = data,
          calibrate.alleles = FALSE,
          parallel.core = parallel.core
        ) %>% 
          dplyr::select(dplyr::one_of(want)) %>% 
          dplyr::mutate(
            GT_MISSING_BINARY = dplyr::if_else(is.na(GT_BIN), 0L, 1L)
            )
      )
    }
    
  } else {
    want <- unique(c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "POP_ID", "GT", "GT_BIN", strata.select))
    tidy.data <- suppressWarnings(
      radiator::tidy_genomic_data(
        data = data,
        strata = strata,
        filename = filename,
        parallel.core = parallel.core,
        vcf.metadata = FALSE,
        vcf.stats = FALSE,
        verbose = FALSE, 
        internal = TRUE
      ) %>% 
        dplyr::select(dplyr::one_of(want))
    )
    
    strata.df <- suppressWarnings(
      dplyr::ungroup(tidy.data) %>%
        dplyr::select(dplyr::one_of(c("INDIVIDUALS", "POP_ID", strata.select))) %>% 
        dplyr::distinct(INDIVIDUALS, .keep_all = TRUE))
    
    # New column with GT_MISSING_BINARY O for missing 1 for not missing...
    if (tibble::has_name(tidy.data, "GT_BIN")) {
      tidy.data <- dplyr::mutate(
        .data = tidy.data,
        GT_MISSING_BINARY = dplyr::if_else(is.na(GT_BIN), "0", "1"),
        GT_MISSING_BINARY = as.numeric(GT_MISSING_BINARY)
      )
    } else {
      tidy.data <- dplyr::mutate(
        .data = tidy.data,
        GT_MISSING_BINARY = dplyr::if_else(GT == "000000", "0", "1"),
        GT_MISSING_BINARY = as.numeric(GT_MISSING_BINARY)
      )
    }
  }# finish importig data
  
  
  # Check if stata have different values
  check.levels <- function(x) nlevels(factor(x)) > 1
  strata.check <- dplyr::select(strata.df, dplyr::one_of(strata.select)) %>% 
    dplyr::distinct(.) %>% 
    dplyr::summarise_all(.tbl = ., .funs = check.levels) %>% 
    purrr::flatten_lgl(.) %>% 
    unique
  if (length(strata.check) > 1 || !isTRUE(strata.check)) {
    stop("more than 1 value in strata groupings required")
  }
  check.levels <- strata.check <- NULL
  
  
  # some statistics  -----------------------------------------------------------
  message("\nInformations:")
  strata.stats <- strata.df %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::tally(.) %>%
    dplyr::mutate(STRATA = stringi::stri_join(POP_ID, n, sep = " = "))
  
  duplicate.id <- nrow(strata.df) - length(unique(strata.df$INDIVIDUALS))
  
  if (tibble::has_name(tidy.data, "CHROM")) {
    n.chrom = dplyr::n_distinct(tidy.data$CHROM)
    n.locus = dplyr::n_distinct(tidy.data$LOCUS)
  }
  n.snp = dplyr::n_distinct(tidy.data$MARKERS)
  
  # output the proportion of missing genotypes
  na.before <- dplyr::summarise(
    .data = tidy.data,
    MISSING = round(length(GT_MISSING_BINARY[GT_MISSING_BINARY == 0])/length(GT_MISSING_BINARY), 6)) %>%
    purrr::flatten_dbl(.) %>% format(., scientific = FALSE)
  
  marker.problem <- radiator::detect_all_missing(data = tidy.data)
  if (marker.problem$marker.problem) {
    tidy.data <- marker.problem$data
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
  if (tibble::has_name(tidy.data, "CHROM")) {
    message("Number of chrom/scaffolds: ", n.chrom)
    message("Number of locus: ", n.locus)
  }
  message("Number of SNPs: ", n.snp)
  
  message("\nProportion of missing genotypes (overall): ", na.before)
  
  
  # missingness per individuals (required now in the IBM with PCoA) ------------
  res$missing.genotypes.ind <- tidy.data %>%
    dplyr::select(MARKERS, INDIVIDUALS, POP_ID, GT_MISSING_BINARY) %>%
    dplyr::group_by(INDIVIDUALS, POP_ID) %>%
    dplyr::summarise(
      MISSING_GENOTYPE = length(GT_MISSING_BINARY[GT_MISSING_BINARY == 0]),
      MARKER_NUMBER = length(MARKERS),
      MISSING_GENOTYPE_PROP = MISSING_GENOTYPE/MARKER_NUMBER,
      PERCENT = round((MISSING_GENOTYPE_PROP)*100, 2), 
      .groups = "keep"
    ) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(POP_ID, INDIVIDUALS)
  
  strata.missing <- dplyr::inner_join(
    strata.df,
    dplyr::select(res$missing.genotypes.ind, INDIVIDUALS, MISSING_GENOTYPE_PERCENT = PERCENT),
    by = "INDIVIDUALS")
  
  # Identity-by-missingness (IBM) analysis -------------------------------------
  # MultiDimensional Scaling analysis (MDS) - Principal Coordinates Analysis (PCoA)
  message("\n\nIdentity-by-missingness (IBM) analysis using\n    Principal Coordinate Analysis (PCoA)...")
  
  input.pcoa <- tidy.data %>%
    dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT_MISSING_BINARY) %>%
    data.table::as.data.table(.) %>%
    data.table::dcast.data.table(
      data = .,
      formula = POP_ID + INDIVIDUALS ~ MARKERS,
      value.var = "GT_MISSING_BINARY"
    ) %>%
    tibble::as_tibble(.) %>% 
    dplyr::ungroup(.)
  
  # we need rownames for this
  rownames.pcoa <- input.pcoa$INDIVIDUALS
  # suppressWarnings(rownames(input.pcoa) <- input.pcoa$INDIVIDUALS)
  input.pcoa %<>% dplyr::select(-POP_ID, -INDIVIDUALS)
  
  radiator::write_rad(data = input.pcoa, path = file.path(path.folder, "input.rda.temp.rad"))
  input.rda <- list.files(path = path.folder, pattern = "input.rda.temp.rad", full.names = TRUE)
  input.pcoa <- data.frame(input.pcoa)
  rownames(input.pcoa) <- rownames.pcoa
  
  
  #Legendre's pcoa in ape
  safe_pcoa <-  purrr::safely(.f = ape::pcoa)
  # ibm <- ape::pcoa(D = stats::dist(x = input.pcoa, method = distance.method))
  ibm <- safe_pcoa(D = stats::dist(x = input.pcoa, method = distance.method))
  
  if (is.null(ibm$error)) {
    ibm <- ibm$result
  } else {
    rlang::abort("error during the PCoA open an issue on grur's github page")
  }
  
  input.pcoa <- NULL
  
  # Should broken_stick values be reported?
  # variance
  variance.component <-  tibble::tibble(EIGENVALUES = ibm$values$Eigenvalues) %>%
    dplyr::mutate(
      VARIANCE_PROP = round(EIGENVALUES/sum(EIGENVALUES), 2)
    )
  
  res$vectors <- dplyr::inner_join(
    strata.missing,
    tibble::as_tibble(x = data.frame(ibm$vectors), rownames = "INDIVIDUALS")
    , by = "INDIVIDUALS"
  )
  
  ibm <- NULL
  
  # adjust pop_id
  res$vectors <- radiator::change_pop_names(data = res$vectors, pop.levels = NULL)
  
  message("Generating Identity by missingness plot")
  pc.to.do <- utils::combn(1:4, 2, simplify = FALSE)
  
  res$ibm.plots <- purrr::map(
    .x = strata.select, 
    .f = generate_pcoa_plot,
    pc.to.do = pc.to.do,
    vectors = res$vectors,
    variance.component = variance.component,
    path.folder = path.folder, write.plot = write.plot) %>%
    purrr::flatten(.)
  
  variance.component <- pc.to.do <- NULL
  
  # RDA missing data analysis --------------------------------------------------
  message("Redundancy analysis...")
  res$rda.analysis <- missing_rda(data = input.rda, strata = strata.df, 
                                  permutations = 1000, parallel.core = parallel.core)
  
  file.remove(input.rda)
  
  # Eric Archer's code here ----------------------------------------------------
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
    data = tidy.data,
    ci = 0.95,
    path.folder = path.folder,
    write.plot = write.plot) %>%
    purrr::flatten(.)
  
  # Heatmap --------------------------------------------------------------------
  res$heatmap <- tidy.data %>%
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
    ggplot2::labs(y = "Markers", x = "Individuals") +
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
    ggplot2::theme_bw() +
    ggplot2::facet_grid(~POP_ID, scales = "free", space = "free_x")
  # res$heatmap
  
  # Missing summary ------------------------------------------------------------
  message("Generating missing information summary tables and plots")
  
  # Individuals-----------------------------------------------------------------
  message("Missingness per individuals")
  # Figures
  axis.title.element.text.fig <- ggplot2::element_text(
    size = 12, family = "Helvetica", face = "bold")
  axis.text.element.text.fig <- ggplot2::element_text(
    size = 10, family = "Helvetica")
  
  # manhattan and violin plots
  res$missing.genotypes.ind.plots <- suppressMessages(ggplot2::ggplot(
    data = res$missing.genotypes.ind,
    ggplot2::aes(x = POP_ID, y = MISSING_GENOTYPE_PROP, colour = POP_ID)) +
      ggplot2::geom_jitter(alpha = 0.5) +
      ggplot2::geom_violin(trim = TRUE, fill = NA) +
      ggplot2::geom_boxplot(width = 0.1, fill = NA, outlier.colour = NA, outlier.fill = NA) +
      ggplot2::labs(x = "Populations", 
                    y = "Individual's missing genotypes (proportion)",
                    colour = "Populations") +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        axis.title.x = axis.title.element.text.fig,
        axis.text.x = axis.text.element.text.fig,
        axis.title.y = axis.title.element.text.fig,
        axis.text.y = axis.text.element.text.fig
      ) +
      ggplot2::theme_bw() +
      ggplot2::coord_flip())
  # res$missing.genotypes.ind.plots
  
  # histogram
  res$missing.genotypes.ind.histo <- suppressMessages(ggplot2::ggplot(
    data = res$missing.genotypes.ind,
    ggplot2::aes(x = MISSING_GENOTYPE_PROP)) +
      ggplot2::geom_histogram() +
      ggplot2::labs(x = "Individual's missing genotypes (proportion)",
                    y = "Individuals (number)") +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = axis.title.element.text.fig,
        axis.title.y = axis.title.element.text.fig,
        axis.text.x = axis.text.element.text.fig,
        axis.text.y = axis.text.element.text.fig
      ) +
      ggplot2::theme_bw()
  )
  # res$missing.genotypes.ind.histo
  
  # helper plot for individual's genotyped threshold
  res$ind.genotyped.helper.plot <- ind_genotyped_helper(res$missing.genotypes.ind)
  
  
  # merge plots
  plots <- cowplot::align_plots(
    res$missing.genotypes.ind.plots, res$ind.genotyped.helper.plot, align = "v", axis = "l")
  top.row.plot <- suppressMessages(cowplot::plot_grid(
    plots[[1]], res$missing.genotypes.ind.histo, labels = c("A", "B"), align = "h"))
  
  res$missing.genotypes.ind.combined.plots <- suppressMessages(cowplot::plot_grid(
    top.row.plot, plots[[2]],
    ncol = 1, nrow = 2, labels = c("", "C"), rel_heights = c(1.5, 1)))
  plots <- top.row.plot <- NULL
  
  if (write.plot) {
    
    # cowplot author says it's better to use his save_plot... will try
    cowplot::save_plot(
      filename = stringi::stri_join(path.folder, "/missing.genotypes.ind.combined.plots.pdf"),
      plot = res$missing.genotypes.ind.combined.plots, 
      base_height = n.pop * 1,
      base_aspect_ratio = 2.5, ncol = 2, nrow = 2, limitsize = FALSE)
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
    tidyr::pivot_longer(
      data = .,
      cols = tidyselect::everything(),
      names_to = "BLACKLIST",
      values_to = "n"
    ) %>%
    dplyr::transmute(BLACKLIST = stringi::stri_join(BLACKLIST, n, sep = " = "))
  message("    Number of individual(s) blacklisted per blacklist generated:\n", stringi::stri_join("    ", blacklists.stats$BLACKLIST, collapse = "\n"))
  res <- c(res, blacklists)
  blacklists.stats <- blacklists <- NULL
  
  # FH -------------------------------------------------------------------------
  message("Calculation of FH: a measure of IBDg")
  fh <- radiator::ibdg_fh(data = tidy.data, path.folder = path.folder, verbose = FALSE)
  res$missing.genotypes.ind.fh <- suppressWarnings(
    dplyr::full_join(
      res$missing.genotypes.ind,
      fh$fh
      , by = c("INDIVIDUALS", "POP_ID")
    )
  )
  # res$fh.manhattan.box.plot <- fh$fh.manhattan.box.plot
  
  res$missing.genotypes.ind.fh.plots <- ggplot2::ggplot(
    res$missing.genotypes.ind.fh, ggplot2::aes(y = FH, x = MISSING_GENOTYPE_PROP)) +
    ggplot2::geom_point() +
    ggplot2::stat_smooth(method = stats::lm, level = 0.99) +
    # labs(title = "Correlation between missingness and inbreeding coefficient") +
    ggplot2::labs(y = "Individual IBDg (FH)") +
    ggplot2::labs(x = "Missing genotype (proportion)") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    ) +
    ggplot2::theme_bw()
  
  # merge plots
  plots <- cowplot::align_plots(
    res$missing.genotypes.ind.fh.plots, fh$fh.box.plot, align = "v", axis = "l")
  bottom.row.plot <- suppressMessages(cowplot::plot_grid(
    plots[[2]], fh$fh.distribution.plot, labels = c("B", "C"), align = "h"))
  
  res$missing.genotypes.ind.fh.combined.plots <- cowplot::plot_grid(
    plots[[1]], bottom.row.plot,
    ncol = 1, nrow = 2, labels = c("A", ""), rel_heights = c(1.5, 1))
  # res$missing.genotypes.ind.fh.plots
  fh <- plots <- bottom.row.plot <- NULL
  
  if (write.plot) {
    cowplot::save_plot(
      filename = stringi::stri_join(path.folder, "/missing.genotypes.ind.fh.combined.plots.pdf"),
      plot = res$missing.genotypes.ind.fh.combined.plots, 
      base_height = n.pop * 1,
      base_aspect_ratio = 1, ncol = 2, nrow = 2, limitsize = FALSE)
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
  
  want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "GT_MISSING_BINARY")
  tidy.col <- colnames(tidy.data)
  markers.meta <- purrr::keep(
    .x = tidy.col,
    .p = tidy.col %in% c("MARKERS", "CHROM", "LOCUS", "POS"))
  res$missing.genotypes.markers.overall <- suppressWarnings(
    tidy.data %>%
      dplyr::select(dplyr::one_of(want)) %>% 
      dplyr::group_by_if(.tbl = ., .predicate = colnames(x = .) %in% markers.meta) %>% 
      dplyr::summarise(
        MISSING_GENOTYPE = length(GT_MISSING_BINARY[GT_MISSING_BINARY == 0]),
        INDIVIDUALS_NUMBER = length(INDIVIDUALS),
        MISSING_GENOTYPE_PROP = MISSING_GENOTYPE / INDIVIDUALS_NUMBER,
        PERCENT = round(MISSING_GENOTYPE_PROP * 100, 2)
      ) %>%
      dplyr::ungroup(.) %>%
      dplyr::arrange(MARKERS))
  
  res$missing.genotypes.markers.pop <- tidy.data %>%
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
  
  # # Save tidy Note to myself: think it's saved when people use filename
  # saved via tidy_genomic_data.
  # fst::write.fst(x = tidy.data, path = file.path(path.folder, "tidy.data.rad"))
  tidy.data <- NULL
  
  markers.missing.geno.threshold <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 ,0.8, 0.9)
  whitelists <- purrr::map(
    .x = markers.missing.geno.threshold,
    .f = whitelists_markers_generator,
    y = res$missing.genotypes.markers.overall,
    path.folder = path.folder) %>%
    purrr::flatten(.)
  
  message("Whitelist(s) of markers generated: ", length(whitelists))
  if (length(whitelists) > 0) {
    whitelists.stats <- purrr::map_df(.x = whitelists, .f = nrow) %>% 
      tidyr::pivot_longer(
        data = .,
        cols = tidyselect::everything(),
        names_to = "WHITELIST",
        values_to = "n"
      ) %>%
      dplyr::transmute(WHITELIST = stringi::stri_join(WHITELIST, n, sep = " = "))
    message("    Number of markers whitelisted per whitelist generated:\n", stringi::stri_join("    ", whitelists.stats$WHITELIST, collapse = "\n"))
  }
  # res <- c(res, whitelists)
  whitelists.stats <- whitelists <- NULL
  
  # Figure markers
  
  # violin plots
  res$missing.genotypes.markers.plots <- suppressMessages(ggplot2::ggplot(
    data = res$missing.genotypes.markers.pop,
    ggplot2::aes(x = POP_ID, y = MISSING_GENOTYPE_PROP, colour = POP_ID)) +
      # ggplot2::geom_jitter(alpha = 0.5) +
      ggplot2::geom_violin(trim = TRUE, fill = NA) +
      ggplot2::geom_boxplot(width = 0.1, fill = NA, outlier.colour = NA, outlier.fill = NA) +
      ggplot2::labs(y = "Marker's missing genotypes (proportion)", x = "Populations",
                    colour = "Populations") +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 10, family = "Helvetica"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.y = ggplot2::element_text(size = 10, family = "Helvetica")
      ) +
      ggplot2::theme_bw() +
      ggplot2::coord_flip())
  # res$missing.genotypes.markers.plots
  
  res$markers.histo <- suppressMessages(ggplot2::ggplot(
    data = res$missing.genotypes.markers.overall, ggplot2::aes(x = MISSING_GENOTYPE_PROP)) +
      ggplot2::geom_histogram() +
      ggplot2::labs(x = "Marker's missing genotypes (proportion)") +
      ggplot2::labs(y = "Markers (number)") +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.text.x = ggplot2::element_text(size = 10, family = "Helvetica"),
        strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
      ) +
      ggplot2::theme_bw()
  )
  # res$markers.histo
  
  # helper plot for markers's genotyped threshold
  res$markers.helper.plot <- markers_genotyped_helper(
    x = res$missing.genotypes.markers.pop, y = res$missing.genotypes.markers.overall)
  
  # merge plots
  plots <- cowplot::align_plots(
    res$missing.genotypes.markers.plots, res$markers.helper.plot, align = "v", axis = "l")
  top.row.plot <- suppressMessages(cowplot::plot_grid(
    plots[[1]], res$markers.histo, labels = c("A", "B"), align = "h"))
  
  res$missing.genotypes.markers.combined.plots <- suppressMessages(
    cowplot::plot_grid(
      top.row.plot, plots[[2]],
      ncol = 1, nrow = 2, labels = c("", "C"), rel_heights = c(1.5, 1)))
  # res$missing.genotypes.markers.plots
  plots <- top.row.plot <- NULL
  
  if (write.plot) {
    cowplot::save_plot(
      filename = stringi::stri_join(
        path.folder, "/missing.genotypes.markers.combined.plots.pdf"),
      plot = res$missing.genotypes.markers.combined.plots, 
      base_height = n.pop * 1,
      base_aspect_ratio = 2.5, ncol = 2, nrow = 2, limitsize = FALSE)
  }
  
  # Results --------------------------------------------------------------------
  return(res)
}
# Internal nested functions ----------------------------------------------------
# Are now moved in file: internal.R
