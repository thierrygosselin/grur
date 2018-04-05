# grur imputations module

#' @name grur_imputations
#' @title Map-independent imputations of missing genotypes
#'
#' @description The goal of this module is to provide a simple solution for
#' a complicated problem: missing genotypes in RADseq genomic datasets.
#' This function will performed \strong{map-independent imputations} of missing
#' genotypes.
#'
#' \strong{Key features:}
#'
#' \itemize{
#' \item \strong{Imputation methods: } several machine learning algorithms;
#' Random forests (predictive and on-the-fly-imputations); 
#' Extreme gradient tree boosting (XGBoost);
#' Fast, distributed, high performance gradient
#' boosting (LightGBM) and
#' Multiple Correspondence Analysis (MCA). Furthermore, the module allows to
#' compare these algorithms with the classic Strawman imputation
#' (~ max/mean/mode: the most frequently observed, i.e. non-missing, genotypes is used).
#' \item \strong{Hierarchical level: } Imputations conducted by strata or globally.
#' \item \strong{SNP linkage/Haplotypes: } Correlation among SNPs is accounted for during
#' rf and tree boosting imputations, i.e. the imputations are 
#' automatically conducted by haplotypes
#' when marker meta-information is available (chromosome, locus and position,
#' usually taken from VCF files). The alternative, considers all the markers 
#' independent and imputation is conducted by SNPs/markers.
#' \item \strong{Genotype likelihood (GL): } The GL info is detected automatically
#' (GL column in FORMAT field of VCF files). Genotypes with higher likelihood
#' will have higher probability during bootstrap samples of trees in Random
#' forests (under devel).
#' \item \strong{Predictive mean matching: } random forest in prediction mode
#' uses a fast k-nearest neighbor (KNN) searching algorithms 
#' (see argument documentation and details below).
#' \item \strong{Optimized for speed: } There's a
#' \href{https://github.com/thierrygosselin/grur#troubleshooting}{tutorial} and a
#' \href{https://github.com/thierrygosselin/grur#vignettes-and-examples}{vignette}
#' detailing the procedure to install the packages from source to enable OpenMP.
#' Depending on algorithm used, a progress bar is sometimes available to see 
#' if you have time for a coffee break!
#' }
#'
#'
#' Before running this function to populate the original dataset with synthetic
#' data, I highly recommend you look for patterns of missingness
#' \code{\link[grur]{missing_visualization}}
#' and explore the reasons for their presence 
#' (see \href{https://www.dropbox.com/s/4zf032g6yjatj0a/vignette_missing_data_analysis.nb.html?dl=0}{vignette}).


#' @param data A tidy genomic dataset.
#' It can be file in the working directory or
#' an object in the global environment.
#' To get a tidy dataset from various genomic format, see
#' \href{https://github.com/thierrygosselin/radiator}{radiator}
#' \code{\link[radiator]{tidy_genomic_data}}.
#' \emph{See details of this function for more info}.

#' @param strata (optional) The strata file is a tab delimited file with
#' 2 columns headers: \code{INDIVIDUALS} and \code{STRATA}. 
#' If a \code{strata} file is specified, the strata argument will have precedence
#' on the population groupings used internally in the dataset (e.g. \code{POP_ID}.
#' The \code{STRATA} column can be any hierarchical grouping.
#' To create a strata file see \code{\link[radiator]{individuals2strata}}.
#' If you have missing data patterns associated with sequencing lanes, use the 
#' different lanes id associated with your samples and store the info in the
#' \code{STRATA} column.
#' Default: \code{strata = NULL}.

#' @param subsample.markers (optional, integer) To speed up computation and rapidly
#' test the function's arguments (e.g. using 200 markers).
#' Default: \code{subsample.markers = NULL}.

#' @param imputation.method (character, optional)
#' Methods available for map-independent imputations of missing genotypes
#' (see details for more info):
#'
#' \enumerate{
#' \item \code{imputation.method = "max"} Strawman imputation,
#' the most frequently observed genotypes (ties are broken at random).
#'
#' \item \code{imputation.method = "rf"} On-the-fly-imputations using
#' Random Forests algorithm.
#'
#' \item \code{imputation.method = "rf_pred"} Random Forests algorithm is used
#' as a prediction problem.
#'
#' \item \code{imputation.method = "xgboost"} extreme gradient boosting trees
#' using depth-wise tree growth.
#' 
#' \item \code{imputation.method = "lightgbm"} for a light and fast 
#' leaf-wise tree growth gradient boosting algorithm (in devel).
#'
#' \item \code{imputation.method = "bpca"} Multiple Correspondence Analysis (in devel).
#'
#' \code{imputation.method = NULL} the function will stop.
#' Default: \code{imputation.method = NULL}.
#' }

#' @param hierarchical.levels (character, optional) \code{c("global", "strata")}.
#' Should the imputations be computed by markers globally or by strata.
#' Historically, this was \code{"populations"}.
#'
#' Note that imputing genotype globally in conjunction with
#' \code{imputation.method = "max" or "strata"} can potentially create huge bias.
#' e.g. by introducing foreign genotypes/haplotypes in some strata/populations
#' (see note for more info).
#' Default: \code{hierarchical.levels = "strata"}.


#' @param pmm (integer, optional) Predictive mean matching used in conjunction with
#' random Forests and lightgbm (\code{imputation.method = "rf_pred"} or
#' \code{imputation.method = "lightgbm"}).
#' Number of candidate non-missing
#' value to sample from during the predictive mean matching step.
#' A fast k-nearest neighbor searching algorithms is used with this approach.
#' \code{pmm = 2} will use 2 neighbors.
#' Default: \code{pmm = 0}, will avoids this step.

#' @param random.seed (integer, optional) For reproducibility, set an integer
#' that will be used to initialize the random generator. With default,
#' a random number is generated. Currently not implemented for XGBoost and LightGBM.
#' Default: \code{random.seed = NULL}.

#' @param verbose (optional, logical) When \code{verbose = TRUE}
#' the function is a little more chatty during execution.
#' Default: \code{verbose = TRUE}.

#' @param parallel.core (optional, integer) The number of core used for parallel.
#' For some algorithms, markers will be imputed in parallel and 
#' strata are processed sequentially.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @param cpu.boost (optional, integer). Number of core for XGBoost and LightGBM.
#' For the best speed, set this to the number of real CPU cores,
#' not the number of threads. 
#' Most CPU using hyper-threading to generate 2 threads per CPU core.
#' Be aware that task manager or any CPU monitoring tool might report cores
#' not being fully utilized. This is normal.
#' Default: \code{cpu_boost = parallel::detectCores() / 2}.


#' @param filename (optional) The function uses \code{\link[fst]{write.fst}},
#' to write the tidy data frame in
#' the working directory. The file extension appended to
#' the \code{filename} provided is \code{.rad}.
#' With default: \code{filename = NULL}, the imputed tidy data frame is
#' in the global environment only (i.e. not written in the working directory...).

#' @param ... (optional) To pass further argument for fine-tuning your
#' imputations. See details below.

#' @return The output in your global environment is the imputed tidy data frame.
#' If \code{filename} is provided, the imputed tidy data frame is also
#' written to the working directory.

#' @details
#' \strong{Predictive mean matching:}
#'
#' Random Forests already behave like a nearest neighbor
#' classifier, with adaptive metric. Now we have the option to conduct
#' predictive mean matching on top of the prediction based missing value
#' imputation. PMM tries to raise the variance in the resulting conditional
#' distributions to a realistic level.
#' The closest k predicted values are identified by a fast
#' k-nearest neighbour approach wrapped in the package
#' \href{https://github.com/mayer79/missRanger}{missRanger}
#' Returned value correspond to the mean value.
#'
#'
#' \strong{haplotype/SNP approach:}
#'
#' The \emph{haplotype approach} is automatically used when markers meta-information
#' is detected (chromosome/CHROM, locus/ID and SNP/POS columns, usually from a VCF file).
#' Missing genotypes from SNPs on the same locus or same RADseq read is undertaken
#' simulteneously to account for the correlation of the linked SNPs. When one or
#' more SNPs on the same read/haplotype is missing, the haplotype is deleted and
#' consequently, imputation might results in different genotype for those SNPs
#' that were not missing. This approach is much safer than potentially creating
#' weird chimeras during haplotype imputations.
#' Alternatively, a \emph{snp approach} is used, and the SNP are considered
#' independent. Imputations of genotypes is then conducted for each marker separately.
#'
#'
#' \strong{Imputing globally or by strata ?}
#' \code{hierarchical.levels = "global"} argument will act differently depending
#' on the \code{imputation.method} selected.
#'
#' \strong{Strawman imputations (~ max/mean/mode) considerations: }
#'
#' With \code{imputation.method = "max"} and \code{hierarchical.levels = "global"}
#' \emph{will likely create bias}.
#'
#' \emph{Example 1 (unbalanced sample size):} Consider 2 populations evolving more
#' by drift than selection: pop1 (n = 36) and pop2 (n = 50).
#' You'll likely have a few polymorphic marker, where pop1 and pop2 are
#' monomorphic for different alleles (pop1 is fixed for the minor/ALT allele and
#' pop2 is fixed for the major/REF allele). Missing genotypes in pop1
#' using the most common filling technique in the literature (using mean/mode/max),
#' will result in pop1 having individuals with the REF allele.
#' Not something you want... unless your population membership is not 100% accurate,
#' (e.g. you might have migrants or wrong assignation),
#' which in this case you still don't want to impute with
#' \code{imputation.method = "max"} (see alternative below).
#'
#' \emph{Example 2 (balanced sample size):} pop1 (n = 100) and pop2 (n = 100).
#' For a particular marker, pop1 as 85 individuals genotyped and pop2 100.
#' Again, if the populations are fixed for different alleles
#' (pop1 = ALT and pop2 = REF), you will end up having REF allele in your pop1,
#' not something you want... unless your population membership is not 100% accurate,
#' (e.g. you might have migrants or wrong assignation),
#' which in this case you still don't want to impute with
#' \code{imputation.method = "max"} (see alternative below).
#'
#' \strong{Random Forests imputations: }
#'
#' Random Forests use machine learning and you can take this into account while
#' choosing argument values. Uncertain of the groupings ? Use random forests with
#' \code{hierarchical.levels = "global"}. Random forests will account for the
#' potential linkage and correlation between
#' markers and genotypes to make the best imputations available. This can potentially
#' results in genotypes for a certain combo population/marker with new groupings
#' (e.g. a new allele). This is much more accurate and not the same thing as
#' the \code{imputation.method = "max"} because the imputed genotype is validated
#' by considering all other variables (all the other markers genotyped for the individual).
#' \emph{Test the option and report bug if you find one.}
#'
#' \emph{random forest with on-the-fly-imputation (rf): }the technique is described
#' in Tang and Ishwaran (2017). Non-missing genotypes are used for
#' the split-statistics. Daughter node assignation membership use random
#' non-missing genotypes from the inbag data. Missing genotypes are imputed at
#' terminal nodes using maximal class rule with out-of-bag non-missing genotypes.
#' Example of computation time: for 1500 individuals, 20 000 markers and
#' 2% of overall missing genotypes using the hierarchical levels globally
#' takes around 10 hours on 6 CPUs to complete the imputations.
#'
#' \emph{random forest as a prediction problem (rf_pred): }markers with
#' missing genotypes are imputed one at a time. The fitted forest is used to
#' predict missing genotypes. Missingness in the response variables are
#' incorporated as attributes for growing the forest.
#' 
#' \emph{boosting trees:} Prediction method is used for both XGBoost and LightGBM.
#' 
#' \href{https://lightgbm.readthedocs.io/en/latest/index.html}{LightGBM documentation}
#' 
#' \href{http://xgboost.readthedocs.io/en/latest/}{XGBoost documentation}
#' 
#' \href{https://sites.google.com/view/lauraepp/parameters}{LightGBM and XGBoost parameters optimization}
#' 
#' 
#' 
#'
#' \strong{... :dot dot dot arguments}
#'
#' The argument is available to tailor your imputations using
#' XGBoost, LightGBM and Random Forests:
#'
#' Available arguments for eXtreme Gradient Boosting tree method:
#' \emph{eta, gamma, max_depth, min_child_weight, subsample, colsample_bytree,
#' num_parallel_tree, nrounds, save_name, early_stopping_rounds}.
#' Refer to \code{\link[xgboost]{xgboost}} for arguments documentation.
#'
#' Available arguments for LightGBM:
#' \emph{boosting, objective, learning_rate, feature_fraction, bagging_fraction,
#' bagging_freq, max_depth, min_data_in_leaf, num_leaves, early_stopping_rounds,
#' nrounds, max_depth, iteration.subsample}. Refer to \code{\link[lightgbm]{lightgbm}}
#' for arguments documentation. \code{iteration.subsample} is the number of iteration
#' to conduct training of the model with new subsamples (default: \code{iteration.subsample = 2}).
#' 
#' Available arguments for Random forests method:
#' \emph{num.tree, nodesize, nsplit, nimpute}.
#' Refer to \code{\link[randomForestSRC]{impute.rfsrc}} for arguments documentation.
#'
#'
#  Multiple Correspondence Analysis option available (upcomming):
# \emph{ncp}.
# Refer to \code{\link[missMDA]{imputeMCA}} for argument documentation.

#' @note
#'
#' \strong{Reference genome or linkage map available ?}
#' Numerous approaches are available and more appropriate, please search
#' the literature
#' (\href{https://online.papersapp.com/collections/05d6e65a-73c9-49e6-9c75-289a818f76f3/share}{references}).
#'
#'
#' \strong{What's the simple imputation message when running the function ?}
#'
#' Before conducting the imputations by strata with the machine learning algorithms, 
#' the data is first screened for markers that are
#' monomorphic WITHIN the strata/populations.
#' Because for those cases, it's clear what the
#' missing genotypes should be, the imputations is very \emph{simple} and missing
#' genotypes are imputed with the only genotype found for the particular strata/populations.
#' There's no need for a fancy method here.
#' The small cost in time is worth it, because the model inside 
#' the machine learning algorithms will benefit from having a more complete and
#' reliable genotype matrix.
#'
#'
#' \strong{Deprecated arguments:}
#'
#' \itemize{
#' \item \code{hierarchical.levels = "populations"} update your codes with "strata".
#' \item \code{imputations.group} is now replaced by \code{hierarchical.levels}
#' \item \code{impute} is no longer available.
#' Imputing using \code{impute = "allele"} option was wrong because it
#' was using F1 genotypes for imputations. Now imputation is only conducted at
#' the genotype level.
#' \item \code{iteration.rf} is no longer used. This argument is now available
#' inside the \code{...} for on-the-fly-imputations (see details). The default
#' is now set to 10 iterations.
#' \item \code{split.number} is automatically set.
#' }

#' @seealso
#' The package \href{https://github.com/imbs-hl/ranger}{ranger}
#' (see Wright and Ziegler, 2016) provides a fast C++ version
#' of the original implementation of rf from Breiman (2001).
#' 
#' \href{https://github.com/mayer79/missRanger}{missRanger}
#' 
#' The package \href{https://github.com/dmlc/xgboost}{randomForestSRC} 
#' (see Tang and Ishwaran, 2017) provides
#' a fast on-the-fly imputation method.
#' 
#' \href{https://github.com/stekhoven/missForest}{missForest}
#' 
#' The package \href{https://github.com/dmlc/xgboost}{XGBoost} 
#' (Chen and Guestrin, 2016) provides
#' a fast C++ implementation for the extreme gradient tree boosting algorithm.
#' 
#' The package \href{https://github.com/Microsoft/LightGBM}{LightGBM} is an
#' exciting new algorithm to conduct tree boosting.




#' @export
#' @rdname grur_imputations
#' @importFrom parallel detectCores
#' @importFrom dplyr distinct group_by ungroup rename arrange tally filter select select_ one_of mutate mutate_all summarise left_join funs bind_rows
#' @importFrom tidyr gather unite drop_na
#' @importFrom purrr map flatten keep discard flatten_chr flatten_dbl flatten_lgl safely
#' @importFrom purrrlyr invoke_rows
#' @importFrom stringi stri_replace_na
#' @importFrom tibble has_name as_data_frame
#' @importFrom stats predict reformulate as.formula
#' @importFrom rlang .data
# @importFrom base split
#' @importFrom readr write_lines write_tsv
#' @importFrom radiator tidy_wide change_alleles detect_biallelic_markers
# @importFrom fst write.fst
#' @importFrom Matrix Matrix
# @importFrom lightgbm lgb.Dataset lgb.train
# @importFrom randomForestSRC impute.rfsrc
# @importFrom ranger ranger
# @importFrom missRanger pmm
# @importFrom xgboost xgb.DMatrix cb.early.stop xgb.train

#' @examples
#' \dontrun{
#' # The simplest way to run when you have a tidy dataset:
#'
#' wolf.imputed <- grur::grur_imputations(
#' data = "wolf.tidy.dataset.rad",
#' imputation.method = "lightgbm")
#'
#' # This will impute the missing genotypes by strata (in this case the strata,
#' is the population), population using lightgbm.
#' # The remaining arguments will be the defaults.
#'
#' # When you start with a vcf file you can use magrittr %>% to `pipe` the
#' # result. Below, an example with more arguments offered by the functions:
#'
#' wolf.imp <- radiator::tidy_genomic_data(
#'     data = "batch_1.vcf",
#'     strata = "strata.wolf.10pop.tsv",
#'     vcf.metadata = TRUE,
#'     whitelist.markers = "whitelist.loci.txt",
#'     verbose = TRUE) %>%
#' grur::grur_imputations(
#'     data = ., imputation.method = "xgboost", parallel.core = 32)
#' }

#' @references Wright, M. N. & Ziegler, A. (2016).
#' ranger: A Fast Implementation of Random Forests for High Dimensional Data
#' in C++ and R.
#' Journal of Statistical Software, in press. http://arxiv.org/abs/1508.04409.
#' @references Breiman, L. (2001). Random forests. Machine learning, 45(1), 5-32.
#' @references Chen T, Guestrin C. (2016).
#' XGBoost: A scalable tree boosting system. arXivorg. 2016.
#' doi:10.1145/2939672.2939785
#' @references Tang F, Ishwaran H. (2017) Random Forest Missing Data Algorithms.
#' arXiorg: 1–24.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

grur_imputations <- function(
  data,
  strata = NULL,
  subsample.markers = NULL,
  imputation.method = NULL,
  hierarchical.levels = "strata",
  pmm = 0,
  verbose = TRUE,
  parallel.core = parallel::detectCores() - 1,
  cpu.boost = parallel::detectCores() / 2,
  random.seed = NULL,
  filename = NULL,
  ...
) {
  opt.change <- getOption("width")
  options(width = 70)
  if (verbose) {
    cat("\n\n")
    cat("#######################################################################\n")
    cat("####################### grur::grur_imputations ########################\n")
    cat("#######################################################################\n")
  }
  timing <- proc.time() #for timing
  
  # Checking for missing and/or default arguments & package required -----------
  if (missing(data)) stop("data file missing")
  
  if (imputation.method == "lightgbm") {
    if (!requireNamespace("lightgbm", quietly = TRUE)) {
      stop("lightgbm needed for this imputation option to work.
Please follow the vignette for install instructions", call. = FALSE)
    }
  }
  
  if (imputation.method == "xgboost") {
    if (!requireNamespace("xgboost", quietly = TRUE)) {
      stop("xgboost needed for this imputation option to work.
Please follow the vignette for install instructions", call. = FALSE)
    }
  }
  
  if (imputation.method == "rf") {
    if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
      stop("randomForestSRC needed for this imputation option to work.
Please follow the vignette for install instructions", call. = FALSE)
    }
  }
  
  if (imputation.method == "rf_pred") {
    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("ranger needed for this imputation option to work.
           Please follow the vignette for install instructions", call. = FALSE)
    }
  }
  
  if (imputation.method == "bpca") stop("Imputations with bpca is under heavy testing... try again soon!")
  # if (imputation.method == "bpca") {
  #   if (!requireNamespace("pcaMethods", quietly = TRUE)) {
  #     stop("pcaMethods needed for this imputation option to work.
  #          Please follow the vignette for install instructions", call. = FALSE)
  #   }
  # }
  
  if (pmm > 0) {
    if (!requireNamespace("missRanger", quietly = TRUE)) {
      stop("missRanger needed for this imputation option to work.
           Please follow the vignette for install instructions", call. = FALSE)
    }
  }
  
  
  # Options selected & Messages ------------------------------------------------
  if (is.null(imputation.method)) {
    message("Imputation method: NULL")
    message("Returning the tidy dataset")
    data.imp <- NULL
  } else {
    message("Imputation method: ", imputation.method)
    message("Hierarchical levels: ", hierarchical.levels)
    # message("Markers linkage: ", markers.linkage, "\n")
    
    # Capture unevaluated ...
    # Inspired by Eric Anderson and Hadley Wickham codes
    # could also importFrom pryr::named_dots
    
    dotslist <- list(...)
    unknowned_param <- setdiff(
      names(dotslist),
      c("eta", "gamma", "max_depth", "min_child_weight", "subsample", "colsample_bytree",
        "num_parallel_tree", "nrounds", "save_name", "early_stopping_rounds",
        "num.tree", "nodesize", "nsplit", "nimpute", "ncp", "boosting",
        "objective", "learning_rate",
        "feature_fraction", "bagging_fraction", "bagging_freq",
        "min_data_in_leaf", "num_leaves",
        "early_stopping_rounds", "iteration.subsample"))
    
    if (length(unknowned_param) > 0) {
      stop("Unknowned \"...\" parameters to grur imputation module: ",
           stringi::stri_join(unknowned_param, collapse = " "))
    }
    
    boost.dots <- dotslist[names(dotslist) %in%
                             c("eta", "gamma", "max_depth", "min_child_weight",
                               "subsample", "colsample_bytree",
                               "num_parallel_tree", "nrounds", "save_name",
                               "early_stopping_rounds", 
                               "num.tree", "nodesize", "nsplit",
                               "nimpute", "ncp",
                               "boosting", "objective", "learning_rate",
                               "feature_fraction", "bagging_fraction", "bagging_freq",
                               "min_data_in_leaf", "num_leaves",
                               "early_stopping_rounds", "iteration.subsample"
                             )]
    
    # randomForestSRC arguments ------------------------------------------------
    
    # @param num.tree (integer, optional) The number of trees to grow
    # when \code{imputation.method = "rf"} or \code{imputation.method = "rf_pred"}.
    # Default: \code{num.tree = 50}.
    
    if (!is.null(boost.dots[["num.tree"]])) {
      num.tree <- boost.dots[["num.tree"]]
    } else {
      num.tree = 50
    }
    
    if (!is.null(boost.dots[["nodesize"]])) {
      nodesize <- boost.dots[["nodesize"]]
    } else {
      nodesize = 1
    }
    if (!is.null(boost.dots[["nsplit"]])) {
      nsplit <- boost.dots[["nsplit"]]
    } else {
      nsplit = 10
    }
    if (!is.null(boost.dots[["nimpute"]])) {
      nimpute <- boost.dots[["nimpute"]]
    } else {
      nimpute = 10
    }
    
    
    # # MCA ----------------------------------------------------------------------
    # 
    # if (!is.null(boost.dots[["ncp"]])) {
    #   ncp <- boost.dots[["ncp"]]
    # } else {
    #   ncp = 2
    # }
    
    # XGBoost arguments --------------------------------------------------------
    # learning rate
    if (!is.null(boost.dots[["eta"]])) {
      eta <- boost.dots[["eta"]]
    } else {
      eta = 0.1
    }
    
    # gamma for minimum loss reduction (larger the number, more conservative the algorithm)
    # prevent overfitting through regularization
    if (!is.null(boost.dots[["gamma"]])) {
      gamma <- boost.dots[["gamma"]]
    } else {
      gamma = 0
    }
    
    # min_child_weight: minimum sum of instance weight (hessian) needed in a child.
    if (!is.null(boost.dots[["min_child_weight"]])) {
      min_child_weight <- boost.dots[["min_child_weight"]]
    } else {
      min_child_weight = 1
    }
    
    # subsample: subsample ratio of the training instance for growing tree and prevent overfitting
    if (!is.null(boost.dots[["subsample"]])) {
      subsample <- boost.dots[["subsample"]]
    } else {
      subsample = 0.5
    }
    
    # colsample_bytree: subsample ratio of columns when growing each tree
    if (!is.null(boost.dots[["colsample_bytree"]])) {
      colsample_bytree <- boost.dots[["colsample_bytree"]]
    } else {
      colsample_bytree = 1
    }
    
    # num_parallel_tree: Experimental parameter.
    # number of trees to grow per round.
    if (!is.null(boost.dots[["num_parallel_tree"]])) {
      num_parallel_tree <- boost.dots[["num_parallel_tree"]]
    } else {
      num_parallel_tree = 1
    }
    
    
    # save_name: the name or path for periodically saved model file
    if (!is.null(boost.dots[["save_name"]])) {
      save_name <- boost.dots[["save_name"]]
    } else {
      save_name = "imputation.model.temp"
    }
    
    # Common to XGBoost and LightGBM arguments ---------------------------------
    
    
    # early_stopping_rounds: If NULL, the early stopping function is not triggered.
    # If set to an integer k, training with a validation set will stop
    # if the performance doesn't improve for k rounds.
    if (!is.null(boost.dots[["early_stopping_rounds"]])) {
      early_stopping_rounds <- boost.dots[["early_stopping_rounds"]]
    } else {
      early_stopping_rounds = 20
    }
    
    # nrounds: the max number of iterations
    if (!is.null(boost.dots[["nrounds"]])) {
      nrounds <- boost.dots[["nrounds"]]
    } else {
      nrounds = 1000
    }
    
    # max_depth
    if (!is.null(boost.dots[["max_depth"]])) {
      max_depth <- boost.dots[["max_depth"]]
    } else {
      if (imputation.method == "xgboost") max_depth = 6
      if (imputation.method == "lightgbm") max_depth = -1# infinite depth
    }
    
    
    # LightGBM arguments -------------------------------------------------------
    
    
    if (!is.null(boost.dots[["boosting"]])) {
      boosting <- boost.dots[["boosting"]]
    } else {
      boosting = "dart"
    }
    
    if (!is.null(boost.dots[["objective"]])) {
      objective <- boost.dots[["objective"]]
    } else {
      objective = "multiclass"
    }
    
    if (!is.null(boost.dots[["learning_rate"]])) {
      learning_rate <- boost.dots[["learning_rate"]]
    } else {
      learning_rate = 0.1
    }
    
    if (!is.null(boost.dots[["feature_fraction"]])) {
      feature_fraction <- boost.dots[["feature_fraction"]]
    } else {
      feature_fraction = 0.9
    }
    
    if (!is.null(boost.dots[["bagging_fraction"]])) {
      bagging_fraction <- boost.dots[["bagging_fraction"]]
    } else {
      bagging_fraction = 0.9
    }
    if (!is.null(boost.dots[["bagging_freq"]])) {
      bagging_freq <- boost.dots[["bagging_freq"]]
    } else {
      bagging_freq = 1
    }
    
    # minimum dat per leaf
    # It is a LightGBM specific feature allowing to choose how many observations
    # must exist before putting them in a leaf during each tree split consideration.
    
    if (!is.null(boost.dots[["min_data_in_leaf"]])) {
      min_data_in_leaf <- boost.dots[["min_data_in_leaf"]]
    } else {
      min_data_in_leaf = 20 # dataset with 100 obs, needs to lower this (in XGBoost = 1)
    }
    # number of leaves in one tree
    if (!is.null(boost.dots[["num_leaves"]])) {
      num_leaves <- boost.dots[["num_leaves"]]
    } else {
      num_leaves = 31
    }
    
    # iteration.subsample
    # number of leaves in one tree
    if (!is.null(boost.dots[["iteration.subsample"]])) {
      iteration.subsample <- boost.dots[["iteration.subsample"]]
    } else {
      iteration.subsample = 2
    }
    
    if (verbose) {
      if (imputation.method == "xgboost") {
        message("Extreme gradient tree boosting options:")
        message("    learning rate: ", eta)
        message("    regularization, minimum loss reduction (gamma): ", gamma)
        message("    maximum depth of tree: ", max_depth)
        message("    minimum sum of instance weight: ", min_child_weight)
        message("    subsample ratio for training when growing trees (prevent overfitting): ", subsample)
        message("    subsample ratio of columns when growing trees: ", colsample_bytree)
        message("    number of trees to grow per round: ", num_parallel_tree)
        message("    max number of iterations: ", nrounds)
        message("    filename for the model for periodical saving: ", save_name)
        message("    early stopping round integer: ", early_stopping_rounds, "\n")
        message("    number of threads (cpu.boost): ", cpu.boost)
      }
      
      if (imputation.method == "lightgbm") {
        message("Light Gradient Boosting Machine options:")
        message("    boosting: ", boosting)
        message("    objective: ", objective)
        message("    learning rate: ", learning_rate)
        message("    feature_fraction: ", feature_fraction)
        message("    bagging_fraction: ", bagging_fraction)
        message("    max_depth: ", max_depth)
        message("    min_data_in_leaf: ", min_data_in_leaf)
        message("    num_leaves: ", num_leaves)
        message("    early_stopping_rounds: ", early_stopping_rounds)
        message("    iteration.subsample: ", iteration.subsample)
        message("    Predictive Mean Matching: ", pmm)
        message("    max number of iterations: ", nrounds)
        message("    number of threads (cpu.boost): ", cpu.boost)
      }
      
      
      if (imputation.method == "rf") {
        message("On-the-fly-imputations options:")
        message("    number of trees to grow: ", num.tree)
        message("    minimum terminal node size: ", nodesize)
        message("    non-negative integer value used to specify random splitting: ", nsplit)
        message("    number of iterations: ", nimpute)
      }
      
      if (imputation.method == "rf_pred") {
        message("Random Forests options:")
        message("    number of trees: ", num.tree)
        message("    minimum terminal node size: ", nodesize)
        message("    non-negative integer value used to specify random splitting: ", nsplit)
        message("    number of iterations: ", nimpute)
        message("    predictive mean matching: ", pmm, "\n")
      }
      
      if (imputation.method == "bpca") {
        # message("Multiple Correspondence Analysis options:")
        # message("    number of dimentions used to predict the missing values: ", ncp)
        # message("    MCA algorithm: Regularized")
      }
      
      message("Number of CPUs used during grur's operations: ", parallel.core)
      if (!is.null(filename)) message("Filename: ", filename)
    }
    
    # Set seed for sampling reproducibility ------------------------------------
    if (is.null(random.seed)) {
      random.seed <- sample(x = 1:1000000, size = 1)
      set.seed(random.seed)
    } else {
      set.seed("num.tree")
    }
    message("\nSeed: ", random.seed)
    
    # Import data --------------------------------------------------------------
    if (is.vector(data)) {
      data <- radiator::tidy_wide(data = data, import.metadata = TRUE)
    }
    data <- dplyr::ungroup(data)
    
    # strata--------------------------------------------------------------------
    if (!is.null(strata)) {
      if (tibble::has_name(data, "POP_ID")) data <- dplyr::select(data, -POP_ID)
      
      if (is.vector(strata)) {
        message("Using strata file: ", strata)
        strata.df <- readr::read_tsv(
          file = strata, col_names = TRUE,
          col_types = readr::cols(.default = readr::col_character()))
      } else {
        message("Using strata in the global environment")
        strata.df <- strata
      }
      
      if (!tibble::has_name(strata.df, "STRATA")) {
        stop("strata file/object requires a columns named: STRATA")
      }
      
      if (!tibble::has_name(strata.df, "INDIVIDUALS")) {
        stop("strata file/object requires a columns named: INDIVIDUALS")
      }
      
      strata.df <- dplyr::distinct(strata.df, INDIVIDUALS, STRATA)
      
      # Remove potential whitespace in pop_id
      strata.df$STRATA <- radiator::clean_pop_names(strata.df$STRATA)
      
      # clean ids
      strata.df$INDIVIDUALS <- radiator::clean_ind_names(strata.df$INDIVIDUALS)
      
      # match ids in data and strata and join STRATA
      strata.df <- dplyr::filter(strata.df, INDIVIDUALS %in% unique(data$INDIVIDUALS))
      data <- dplyr::left_join(data, strata.df, by = "INDIVIDUALS")
    }
    strata.df <- NULL

    if (tibble::has_name(data, "POP_ID") && !tibble::has_name(data, "STRATA")) {
      data <- dplyr::rename(data, STRATA = POP_ID)
    }
    
    # Subsampling markers ------------------------------------------------------
    if (!is.null(subsample.markers)) {
      # subsample.markers <- 500
      sample.markers <- dplyr::distinct(data, MARKERS) %>%
        dplyr::sample_n(tbl = ., size = subsample.markers) %>%
        readr::write_tsv(
          x = .,
          path = stringi::stri_join("subsampled.markers_",
                                    subsample.markers,
                                    "_random.seed_", 
                                    random.seed, ".tsv")) %>%
        purrr::flatten_chr(.)
      data <- dplyr::filter(data, MARKERS %in% sample.markers)
      sample.markers <- NULL
    }
    
    # stats about the dataset --------------------------------------------------
    message("\nNumber of populations: ", dplyr::n_distinct(data$STRATA))
    message("Number of individuals: ", dplyr::n_distinct(data$INDIVIDUALS))
    message("Number of markers: ", dplyr::n_distinct(data$MARKERS))
    
    if (imputation.method == "rf") package.used <- "randomForestSRC"
    if (imputation.method == "xgboost") package.used <- "xgboost"
    if (imputation.method == "lightgbm") package.used <- "lightgbm"
    if (imputation.method %in% c("max", "rf_pref", "bpca")) package.used <- NULL
    
    if (!is.null(package.used)) {
      message("\nNote: If you have speed issues:")
      message("    1. Have you installed, ", package.used," package with OpenMP enabled ?")
      message("    2. Have you followed grur's vignette ?")
    }
    
    # output the proportion of missing genotypes BEFORE imputations
    na.before <- dplyr::summarise(.data = data, MISSING = round(length(GT[GT == "000000"])/length(GT), 6)) %>%
      purrr::flatten_dbl(.) %>% format(., scientific = FALSE)
    message("\nProportion of missing genotypes before imputations: ", na.before)
    
    # necessary steps to make sure we work with unique markers and not duplicated LOCUS
    if (!tibble::has_name(data, "MARKERS") && tibble::has_name(data, "LOCUS")) {
      data <- dplyr::rename(.data = data, MARKERS = LOCUS)
    }
    
    # Generate simple markers id -----------------------------------------------
    # formula in RF difficult with markers containing separators and numbers
    if (imputation.method %in% c("rf", "xgboost", "lightgbm")) {
      markers.meta <- dplyr::distinct(.data = data, MARKERS) %>%
        dplyr::mutate(NEW_MARKERS = stringi::stri_join("M", seq(1, nrow(.))))
      data <- dplyr::full_join(markers.meta, data, by = "MARKERS")
      
      markers.meta <- suppressWarnings(
        dplyr::select(
          data,
          dplyr::one_of(c("NEW_MARKERS", "MARKERS", "CHROM", "LOCUS", "POS"))) %>% 
          dplyr::distinct(NEW_MARKERS, .keep_all = TRUE))
      
      data <- dplyr::select(.data = data, -MARKERS) %>%
        dplyr::rename(MARKERS = NEW_MARKERS)
    } else {
      markers.meta <- suppressWarnings(
        dplyr::select(
          data,
          dplyr::one_of(c("MARKERS", "CHROM", "LOCUS", "POS"))) %>% 
          dplyr::distinct(MARKERS, .keep_all = TRUE))
    }
    
    # scan for REF allele column -----------------------------------------------
    if (tibble::has_name(data, "REF")) {
      ref.column <- TRUE
    } else {
      ref.column <- FALSE
    }
    
    # Check if biallelic  ------------------------------------------------------
    biallelic <- radiator::detect_biallelic_markers(data = data)
    
    # Manage Genotype Likelihood -----------------------------------------------
    if (tibble::has_name(data, "GL")) {
      if (verbose) message("\nGenotype likelihood (GL) column detected")
    }
    have <- colnames(data)
    
    # For haplotype VCF
      if (!biallelic && ref.column && tibble::has_name(data, "GT_VCF_NUC")) {
      haplo.vcf.check <- TRUE
      
      want <- c("MARKERS", "CHROM", "LOCUS", "POS", "STRATA", "INDIVIDUALS", "GT_VCF_NUC", "GL")
      selected.columns <- purrr::keep(.x = have, .p = have %in% want)
      
      data <- dplyr::select(.data = data,
                            dplyr::one_of(selected.columns)) %>%
        dplyr::mutate(GT = replace(GT_VCF_NUC, which(GT_VCF_NUC == "./."), NA)) %>%
        dplyr::select(-GT_VCF_NUC)
    } else {
      haplo.vcf.check <- FALSE
      
      want <- c("MARKERS", "CHROM", "LOCUS", "POS", "STRATA", "INDIVIDUALS", "GT", "GL")
      selected.columns <- purrr::keep(.x = have, .p = have %in% want)
      
      data <- dplyr::select(.data = data,
                            dplyr::one_of(selected.columns)) %>%
        dplyr::mutate(GT = replace(GT, which(GT == "000000"), NA))
    }
    have <- want <- selected.columns <- NULL
    
    # keep bk of stratification ------------------------------------------------
    strata.before <- dplyr::distinct(.data = data, INDIVIDUALS, STRATA)
    
    # SNP/haplotype approach -----------------------------------------------------
    
    # detect the presence of SNP/LOCUS info and combine SNPs on the same RADseq locus together
    if (tibble::has_name(data, "CHROM") && tibble::has_name(data, "LOCUS") && imputation.method != "max") {
      data <- tidyr::unite(data = data, col = CHROM_LOCUS, CHROM, LOCUS) %>%
        dplyr::select(-POS) # no longer necessary info in MARKERS
      
      # Locus with > 1 SNP
      snp.number <- dplyr::distinct(.data = data, MARKERS, CHROM_LOCUS) %>%
        dplyr::count(CHROM_LOCUS)
      
      if (max(snp.number$n) > 1) {
        # Encoding SNPs into haplotypes
        # combining SNPs into groups defined by chromosomes and locus info
        if (verbose) message("Encoding SNPs into haplotypes")
        
        if (tibble::has_name(data, "GL")) {
          keep.gl <- dplyr::ungroup(data) %>%
            dplyr::filter(!is.na(GL)) %>%
            dplyr::distinct(CHROM_LOCUS, STRATA, INDIVIDUALS, GL) %>%
            dplyr::group_by(CHROM_LOCUS, STRATA, INDIVIDUALS) %>%
            dplyr::summarise(GL = mean(GL))
          data <- dplyr::select(.data = data, -GL)
        } else {
          keep.gl <- NULL
        }
        
        locus.multiple.snp <- dplyr::filter(.data = snp.number, n > 1) %>%
          dplyr::select(CHROM_LOCUS) %>%
          purrr::flatten_chr(.)
        data.multiple.snp <- dplyr::filter(.data = data, CHROM_LOCUS %in% locus.multiple.snp)
        data.one.snp <- dplyr::filter(.data = data, !CHROM_LOCUS %in% locus.multiple.snp)
        if (length(locus.multiple.snp) > 100) {
          # parallel
          data <- list()
          data <- .grur_parallel_mc(
            X = locus.multiple.snp,
            FUN = encoding_snp,
            mc.cores = parallel.core,
            data = data.multiple.snp
          ) %>% dplyr::bind_rows(.) %>%
            dplyr::bind_rows(data.one.snp)
        } else {
          data <- purrr::map(.x = locus.multiple.snp,
                             .f = encoding_snp,
                             data = data.multiple.snp) %>%
            dplyr::bind_rows(.) %>%
            dplyr::bind_rows(data.one.snp)
        }
        
        # include GL back and use relative measure group_by locus
        if (!is.null(keep.gl)) {
          data <- dplyr::filter(.data = data, !is.na(GT)) %>%
            dplyr::select(CHROM_LOCUS, STRATA, INDIVIDUALS) %>%
            dplyr::left_join(keep.gl, by = c("CHROM_LOCUS", "STRATA", "INDIVIDUALS")) %>%
            dplyr::right_join(data, by = c("CHROM_LOCUS", "STRATA", "INDIVIDUALS")) %>%
            dplyr::select(MARKERS, CHROM_LOCUS, STRATA, INDIVIDUALS, GT, GL) %>%
            dplyr::group_by(CHROM_LOCUS) %>%
            dplyr::mutate(GL = GL/max(GL, na.rm = TRUE)) %>%
            dplyr::ungroup(.)
        }
        separate.haplo <- TRUE
        data.multiple.snp <- data.one.snp <- snp.number <- keep.gl <- NULL
      } else {
        separate.haplo <- FALSE
      }
      # End of grouping SNPs
      
      # remove unnecessary column
      data <- dplyr::select(.data = data, -CHROM_LOCUS)
    } else {
      separate.haplo <- FALSE
    }#End preparing SNP/haplo approach
    
    # Strawman imputations (max) -------------------------------------------------
    if (imputation.method == "max") {
      if (hierarchical.levels == "strata") {
        if (verbose) message("Using the most observed genotype per marker/strata for imputations")
        if (tibble::has_name(data, "GL")) {
          data.imp <- dplyr::select(data, MARKERS, STRATA, INDIVIDUALS, GT, GL) %>%
            dplyr::group_by(MARKERS, STRATA) %>%
            dplyr::mutate(
              GT = stringi::stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)),
              GT = replace(GT, which(GT == "NA"), NA),
              GL = stringi::stri_replace_na(GL, replacement = mean(GL, na.rm = TRUE)),
              GL = replace(GL, which(GL == "NA"), NA),
              GL = as.numeric(GL)) %>%
            dplyr::ungroup(.)
        } else {
          data.imp <- dplyr::select(data, MARKERS, STRATA, INDIVIDUALS, GT) %>%
            dplyr::group_by(MARKERS, STRATA) %>%
            dplyr::mutate(GT = stringi::stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)),
                          GT = replace(GT, which(GT == "NA"), NA)) %>%
            dplyr::ungroup(.)
        }
        data <- NULL
        # # detect remaining NA
        # # e.g. if one strata is missing all GT... when not using common markers
        # if (anyNA(data.one.snp)) {
        #   warning("Missing data is still present in the dataset",
        #           "\n    2 options:",
        #           "\n    run the function again with hierarchical.levels = 'global'",
        #           "\n    use common.markers = TRUE when using hierarchical.levels = 'strata'")
        # }
      }# End imputation max populations
      
      # global
      if (hierarchical.levels == "global") {
        if (verbose) message("Using the most observed genotype per marker for imputations")
        if (tibble::has_name(data, "GL")) {
          data.imp <- dplyr::select(data, MARKERS, STRATA, INDIVIDUALS, GT, GL) %>%
            dplyr::group_by(MARKERS) %>%
            dplyr::mutate(
              GT = stringi::stri_replace_na(str = GT, replacement = max(GT, na.rm = TRUE)),
              GL = as.numeric(stringi::stri_replace_na(str = GL, replacement = mean(GL, na.rm = TRUE)))
            ) %>%
            dplyr::ungroup(.)
        } else {
          data.imp <- dplyr::select(data, MARKERS, STRATA, INDIVIDUALS, GT) %>%
            dplyr::group_by(MARKERS) %>%
            dplyr::mutate(GT = stringi::stri_replace_na(str = GT, replacement = max(GT, na.rm = TRUE))) %>%
            dplyr::ungroup(.)
        }
        data <- NULL
      } # End imputation max global
    }#End imputation max
    
    # Imputation with Random Forests and tree boosting ---------------------------
    ### Note to myself: Need to add CHROM hierarchy
    ### to impute one chromosome at a time
    
    if (imputation.method %in% c("rf", "rf_pred", "xgboost", "lightgbm", "bpca")) {
      # Vector of markers
      markers.list <- unique(data$MARKERS)
      
      # Simple imputation for monomorphic markers in some pop ------------------
      
      # Problem encountered:
      # Merged SNPs might end up be missing if one of the SNP on the read is missing
      # It's usually safer to delete the whole genotype and impute it back,
      # instead of creating weird chimeras
      
      # In some dataset I've seen markers becomming monomorphic after this change
      # These markers have very very very low polymorphism and further test are
      # required to validate this technique. 
      # remedy : better MAF filtering can remove those problem...
      
      if (hierarchical.levels == "strata") {
        # 1) dont' waist time imputing
        # 2) do some screening first
        # 3) the small cost in time is worth it,
        # 4) because machine learning algorithm  in RF and xgboost will benefit having more complete 
        #    and reliable genotypes
        
        if (verbose) message("Scanning dataset for strata(s) with monomorphic marker(s)")
        # scanning for populations with one genotype group
        scan.pop <- dplyr::group_by(.data = data, MARKERS, STRATA, GT) %>%
          dplyr::tally(.)
        
        markers.pop.na <- dplyr::filter(.data = scan.pop, is.na(GT)) %>%
          dplyr::select(MARKERS, STRATA)
        
        if (nrow(markers.pop.na) > 1) {
          simple.imputation <- dplyr::filter(.data = scan.pop, !is.na(GT)) %>%
            dplyr::select(-n) %>%
            dplyr::group_by(MARKERS, STRATA) %>%
            dplyr::tally(.) %>%
            dplyr::filter(n == 1) %>%
            dplyr::select(MARKERS, STRATA) %>%
            dplyr::semi_join(markers.pop.na, by = c("MARKERS", "STRATA"))
          
          simple.imputation.number <- nrow(simple.imputation)
          
          if (simple.imputation.number > 1) {
            if (verbose) message("    Simple strawman imputations conducted on ", simple.imputation.number, " markers/pops combo")
            # update markers.list
            markers.list <- dplyr::anti_join(markers.pop.na, simple.imputation, by = c("MARKERS", "STRATA")) %>%
              dplyr::ungroup(.) %>%
              dplyr::distinct(MARKERS) %>%
              purrr::flatten_chr(.)
            
            data.imp <- dplyr::inner_join(
              data, simple.imputation, by = c("MARKERS", "STRATA")) %>%
              dplyr::group_by(MARKERS, STRATA) %>%
              dplyr::mutate(
                GT = stringi::stri_replace_na(
                  GT, replacement = max(GT, na.rm = TRUE))) %>%
              dplyr::ungroup(.)
            
            # if GL is present give the mean value for the imputed genotype
            # note to myself: update doc to say what you're doing here...
            if (tibble::has_name(data.imp, "GL")) {
              data.imp <- data.imp %>%
                dplyr::group_by(MARKERS, STRATA) %>%
                dplyr::mutate(GL = as.numeric(stringi::stri_replace_na(GL, replacement = mean(GL, na.rm = TRUE)))) %>%
                dplyr::ungroup(.)
            }
            # update dataframe
            data <- dplyr::anti_join(
              data, simple.imputation, by = c("MARKERS", "STRATA")) %>%
              dplyr::bind_rows(data.imp) %>%
              dplyr::arrange(MARKERS, STRATA, INDIVIDUALS)
          }
        }
        
        # removed unused object
        data.imp <- simple.imputation <- simple.imputation.number <- NULL
        scan.pop <- markers.pop.na <- NULL
      }# End imputation prep by pop
      
      
      # after some testing the next step is not warranted 
      # prefer to let the algorithm decide how to impute those one. 
      if (hierarchical.levels == "not_needed") {
        # if (hierarchical.levels == "global") {
        scan.markers.na <- dplyr::filter(data, is.na(GT)) %>%
          dplyr::distinct(MARKERS) %>%
          purrr::flatten_chr(.)
        
        if (length(scan.markers.na) < length(markers.list)) {
          
          simple.imputation <- dplyr::group_by(.data = data, MARKERS, GT) %>%
            dplyr::tally(.) %>%
            dplyr::filter(!is.na(GT)) %>%
            dplyr::select(-n) %>%
            dplyr::group_by(MARKERS) %>%
            dplyr::tally(.) %>%
            dplyr::filter(n == 1) %>%
            dplyr::select(MARKERS) %>%
            purrr::flatten_chr(.)
          
          simple.imputation.number <- length(simple.imputation)
          if (simple.imputation.number >= 1) {
            if (verbose) message("Simple imputations conducted on ", simple.imputation.number, " markers")
            
            data.imp <- dplyr::filter(data, MARKERS %in% simple.imputation) %>%
              dplyr::group_by(MARKERS) %>%
              dplyr::mutate(
                GT = stringi::stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
              dplyr::ungroup(.)
            
            # if GL is present give the mean value for the imputed genotype
            if (tibble::has_name(data.imp, "GL")) {
              data.imp <- data.imp %>%
                dplyr::group_by(MARKERS) %>%
                dplyr::mutate(
                  GL = as.numeric(stringi::stri_replace_na(GL, replacement = mean(GL, na.rm = TRUE)))) %>%
                dplyr::ungroup(.)
            }
            
            data <- dplyr::filter(data, !MARKERS %in% simple.imputation) %>%
              dplyr::bind_rows(data.imp)
            
            data.imp <- NULL
            
            markers.list <- dplyr::filter(data, is.na(GT)) %>%
              dplyr::distinct(MARKERS) %>%
              purrr::flatten_chr(.)
            
            # data.imp.bk <- dplyr::filter(data, !MARKERS %in% markers.list)
          } else {
            markers.list <- scan.markers.na
            # data.imp.bk <- dplyr::filter(data, !MARKERS %in% markers.list)
          }
        }
        # else {
        # data.imp.bk <- NULL
        # }
        scan.markers.na <- NULL
      }
      
      
      
      # On-the-fly-imputations using Random Forests ------------------------------
      if (imputation.method == "rf") {
        if (verbose) message("On-the-fly-imputations using Random Forests algorithm")
        
        # Parallel computations options for randomForestSRC
        rf.option <- getOption("rf.cores")
        mc.option <- getOption("mc.cores")
        options(rf.cores = parallel.core, mc.cores = parallel.core)
        
        # prepare data
        data.imp <- dplyr::select(data, MARKERS, STRATA, INDIVIDUALS, GT) %>%
          dplyr::group_by(STRATA, INDIVIDUALS) %>%
          tidyr::spread(data = ., key = MARKERS, value = GT) %>%
          dplyr::ungroup(.) %>%
          dplyr::mutate_all(.tbl = ., .funs = factor)
        
        # Random Forest by pop
        if (hierarchical.levels == "strata") {
          if (verbose) message("    Imputations computed by strata, take a break...")
          data.imp <- split(x = data.imp, f = data.imp$STRATA) %>%
            purrr::map_df(
              .x = ., .f = impute_rf,
              num.tree = num.tree, nodesize = nodesize, nsplit = nsplit,
              nimpute = nimpute,
              verbose = FALSE,
              hierarchical.levels = "strata") %>%
            dplyr::mutate_all(.tbl = ., .funs = as.character) %>%
            tidyr::gather(data = ., key = MARKERS, value = GT, -INDIVIDUALS) %>%
            dplyr::right_join(strata.before, by = "INDIVIDUALS") %>%
            dplyr::arrange(MARKERS, STRATA, INDIVIDUALS)
        }#End RF by pop
        
        # Random Forests global
        if (hierarchical.levels == "global") { # Globally/overall
          if (verbose) message("    Imputations computed globally, take a break...")
          data.imp <- impute_rf(
            x = data.imp,
            num.tree = num.tree, nodesize = nodesize, nsplit = nsplit,
            nimpute = nimpute, verbose = FALSE,
            hierarchical.levels = "global") %>%
            dplyr::mutate_all(.tbl = ., .funs = as.character) %>%
            tidyr::gather(data = ., key = MARKERS, value = GT, -c(STRATA, INDIVIDUALS)) %>%
            dplyr::mutate(STRATA = factor(STRATA)) %>% 
            dplyr::arrange(MARKERS, STRATA, INDIVIDUALS)
        } #End RF global
        
        
        # separate the haplotypes/snp group
        # separating SNPs on the same locus and chromosome
        # back to original data format
        if (separate.haplo) {
          if (verbose) message("Decoding haplotypes")
          data.imp <- decoding_haplotypes(
            data = data.imp, parallel.core = parallel.core)
        }
        
        
        # Restore original Parallel computations options
        options(rf.cores = rf.option, mc.cores = mc.option)
      }# End on-the-fly-imputations using RF
      
      # Random Forest imputation as a prediction problem -------------------------
      if (imputation.method == "rf_pred") {
        if (verbose) message("Using Random Forests algorithm as a prediction problem, take a break...")
        
        if (hierarchical.levels == "strata") {
          data.imp <- purrr::map(.x = data,
                                 .f = grur_imputer,
                                 hierarchical.levels = hierarchical.levels,
                                 num.tree = num.tree,
                                 pmm = pmm,
                                 random.seed = random.seed,
                                 parallel.core = parallel.core) %>%
            dplyr::bind_rows(.)
        }
        
        # Random Forests global
        if (hierarchical.levels == "global") { # Globally/overall
          # if (verbose) message("Imputations computed globally, take a break...")
          data.rf.imp <- list() # to store results
          data.rf.imp <- grur_imputer(data = data,
                                      hierarchical.levels = hierarchical.levels,
                                      num.tree = num.tree,
                                      pmm = pmm,
                                      random.seed = random.seed,
                                      parallel.core = parallel.core)
        } # End imputation RF global
        
        
      }# End rf_pred
      
      # XGBoost & LightGBM ------------------------------------------------------
      if (imputation.method %in% c("xgboost", "lightgbm")) {
        # prep data
        if (hierarchical.levels == "strata") {
          data <- dplyr::ungroup(data) %>%
            dplyr::arrange(MARKERS, STRATA, INDIVIDUALS) %>%
            dplyr::mutate(
              # GT_N = as.numeric(factor(replace(GT, which(is.na(GT)), 0))),
              POP_ID_N = as.numeric(factor(STRATA)),
              INDIVIDUALS_N = as.numeric(factor(INDIVIDUALS))
            ) %>%
            dplyr::group_by(MARKERS) %>%
            dplyr::mutate(GT_N = factorize_gt(GT)) %>%
            dplyr::ungroup(.)
            readr::write_tsv(x = data, path = "imputation_factor_dictionary.rad")
          #fst::write.fst(x = data, path = "imputation_factor_dictionary.rad")
          data.boost <- data %>%
            dplyr::select(MARKERS, STRATA = POP_ID_N, INDIVIDUALS = INDIVIDUALS_N, GT = GT_N) %>%
            dplyr::arrange(MARKERS, STRATA, INDIVIDUALS) %>%
            dplyr::group_by(STRATA, INDIVIDUALS) %>%
            tidyr::spread(data = ., key = MARKERS, value = GT) %>%
            dplyr::ungroup(.) %>% 
            as.matrix(.) %>% Matrix::Matrix(., sparse = TRUE)
        }
        
        if (hierarchical.levels == "global") {
          data <- dplyr::ungroup(data) %>%
            dplyr::arrange(MARKERS, STRATA, INDIVIDUALS) %>%
            dplyr::select(-STRATA) %>%
            dplyr::mutate(INDIVIDUALS_N = as.numeric(factor(INDIVIDUALS))) %>%
            dplyr::group_by(MARKERS) %>%
            dplyr::mutate(GT_N = factorize_gt(GT)) %>%
            dplyr::ungroup(.)
          #fst::write.fst(x = data, path = "imputation_factor_dictionary.rad")
          
          data.boost <- data %>%
            dplyr::select(MARKERS, INDIVIDUALS = INDIVIDUALS_N, GT = GT_N) %>%
            dplyr::arrange(MARKERS, INDIVIDUALS) %>%
            dplyr::group_by(INDIVIDUALS) %>%
            tidyr::spread(data = ., key = MARKERS, value = GT) %>%
            dplyr::ungroup(.) %>% 
            as.matrix(.) %>% 
            Matrix::Matrix(., sparse = TRUE)
        }
        data <- gl.wide <- NULL
        
        if (imputation.method == "xgboost") {
          if (verbose) message("Using extreme gradient tree boosting algorithm, take a break...")
          markers.df <- tibble::data_frame(MARKERS = markers.list) %>% 
            dplyr::mutate(SPLIT_VEC = split_imp(x = ., cpu.rounds = 20, parallel.core = cpu.boost))
          boost.split <- unique(markers.df$SPLIT_VEC)
          # single thread test e.g.
          # cpu.boost <- 1
          # cpu.boost <- 2
          # cpu.boost <- 4
          # cpu.boost <- 8
          # system.time(data.imp <- purrr::map_df(
          #   .x = boost.split,
          #   .f = grur_boost_imputer,
          #   markers.df = markers.df,
          #   data.xgb = data.boost,
          #   gl.wide = gl.wide,
          #   eta = eta,
          #   gamma = gamma,
          #   max_depth = max_depth,
          #   min_child_weight = min_child_weight,
          #   subsample = subsample,
          #   colsample_bytree = colsample_bytree,
          #   num_parallel_tree = num_parallel_tree,
          #   cpu.boost = cpu.boost,
          #   nrounds = nrounds,
          #   early_stopping_rounds = early_stopping_rounds,
          #   save_name = save_name))
          
          # parallel testing core option
          # cpu.boost <- 1
          # cpu.boost <- 2
          # cpu.boost <- 4
          # cpu.boost <- 8
          # parallel.core <- 2
          # parallel.core <- 4
          # parallel.core <- 8
          
          markers.df <- tibble::data_frame(MARKERS = markers.list) %>%
            dplyr::mutate(SPLIT_VEC = split_imp(x = ., cpu.rounds = 2, 
                                                parallel.core = parallel.core))
          boost.split <- unique(markers.df$SPLIT_VEC)
          
          data.imp <- list()
          data.imp <- .grur_parallel_mc(
            X = boost.split,
            FUN = grur_boost_imputer,
            mc.cores = parallel.core,
            markers.df = markers.df,
            data.xgb = data.boost,
            gl.wide = gl.wide,
            eta = eta,
            gamma = gamma,
            max_depth = max_depth,
            min_child_weight = min_child_weight,
            subsample = subsample,
            colsample_bytree = colsample_bytree,
            num_parallel_tree = num_parallel_tree,
            cpu.boost = cpu.boost,
            nrounds = nrounds,
            early_stopping_rounds = early_stopping_rounds,
            save_name = save_name) %>% 
            dplyr::bind_rows(.)
        }
        if (imputation.method == "lightgbm") {
          if (verbose) message("Using Light Gradient Boosting Machine algorithm, take a break...")
          markers.df <- tibble::data_frame(MARKERS = markers.list) %>%
            dplyr::mutate(
              SPLIT_VEC = as.integer(floor((10 * (1:nrow(.) - 1) / 
                                              nrow(.)) + 1)))
          boost.split <- unique(markers.df$SPLIT_VEC)
          
          message("    Imputations conducted in 10 rounds")
          # testing serial
          # cpu.boost <- 4
          # boosting <- "dart"
          # learning_rate <- 0.1
          # bagging_freq <- 1
          # max_depth <- 9
          # num_leaves <- 512
          
          data.imp <- purrr::map_df(
            .x = boost.split,
            .f = grur_lgbm_imputer,
            markers.df = markers.df,
            data.gbm = data.boost,
            boosting = boosting,
            objective = objective,
            learning_rate = learning_rate,
            feature_fraction = feature_fraction,
            bagging_fraction = bagging_fraction,
            bagging_freq = bagging_freq,
            max_depth = max_depth,
            min_data_in_leaf = min_data_in_leaf,
            num_leaves = num_leaves,
            cpu.boost = cpu.boost,
            early_stopping_rounds = early_stopping_rounds,
            nrounds = nrounds,
            iteration.subsample = iteration.subsample,
            pmm = pmm)
          
          
          # LightGBM doesn't work well in parallel for features
          # data.imp <- list()
          # data.imp <- .grur_parallel_mc(
          #   X = boost.split,
          #   FUN = grur_lgbm_imputer,
          #   mc.cores = parallel.core,
          #   markers.df = markers.df,
          #   data.gbm = data.boost,
          #   boosting = boosting,
          #   objective = objective,
          #   learning_rate = learning_rate,
          #   feature_fraction = feature_fraction,
          #   bagging_fraction = bagging_fraction,
          #   max_depth = max_depth,
          #   min_data_in_leaf = min_data_in_leaf,
          #   num_leaves = num_leaves,
          #   cpu.boost = 2,
          #   early_stopping_rounds = early_stopping_rounds,
          #   nrounds = nrounds) %>% 
          #   dplyr::bind_rows(.)
        }
        
        data.boost <- NULL# remove unused objects
        
        # save.image("imputation.test.RData")
        #Note to myself: if you want to further explore error
        # the column need parsing...
        # error.lightgbm <-  dplyr::select(data.imp, MARKERS, ERROR)
        data.imp <- data.imp %>% dplyr::select(IMPUTED_DATA) %>% tidyr::unnest(.)
        
        # Remove factor/integer type genotype (revert back to original)
        data.imp <- defactorize_gt(data.to.change = data.imp)
        
        # Reintroduce the stratification (check if required)
        data.imp <- dplyr::left_join(strata.before, data.imp, by = "INDIVIDUALS") %>%
          dplyr::arrange(MARKERS, STRATA, INDIVIDUALS)
        strata.before <- NULL# remove unused objects
        
        if (separate.haplo) {
          # separate the haplotypes/snp group
          # Decoding haplotypes
          # separating SNPs on the same locus and chromosome, back to original data format
          if (verbose) message("Decoding haplotypes")
          data.imp <- decoding_haplotypes(
            data = data.imp, parallel.core = parallel.core)
        }
        
        
        # Impute GL
        #Note: will no longer work with stacks
        # if (tibble::has_name(data.imp, "GL") && hierarchical.levels == "strata") {
        #   message("Imputing GL with mean value per populations")
        #   data.imp <- dplyr::group_by(.data = data.imp, MARKERS, STRATA, GT) %>%
        #     dplyr::mutate(
        #       GL = stringi::stri_replace_na(GL, replacement = mean(GL, na.rm = TRUE)),
        #       GL = replace(GL, which(GL %in% c("NA", "NaN")), NA),
        #       GL = as.numeric(GL)) %>%
        #     dplyr::group_by(MARKERS, GT) %>%
        #     dplyr::mutate(
        #       GL = as.numeric(stringi::stri_replace_na(GL, replacement = mean(GL, na.rm = TRUE)))
        #     ) %>%
        #     dplyr::ungroup(.)
        # }
        # if (tibble::has_name(data.imp, "GL") && hierarchical.levels == "global") {
        #   message("Imputing GL with mean overall value")
        #   data.imp <- dplyr::group_by(.data = data.imp, MARKERS, GT) %>%
        #     dplyr::mutate(
        #       GL = as.numeric(stringi::stri_replace_na(GL, replacement = mean(GL, na.rm = TRUE)))
        #     ) %>%
        #     dplyr::ungroup(.)
        # }
      }# End boost
      
      
      # Baysian PCA here -------------------------------------------------------
      # if (imputation.method == "bpca") {
      # 
      # }
      
      
      # }# End mca
      
    } # End imputation RF, xgboost lightgbm and MCA
    
    # prep results ---------------------------------------------------------------
    
    # Replace NA by 000000 in GT column if found
    if (anyNA(data.imp)) {
      warning("Missing data is still present in the dataset",
              "\n    2 options:",
              "\n    run the function again with hierarchical.levels = 'global'",
              "\n    use common.markers = TRUE when using hierarchical.levels = 'strata'")
      if (haplo.vcf.check) {
        data.imp$GT <- stringi::stri_replace_na(str = data.imp$GT, replacement = "./.")
      } else {
        data.imp$GT <- stringi::stri_replace_na(str = data.imp$GT, replacement = "000000")
      }
    }
    
    if (haplo.vcf.check) {
      data.imp <- dplyr::rename(data.imp, GT_VCF_NUC = GT)
    }
    
    # Compute REF/ALT allele... might have change depending on prop of missing values
    if (verbose) message("Adjusting REF/ALT alleles to account for imputations...")
    data.imp <- radiator::change_alleles(
      data = dplyr::rename(data.imp, POP_ID = STRATA),
      biallelic = biallelic,
      parallel.core = parallel.core,
      verbose = verbose)$input
    
    if (tibble::has_name(data.imp, "POLYMORPHIC.x")) data.imp <- dplyr::select(data.imp, -POLYMORPHIC.x)
    if (tibble::has_name(data.imp, "POLYMORPHIC.y")) data.imp <- dplyr::select(data.imp, -POLYMORPHIC.y)
    
    data.imp <- dplyr::rename(data.imp, STRATA = POP_ID)
    # Integrate markers.meta columns and sort
    if (!is.null(markers.meta)) {
      want <- c( "MARKERS", "CHROM", "LOCUS", "POS", "STRATA",
                 "INDIVIDUALS", "REF", "ALT", "GT", "GT_VCF",
                 "GT_VCF_NUC", "GT_BIN", "GL")
      
      if (tibble::has_name(markers.meta, "NEW_MARKERS")) {
        data.imp <- suppressWarnings(dplyr::left_join(
          dplyr::rename(data.imp, NEW_MARKERS = MARKERS),
          markers.meta, by = "NEW_MARKERS") %>%
            dplyr::select(dplyr::one_of(want)))
      } else {
        data.imp <- suppressWarnings(dplyr::left_join(
          data.imp, markers.meta, by = "MARKERS") %>%
            dplyr::select(dplyr::one_of(want)))
      }
      want <- markers.meta <- NULL
    } else {
      data.imp <- dplyr::arrange(.data = data.imp, MARKERS, STRATA, INDIVIDUALS)
    }
    
    # Write to working directory
    if (!is.null(filename)) {
      if (is.null(subsample.markers)) {
        tidy.name <- stringi::stri_join(filename, ".rad")
      } else {
        tidy.name <- stringi::stri_join(
          filename, "_subsample.markers_", subsample.markers, ".rad")
      }
      if (verbose) message("Writing the imputed data: \n", tidy.name)
      #fst::write.fst(x = data.imp, path = tidy.name, compress = 85)
    }
    
    # Missing after imputation:
    na.after <- dplyr::summarise(.data = data.imp, MISSING = round(length(GT[GT == "000000"])/length(GT), 6)) %>%
      purrr::flatten_dbl(.) %>% format(., scientific = FALSE)
    message("\nProportion of missing genotypes after imputations: ", na.after)
    
    # Error notices
    if (imputation.method == "xgboost" && file.exists("grur_imputations_error.txt")) {
      message("Error notice: Tree boosting imputations encountered error(s),
              please look in the working directory for a file:
              grur_imputations_error.txt
              email the problem to the author: thierrygosselin@icloud.com")
    }
  }
  if (verbose) {
    # output the proportion of missing genotypes after imputations
    timing <- proc.time() - timing
    message("\nComputation time: ", round(timing[[3]]), " sec")
    cat("################## grur::grur_imputations completed ###################\n")
  }
  options(width = opt.change)
  return(data.imp)
} # End imputations

# Internal nested functions ----------------------------------------------------
# on-the-fly-imputations with randomForestSRC package --------------------------

#' @title impute_rf
#' @description on-the-fly-imputations using randomForestSRC package
#' @rdname impute_rf
#' @keywords internal
#' @export

impute_rf <- function(
  x,
  num.tree = 10,
  nodesize = 1,
  # splitrule = "random",
  nsplit = 10,
  nimpute = 10,
  verbose = FALSE,
  hierarchical.levels = "strata") {
  
  if (hierarchical.levels == "strata") {
    message("        Imputations for strata: ", unique(x$STRATA))
    x <- dplyr::select(x, -STRATA)
  }
  
  res <- randomForestSRC::impute.rfsrc(
    data = data.frame(x),
    ntree = num.tree,
    nodesize = nodesize,
    # splitrule = "random",
    nsplit = nsplit, #split.number,
    nimpute = nimpute, #iteration.rf,
    do.trace = verbose)
  return(res)
} # End on-the-fly imputation function


# grur_imputer ---------------------------------------------------------------
#' @title grur_imputer
#' @description imputations using Ranger package and predictive mean matching
#' @rdname grur_imputer
#' @keywords internal
#' @export

grur_imputer <- function(
  data,
  num.tree = 100,
  pmm = 0,
  random.seed = NULL,
  parallel.core = parallel::detectCores() - 1,
  # markers.linkage = "multivariate",
  hierarchical.levels = "strata",
  markers.list = markers.list,
  verbose = verbose
) {
  # data <- data #test
  
  data <- data %>%
    dplyr::select(STRATA, INDIVIDUALS, MARKERS, GT) %>%
    dplyr::group_by(INDIVIDUALS, STRATA) %>%
    dplyr::mutate(GT = replace(GT, which(is.na(GT)), "missing")) %>%
    tidyr::spread(data = ., key = MARKERS, value = GT) %>%
    dplyr::ungroup(.)
  
  # data[is.na(data.model)] <- "missing"
  
  
  data.gl <- NULL
  # if (tibble::has_name(data.pop, "GL")) {
  #   # not sure useful
  #   data.gl <- data.pop %>%
  #     dplyr::select(STRATA, INDIVIDUALS, MARKERS, GL_RF) %>%
  #     dplyr::group_by(INDIVIDUALS, STRATA) %>%
  #     tidyr::spread(data = ., key = MARKERS, value = GL_RF) %>%
  #     dplyr::ungroup(.)
  # } else {
  #   data.gl <- NULL
  # }
  
  
  data.na <- NULL # remove after test
  
  
  # initiate while loop
  i <- 1
  pred.error <- rep(1, length(markers.list))
  names(pred.error) <- markers.list
  oob.error <- TRUE
  maxiter <- 10000
  data.imp <- tibble::data_frame(character(0))
  
  # imp <- list()
  while (oob.error && i <= maxiter) {
    data.last <- data
    pred.error.last <- pred.error
    
    data.rf <- list()
    data.rf <- .grur_parallel_mc(
      X = markers.list,
      FUN = impute_genotypes,
      mc.cores = parallel.core,
      data = data,
      data.na = data.na,
      data.gl = data.gl,
      num.tree = num.tree,
      pmm = pmm,
      random.seed = random.seed,
      parallel.core = parallel.core,
      hierarchical.levels = hierarchical.levels,
      # markers.linkage = markers.linkage,
      pred.error = pred.error
    ) #%>%
    # dplyr::bind_rows(.)
    
    system.time(test <- purrr::map(
      .x = markers.list, .f = impute_genotypes,
      data = data,
      data.na = data.na,
      data.gl = data.gl,
      data.imp = data.imp,
      num.tree = num.tree,
      pmm = pmm,
      random.seed = random.seed,
      parallel.core = parallel.core,
      hierarchical.levels = hierarchical.levels,
      # markers.linkage = markers.linkage,
      pred.error = pred.error))
    
    
    test <- dplyr::bind_rows(data.imp)
    
    # update error
    oob.error <- mean(pred.error) < mean(pred.error.last)
    i <- i + 1 # update iteration
  } # End of loop
  
  if (i == maxiter && oob.error || i == 2) {
    imputed.dataset <- data
  } else {
    imputed.dataset <- data.last
  }
  
  imputed.dataset <- dplyr::mutate_all(.tbl = imputed.dataset,
                                       .funs = as.character, exclude = NA)
  
  
  # results --------------------------------------------------------------------
  
  
  return(data.imp)
} #End grur_imputer

# impute_genotypes -------------------------------------------------------------
#' @title impute_genotypes
#' @description imputations using Ranger package and predictive mean matching of missRanger
#' @rdname impute_genotypes
#' @keywords internal
#' @export

impute_genotypes <- function(
  markers.list,
  data,
  data.na,
  data.gl = NULL,
  data.imp,
  num.tree = 100,
  pmm = 0,
  random.seed = NULL,
  parallel.core = parallel::detectCores() - 1,
  hierarchical.levels = "strata",
  # markers.linkage = "multivariate",
  pred.error = pred.error
) {
  # m <- "BINDED_M7114_M7115_M7116_M7117"
  # m <- "BINDED_M1_M2_M3_M4_M5"
  # m <- "BINDED_M101_M102"
  # m <- "M993"
  m <- markers.list
  message("Marker: ", m)# for diagnostic
  
  # Handling complete and missing data ---------------------------------------
  data.model <- dplyr::filter(.data = data, rlang::.data[[m]] != "missing") %>%
    dplyr::mutate_all(.tbl = ., .funs = factor)
  data.missing <- dplyr::filter(.data = data, rlang::.data[[m]] == "missing") %>%
    dplyr::select(-dplyr::one_of(m))
  
  # If all missing screening # this should be done with markers.list before all this
  # if (nrow(data.missing) > 0) {
  
  # GL
  data.gl <- NULL
  
  data.complete <- NULL # remove after test
  
  if (!is.null(data.gl)) {
    # mean GL per sample
    case.weights <- suppressWarnings(
      dplyr::select(.data = data.complete, INDIVIDUALS) %>%
        dplyr::left_join(data.gl, by = "INDIVIDUALS") %>%
        dplyr::ungroup(.) %>%
        dplyr::select(-dplyr::one_of(c("STRATA", "INDIVIDUALS"))) %>%
        purrrlyr::invoke_rows(.f = purrr::lift_vd(mean), .to = "GL", .collate = "cols") %>%
        dplyr::select(GL) %>%
        purrr::flatten_dbl(.)
    )
    # mean GL per markers
    split.select.weights <- suppressWarnings(
      dplyr::select(.data = data.complete, INDIVIDUALS) %>%
        dplyr::left_join(data.gl, by = "INDIVIDUALS") %>%
        dplyr::ungroup(.) %>%
        dplyr::select(-dplyr::one_of(c(m, "STRATA", "INDIVIDUALS"))) %>%
        dplyr::summarise_all(.tbl = ., .funs = mean) %>%
        purrr::flatten_dbl(.)
    )
  } else {
    case.weights <- NULL
    split.select.weights <- NULL
  }
  message("Data preparation: ok")# for diagnostic
  
  # Formula ------------------------------------------------------------------
  # if (markers.linkage == "multivariate") {
  # if (hierarchical.levels == "strata") {
  # discard.columns <- c(m, "STRATA", "INDIVIDUALS")
  # discard.columns <- c(m, "STRATA")
  # model.columns <- setdiff(colnames(data.complete), discard.columns)
  model.columns <- setdiff(colnames(data.model), m)
  rf.formula <- stats::as.formula(
    stringi::stri_join(m, " ~ ",
                       stringi::stri_join(model.columns, collapse = "+")))
  always.split.variables <- NULL
  always.split.variables <- "STRATA"
  # } else {
  #   discard.columns <- c(m, "INDIVIDUALS")
  #   model.columns <- setdiff(colnames(data.complete), discard.columns)
  #   model.columns <- setdiff(colnames(data.complete), m)
  #   rf.formula <- stats::as.formula(
  #     stringi::stri_join(m, " ~ ",
  #                        stringi::stri_join(model.columns, collapse = "+")))
  #   # rf.formula <- stats::reformulate(termlabels = "STRATA", response = m)
  #   always.split.variables <- c("STRATA")
  # }
  # } else {#univariate (one marker at a time)
  #   if (hierarchical.levels == "strata") {
  #     rf.formula <- stats::reformulate(termlabels = ".", response = m)
  #     always.split.variables <- NULL
  #   } else {
  #     rf.formula <- stats::reformulate(termlabels = "STRATA", response = m)
  #     always.split.variables <- c("STRATA")
  #   }
  # }#End composing formula
  message("Formula: ok")# for diagnostic
  # RF -----------------------------------------------------------------------
  ranger.res <- ranger::ranger(
    formula = rf.formula,
    data = data.model,
    # num.trees = 1000,
    num.trees = num.tree,
    case.weights = case.weights,
    split.select.weights = split.select.weights,
    always.split.variables = always.split.variables,
    num.threads = 1,
    # num.threads = parallel.core,
    seed = random.seed)
  
  # ranger.res
  predicted <- stats::predict(ranger.res$forest, data.missing)$predictions
  
  message("RF: ok")# for diagnostic
  # predictive mean matching ---------------------------------------------
  if (pmm > 0) {
    
    # ytrain <- dplyr::select(data.complete, GT) %>%
    #   dplyr::mutate(GT = as.character(GT)) %>%
    #   purrr::flatten_chr(.)
    
    # wide format
    # ytrain <- dplyr::select(.data = data.complete, dplyr::one_of(m)) %>%
    #   dplyr::ungroup(.) %>%
    #   dplyr::mutate_all(.tbl = ., .funs = as.character) %>%
    #   purrr::flatten_chr(.)
    
    # long format (doesn't give reliable results with missRanger...)
    # ytrain <- dplyr::select(.data = data.complete, GT_IMP) %>%
    #   dplyr::ungroup(.) %>%
    #   dplyr::mutate_all(.tbl = ., .funs = as.character) %>%
    #   purrr::flatten_chr(.)
    
    # To pass Travis
    # predicted <- missRanger::pmm(
    #   xtrain = ranger.res$predictions,
    #   xtest = predicted,
    #   ytrain = ytrain,
    #   k = pmm)
    predicted <- NULL
  }
  
  message("pred.mean.matching: ok")# for diagnostic
  
  # data.missing[,m] <- predicted
  # data <- suppressWarnings(dplyr::bind_rows(data.complete, data.missing) %>% dplyr::arrange(INDIVIDUALS))
  # data.imp <- data[,m]
  
  
  # imp <- dplyr::select(.data = data2, dplyr::one_of(c("STRATA", "INDIVIDUALS", m))) %>%
  #   decoding_haplotypes(parallel.core = parallel.core)
  # imp[[m]] <- decoding_haplotypes(data = data.imp, parallel.core = parallel.core)
  
  data.imp <- dplyr::select(.data = data.missing, INDIVIDUALS) %>%
    dplyr::mutate(
      MARKERS = rep(m, nrow(data.missing)),
      GT = predicted
    )
  
  
  
  pred.error[[m]] <- ranger.res$prediction.error
  if (is.nan(pred.error[[m]])) pred.error[[m]] <- 0
  
  res <- list(data = data, pred.error = pred.error, data.imp = data.imp)
  # } else {
  #   pred.error[[m]] <- 0
  #   data.imp <- data[,m]
  #   res <- list(data = data, pred.error = pred.error, data.imp = data.imp)
  # }
  # difference with missForest, missRanger and randomForestSRC:
  # other package are updating the dataset with the imputed value and continue
  # looping through the columns, creating a bias or differences between columns
  # imputed first through last as none have the same predictor columns missigness
  # e.g. missRanger: completed <- union(completed, v)
  # I think it's preferable to leave the data frame as is and merge columns in the end
  
  # if (separate.haplo) {
  #   imputed.dataset <- imputed.dataset %>%
  #     tidyr::separate(data = ., col = GT, into = haplo.meta, sep = "-", remove = TRUE) %>%
  #     tidyr::gather(data = ., key = MARKERS, value = GT, -c(CHROM_LOCUS, STRATA, INDIVIDUALS))
  # }
  return(res)
  message("results: ok")# for diagnostic
} #End impute_genotypes

# grur_boost_imputer ---------------------------------------------------------------
#' @title grur_boost_imputer
#' @description imputations using eXtreme Gradient Boosting Tree
#' @rdname grur_boost_imputer
#' @keywords internal
#' @export
grur_boost_imputer <- function(
  boost.split = NULL,
  markers.df = NULL,
  data.xgb = NULL,
  gl.wide = NULL,
  eta = 0.2,
  gamma = 0,
  max_depth = 6,
  min_child_weight = 1,
  subsample = 0.8,
  colsample_bytree = 1,
  num_parallel_tree = 1,
  cpu.boost = 1,
  nrounds = 200,
  early_stopping_rounds = 20,
  save_name = "imputation.model.temp"
) {
  # boost.split <- 3
  # markers.list <- dplyr::distinct(markers.df, MARKERS) %>% purrr::flatten_chr(.)
  markers.list <- dplyr::filter(markers.df, SPLIT_VEC == boost.split) %>%
    dplyr::select(MARKERS) %>% purrr::flatten_chr(.)
  
  xgb_imp <- function(
    markers.list = NULL, data.xgb = NULL, gl.wide = NULL,
    eta = 0.2, gamma = 0, max_depth = 6, min_child_weight = 1,
    subsample = 0.8, colsample_bytree = 1, num_parallel_tree = 1,
    nthread = 1, nrounds = 200, early_stopping_rounds = 20,
    save_name = "imputation.model.temp"
  ) {
    # m <- markers.list <- "BINDED_M10642_M10643"
    # m <- markers.list <- "BINDED_M10638_M10639_M10640"
    # m <- markers.list <- "M29"
    m <- markers.list
    message("Imputation of marker: ", m)
    # preparing data
    
    all.var <- colnames(data.xgb)
    train.var <- !all.var %in% c(m)
    data.na <- is.na(data.xgb[, all.var, drop = FALSE])
    select.na <- data.na[, m]
    all.var <- data.na <- NULL
    
    train.data <- data.xgb[!select.na, train.var]
    train.label <- data.xgb[!select.na, m]
    
    data.xgb <- xgboost::xgb.DMatrix(data = train.data, label = train.label, missing = NA)
    train.data <- NULL
    
    train.missing <- data.xgb[select.na, train.var, drop = FALSE]
    if (nrow(train.missing) == 0) stop("code error: email author")
    train.missing.label <- data.xgb[select.na, m]
    train.var <- select.na <- NULL
    id.string <- train.missing[, "INDIVIDUALS"]
    train.missing <- xgboost::xgb.DMatrix(
      data = train.missing, label = train.missing.label, missing = NA)
    train.missing.label <- NULL
    # prepare table
    res <- tibble::data_frame(
      INDIVIDUALS = id.string,
      MARKERS = rep(m, length(id.string)))
    
    params <- list(
      booster = "dart",
      silent = 0,
      eta = eta,
      gamma = gamma,
      max_depth = max_depth,
      min_child_weight = min_child_weight,
      subsample = subsample,
      colsample_bytree = colsample_bytree,
      num_parallel_tree = num_parallel_tree,
      objective = "multi:softmax",
      num_class = length(unique(train.label)), #max(data.label) + 1,
      nthread = nthread #XGBoost is actually faster when set to 1 for most dataset...
    )
    
    watchlist <- list(train = data.xgb) # not an argument of xgboost::xgboost
    
    callbacks <- list(xgboost::cb.early.stop(
      stopping_rounds = early_stopping_rounds,
      metric_name = "train_merror", verbose = FALSE))
    
    # Catch error while tree boosting
    safe_boost <- purrr::safely(.f = xgboost::xgb.train)
    
    model <- safe_boost(
      data = data.xgb,
      missing = NA,
      weight = NULL,
      params = params,
      nrounds = nrounds,
      verbose = 0,
      watchlist = watchlist,
      callbacks = callbacks,
      save_period = 0,
      save_name = save_name)
    
    if (is.null(model$error)) {
      res$GT <- stats::predict(
        model$result, train.missing,
        ntreelimit = model$result$best_ntreelimit,
        missing = NA)
    } else {
      readr::write_lines(x = boost.res$error, path = "grur_imputations_error.txt", append = TRUE)
      res$GT <- as.numeric(rep(NA, nrow(res)))
    }
    
    # Unused arguments
    m <- params <- watchlist <- callbacks <- model <- train.missing <- model <- NULL
    return(res)
  }#End xgb_imp
  
  boost.res <- purrr::map_df(
    .x = markers.list, .f = xgb_imp,
    data.xgb = data.xgb, gl.wide = gl.wide,
    eta = eta, gamma = gamma,
    max_depth = max_depth, min_child_weight = min_child_weight,
    subsample = subsample, colsample_bytree = colsample_bytree,
    num_parallel_tree = num_parallel_tree,
    nthread = cpu.boost, nrounds = nrounds,
    early_stopping_rounds = early_stopping_rounds,
    save_name = save_name)
  
  return(boost.res)
  
}#End xgboost


# grur_lgbm_imputer ------------------------------------------------------------
#' @title grur_lgbm_imputer
#' @description imputations using Ligh Gradient Boosting Machine
#' @rdname grur_lgbm_imputer
#' @keywords internal
#' @export

grur_lgbm_imputer <- function(
  boost.split = NULL,
  markers.df = NULL,
  data.gbm = NULL,
  boosting = "dart",
  objective = "multiclass",
  learning_rate = 0.1,
  feature_fraction = 0.9,
  bagging_fraction = 0.9,
  bagging_freq = 1,
  max_depth = -1,
  min_data_in_leaf = 20,#default
  num_leaves = 31,#default
  cpu.boost = parallel::detectCores() / 2,
  early_stopping_rounds = 10,
  nrounds = 200,
  iteration.subsample = 2,
  pmm = 2
) {
  timing.imp <- proc.time() #for timing
  # boost.split <- 1 #test
  # data.gbm <- data.boost #test
  markers.list <- dplyr::filter(markers.df, SPLIT_VEC == boost.split) %>%
    dplyr::select(MARKERS) %>% purrr::flatten_chr(.)
  
  lightbgm_imp <- function(
    markers.list = NULL,
    data.gbm = NULL,
    boosting = "dart",
    objective = "multiclass",
    learning_rate = 0.1,
    feature_fraction = 0.9,
    bagging_fraction = 0.9,
    bagging_freq = 1,
    max_depth = -1,
    min_data_in_leaf = 20,#default
    num_leaves = 31,#default
    cpu.boost = parallel::detectCores() / 2,
    early_stopping_rounds = 10,
    nrounds = 200,
    iteration.subsample = 2,
    pmm = 2
  ) {
    # message("Imputations markers: ", markers.list)
    # preparing data
    # markers.list <- "BINDED_M667_M668"
    all.var<- colnames(data.gbm)
    data.var <- !all.var %in% c(markers.list)
    data.na <- is.na(data.gbm[, all.var, drop = FALSE])
    select.na <- data.na[, markers.list]
    
    # use train.gbm but here it's just to recycle object name below
    # here it's all the data
    xtrain <- data.gbm <- data.gbm[!select.na, data.var] #xtrain for PMM
    data.label <- data.gbm[!select.na, markers.list]
    ytrain <- unname(data.label)#ytrain for PMM
    
    # Note to myself: sampling for training and test set was not a good choice
    # with genomic data. Some training sets didn't contain all the potential genotypes.
    # tried sample with prob but trying a stratified sampling below
    
    # before
    # sample 90% of rows (samples) for training
    # train.row <- sort(sample(x = 1:nrow(train.gbm), size = ceiling(0.9 * nrow(train.gbm)), replace = FALSE))
    # test.row <- purrr::discard(.x = 1:nrow(train.gbm), .p = 1:nrow(train.gbm) %in% train.row)
    
    subsampling_gbm <- function(
      iteration.subsample, data.label, data.gbm,
      boosting = "dart",
      objective = "multiclass",
      learning_rate = 0.1,
      feature_fraction = 0.9,
      bagging_fraction = 0.9,
      bagging_freq = 1,
      max_depth = -1,
      min_data_in_leaf = 20,#default
      num_leaves = 31,#default
      cpu.boost = parallel::detectCores() / 2,
      early_stopping_rounds = 10,
      nrounds = 200) {
      # now using this function to make sure that label as all the potential num_classes
      stratified_sampling <- function(x, size = 0.9) {
        # x <- train.label
        # size <-  0.9
        data <- tibble::data_frame(LABEL = x) %>%
          dplyr::mutate(ROWS = seq(1, n(), by = 1))
        
        sampled <- data %>% 
          dplyr::group_by(LABEL) %>% 
          dplyr::sample_frac(tbl = ., size = size) %>% 
          dplyr::arrange(ROWS)
        
        train.row <- sampled$ROWS
        train.label <- sampled$LABEL
        test.label <- data %>% dplyr::filter(!ROWS %in% sampled$ROWS) %>% 
          dplyr::select(LABEL) %>% purrr::flatten_dbl(.)
        return(res = list(
          train.row = train.row, 
          train.label = train.label,
          test.label = test.label))
      }#End stratified_sampling
      
      sampled <- stratified_sampling(data.label)
      
      train.row <- sampled$train.row
      test.gbm <- data.gbm[-train.row,]
      test.label <- sampled$test.label
      train.gbm <- data.gbm[train.row,]
      train.label <- sampled$train.label
      sampled <- NULL
      # train.gbm <- lightgbm::lgb.Dataset(data = train.gbm, label = train.label)
      
      # test.gbm <- lightgbm::lgb.Dataset.create.valid(dataset = train.gbm,
      # data = test.gbm, label = test.label)
      valids <- list(test = test.gbm)
      
      # parameters
      params <- list(
        boosting = boosting, #"gbdt"#"dart",
        # boosting = "gbdt",
        objective = objective,
        num_classes = length(unique(train.label)),
        learning_rate = learning_rate,
        feature_fraction = feature_fraction,
        bagging_fraction = bagging_fraction,
        bagging_freq = bagging_freq,
        max_depth = max_depth,
        min_data_in_leaf = min_data_in_leaf,
        num_leaves = num_leaves,
        device = "cpu",
        num_threads = cpu.boost,
        early_stopping_rounds = early_stopping_rounds,
        # metric = "multi_logloss")
        metric = "multi_error")
      
      
      # model <- lightgbm::lgb.train(
      #   # init_model = model,
      #   params = params,
      #   data = train.gbm,
      #   nrounds = nrounds,
      #   valids = valids,
      #   pred_early_stop = TRUE,
      #   verbose = -1,
      #   reset_data = TRUE
      # )
      
      model.df <- tibble::data_frame(ITERATIONS = iteration.subsample, MODEL = list(model), SCORE = model$best_score)
      gc(verbose = FALSE)
      return(model.df)
    } # End of loop
    
    # #testing
    # learning_rate <- 0.05
    # feature_fraction <- 0.5
    # bagging_fraction <- 0.5
    # bagging_freq <- 1
    # min_data_in_leaf <- 20
    # num_leaves <- 31
    # max_depth <- -1
    # early_stopping_rounds <- 20
    # iteration.subsample <- 20
    
    model <- purrr::map_df(
      .x = 1:iteration.subsample,
      .f = subsampling_gbm,
      data.label = data.label, 
      data.gbm = data.gbm,
      boosting = boosting,
      objective = objective,
      learning_rate = learning_rate,
      feature_fraction = feature_fraction,
      bagging_fraction = bagging_fraction,
      bagging_freq = bagging_freq,
      max_depth = max_depth,
      min_data_in_leaf = min_data_in_leaf,
      num_leaves = num_leaves,
      cpu.boost = cpu.boost,
      early_stopping_rounds = early_stopping_rounds,
      nrounds = nrounds
    ) %>% 
      dplyr::filter(SCORE == min(SCORE)) %>% 
      dplyr::distinct(SCORE, .keep_all = TRUE) %>% 
      dplyr::select(MODEL) %>% 
      purrr::flatten(.)
    model <- model$MODEL
    error <- model$best_score
    
    train.missing <- data.gbm[select.na, data.var, drop = FALSE]
    id.string <- train.missing[, "INDIVIDUALS"]
    
    all.var <- data.na <- train.gbm <- data.var <- select.na <- data.gbm <- data.gbm <- NULL
    
    res <- model$predict(train.missing, reshape = TRUE) %>% #using best_iter by default
      tibble::as_data_frame(.) %>% 
      `colnames<-`(sort(unique(data.label))) %>%
      dplyr::mutate(
        INDIVIDUALS = id.string,
        MARKERS = rep(markers.list, n())) %>% 
      tidyr::gather(data = ., key = GT, value = SCORE, -c(INDIVIDUALS, MARKERS)) %>%
      dplyr::group_by(INDIVIDUALS, MARKERS) %>% 
      dplyr::filter(SCORE == max(SCORE)) %>%
      dplyr::distinct(INDIVIDUALS, MARKERS, .keep_all = TRUE) %>% 
      dplyr::arrange(INDIVIDUALS) %>%
      dplyr::mutate(GT = as.numeric(GT)) %>% 
      dplyr::select(-SCORE) %>% dplyr::ungroup(.)
    
    train.missing <- NULL
    
    # Predictive Mean Matching
    # pmm = 2 # test
    if (pmm > 0) {
      xtrain <- tibble::as_data_frame(model$predict(xtrain, reshape = TRUE)) %>% 
        `colnames<-`(sort(unique(data.label))) %>%
        dplyr::mutate(INDIVIDUALS = unique(xtrain[, "INDIVIDUALS"])) %>% 
        tidyr::gather(data = ., key = GT, value = SCORE, -INDIVIDUALS) %>%
        dplyr::group_by(INDIVIDUALS) %>% 
        dplyr::filter(SCORE == max(SCORE)) %>%
        dplyr::distinct(INDIVIDUALS, .keep_all = TRUE) %>% 
        dplyr::arrange(INDIVIDUALS) %>%
        dplyr::mutate(GT = as.numeric(GT)) %>% 
        dplyr::select(-SCORE) %>%
        dplyr::ungroup(.) %>% 
        dplyr::select(GT) %>% 
        purrr::flatten_dbl(.)
      res$GT <- missRanger::pmm(xtrain, xtest = res$GT, ytrain, k = pmm, seed = NULL)
    }
    res <- tibble::data_frame(IMPUTED_DATA = list(res), MARKERS = markers.list, ERROR = error)
    return(res)
  }# End lightbgm_imp
  
  #testing
  # learning_rate <- 0.1
  # feature_fraction <- 0.9
  # bagging_fraction <- 0.9
  # bagging_freq <- 1
  # min_data_in_leaf <- 20
  # num_leaves <- 31
  # max_depth <- -1
  # early_stopping_rounds <- 20
  # iteration.subsample <- 2
  # pmm <- 2
  
  res <- purrr::map_df(
    .x = markers.list, .f = lightbgm_imp,
    data.gbm = data.gbm,
    boosting = boosting, objective = objective,
    learning_rate = learning_rate,
    feature_fraction = feature_fraction,
    bagging_fraction = bagging_fraction,
    max_depth = max_depth,
    min_data_in_leaf = min_data_in_leaf,
    num_leaves = num_leaves,
    cpu.boost = cpu.boost,
    early_stopping_rounds = early_stopping_rounds,
    nrounds = nrounds,
    iteration.subsample = 2, pmm = pmm)
  timing.imp <- proc.time() - timing.imp
  message("    Imputations round ", boost.split, "/10 conducted in ", round(timing.imp[[3]]), " sec")
  return(res)
}#End grur_lgbm_imputer

# encoding_snp --------------------------------------------------------------------
#' @title encoding_snp
#' @description bind snp found on the same locus
#' @rdname encoding_snp
#' @keywords internal
#' @export
encoding_snp <- function(locus.list = NULL, data = NULL) {
  # locus.list <- "1_135"
  res <- dplyr::filter(.data = data, CHROM_LOCUS %in% locus.list)
  binded.markers <- dplyr::distinct(.data = res, MARKERS) %>%
    purrr::flatten_chr(.) %>%
    stringi::stri_join(., collapse = "_")
  binded.markers <- stringi::stri_join("BINDED_", binded.markers)
  
  res <- res %>%
    dplyr::group_by(CHROM_LOCUS, STRATA, INDIVIDUALS) %>%
    tidyr::spread(data = ., key = MARKERS, value = GT) %>%
    tidyr::unite(data = ., col = GT, -CHROM_LOCUS, -STRATA, -INDIVIDUALS, sep = "_", remove = TRUE) %>%
    dplyr::mutate(#Haplotype with combination of SNP and NA = NA (up for an argument?)
      GT = stringi::stri_replace_all_fixed(
        str = GT, pattern = "NA", replacement = NA, vectorize_all = FALSE),
      MARKERS = rep(binded.markers, n())
    ) %>%
    dplyr::ungroup(.) %>%
    dplyr::select(MARKERS, CHROM_LOCUS, STRATA, INDIVIDUALS, GT)
  
  return(res)
}#End encoding_snp

# decoding_haplotypes------------------------------------------------------------
#' @title decoding_haplotypes
#' @description separate snp group merged with encoding_snp
#' @rdname decoding_haplotypes
#' @keywords internal
#' @export

decoding_haplotypes <- function(data = NULL, parallel.core = parallel::detectCores() - 1) {
  # data <- data.imp.bk#test
  # data <- data.imp#test
  
  # find the markers that need sep.
  if (tibble::has_name(data, "MARKERS")) {
    markers.sep <- dplyr::distinct(data, MARKERS) %>%
      dplyr::filter(stringi::stri_detect_regex(str = MARKERS, pattern = "^BINDED")) %>%
      purrr::flatten_chr(.)
  } else {
    col.names.data <- colnames(data)
    markers.sep <- purrr::keep(
      .x = col.names.data,
      .p = stringi::stri_detect_regex(str = col.names.data, pattern = "^BINDED"))
    col.names.data <- NULL
  }
  # nested function required ---------------------------------------------------
  separate_locus <- function(binded.markers = NULL, data = NULL) {
    # binded.markers <- markers.sep[[1]]
    
    col.replace <- stringi::stri_replace_all_fixed(
      str = binded.markers, pattern = "BINDED_", replacement = "", vectorize_all = FALSE)
    
    if (tibble::has_name(data, "MARKERS")) {
      data.sep <- dplyr::filter(.data = data, MARKERS %in% binded.markers) %>%
        dplyr::select(-MARKERS) %>%
        tidyr::separate_(
          data = .,
          col = "GT",
          into = stringi::stri_split_fixed(str = col.replace, pattern = "_", simplify = TRUE),
          sep = "_", extra = "drop")
      
      if (tibble::has_name(data.sep, "GL")) {
        data.sep <- tidyr::gather(data = data.sep, key = MARKERS, value = GT, -STRATA, -INDIVIDUALS, -GL) %>%
          dplyr::select(STRATA, INDIVIDUALS, MARKERS, GT, GL) %>%
          dplyr::arrange(STRATA, INDIVIDUALS, MARKERS, GT, GL)
      } else {
        data.sep <- tidyr::gather(data = data.sep, key = MARKERS, value = GT, -STRATA, -INDIVIDUALS) %>%
          dplyr::select(STRATA, INDIVIDUALS, MARKERS, GT) %>%
          dplyr::arrange(STRATA, INDIVIDUALS, MARKERS, GT)
      }
    } else {
      data.sep <- dplyr::select(.data = data, dplyr::one_of(c("STRATA", "INDIVIDUALS", binded.markers)))
      colnames(data.sep) <- c("STRATA", "INDIVIDUALS", col.replace)
      data.sep <- tidyr::separate_(
        data = data.sep,
        col = col.replace,
        into = stringi::stri_split_fixed(str = col.replace, pattern = "_", simplify = TRUE),
        sep = "_", extra = "drop") %>%
        tidyr::gather(data = ., key = MARKERS, value = GT, -STRATA, -INDIVIDUALS)
    }
    
    return(data.sep)
  }#End separate_locus
  
  if (length(markers.sep) > 0) {
    if (length(markers.sep) > 100) {
      data.sep <- list()
      data.sep <- .grur_parallel_mc(
        X = markers.sep,
        FUN = separate_locus,
        mc.cores = parallel.core,
        data = data
      ) %>% dplyr::bind_rows(.)
    } else {
      data.sep <- purrr::map(.x = markers.sep, .f = separate_locus, data = data) %>%
        dplyr::bind_rows(.)
    }
    
    # Include markers 1 snp/read
    if (tibble::has_name(data, "MARKERS")) {
      markers.no.sep <- dplyr::distinct(data, MARKERS) %>%
        dplyr::filter(!stringi::stri_detect_regex(str = MARKERS, pattern = "^BINDED")) %>%
        purrr::flatten_chr(.)
    } else {
      col.names.data <- colnames(data)
      markers.no.sep <- purrr::discard(
        .x = col.names.data,
        .p = stringi::stri_detect_regex(str = col.names.data, pattern = "^BINDED"))
      markers.no.sep <- purrr::discard(
        .x = markers.no.sep,
        .p = markers.no.sep %in% c("STRATA", "INDIVIDUALS"))
      col.names.data <- NULL
    }
    
    if (length(markers.no.sep) > 0) {
      if (tibble::has_name(data, "MARKERS")) {
        data.no.sep <- suppressWarnings(
          dplyr::filter(.data = data, MARKERS %in% markers.no.sep) %>%
            dplyr::select(dplyr::one_of(c("STRATA", "INDIVIDUALS", "MARKERS", "GT", "GL")))
        )
      } else {
        data.no.sep <- suppressWarnings(
          dplyr::select(
            .data = data,
            dplyr::one_of(c("STRATA", "INDIVIDUALS", markers.no.sep))) %>%
            tidyr::gather(data = ., key = MARKERS, value = GT, -STRATA, -INDIVIDUALS)
        )
      }
    } else {
      data.no.sep <- NULL
    }
    
    # combined data
    if (!is.null(data.no.sep)) {
      data.sep <- dplyr::bind_rows(data.sep, data.no.sep) %>%
        dplyr::arrange(STRATA, INDIVIDUALS, MARKERS)
    }
  } else {
    data.sep <- data
  }
  return(data.sep)
}#End decoding_haplotypes

# factorize_gt------------------------------------------------------------
#' @title factorize_gt
#' @description Necessary to factorize by markers in tidy format.
#' XGBoost needs numbering to start at , hence the codes below.
#' @rdname factorize_gt
#' @keywords internal
#' @export

factorize_gt <- function(x) {
  x <- as.numeric(factor(x)) - 1
}#End factorize_gt


# defactorize_gt------------------------------------------------------------
#' @title defactorize_gt
#' @description Function to "defactorize/decode" the imputed data back to original.
#' @rdname defactorize_gt
#' @keywords internal
#' @export

defactorize_gt <- function(data.to.change,
                           data.with.info = "imputation_factor_dictionary.rad") {
  #data.with.info <- fst::read.fst(path = data.with.info)
  data.with.info <- readr::read_tsv(data.with.info)
  
  clean.id <- dplyr::distinct(.data = data.with.info, INDIVIDUALS, INDIVIDUALS_N)
  clean.gt <- dplyr::distinct(.data = data.with.info, MARKERS, GT, GT_N) %>%
    tidyr::drop_na(.)
  
  if (tibble::has_name(data.to.change, "STRATA")) {
    clean.pop <- dplyr::distinct(.data = data.with.info, STRATA, POP_ID_N)
    res <- suppressWarnings(
      dplyr::arrange(.data = data.to.change, STRATA, INDIVIDUALS) %>%
        dplyr::rename(POP_ID_N = STRATA, INDIVIDUALS_N = INDIVIDUALS, GT_N = GT) %>%
        dplyr::inner_join(clean.id, by = "INDIVIDUALS_N") %>%
        dplyr::select(-INDIVIDUALS_N) %>%
        dplyr::inner_join(clean.pop, by = "POP_ID_N") %>%
        dplyr::select(-POP_ID_N) %>%
        dplyr::left_join(clean.gt, by = c("MARKERS", "GT_N")) %>%
        dplyr::select(STRATA, INDIVIDUALS, MARKERS, GT) %>%
        dplyr::bind_rows(
          tidyr::drop_na(
            data = dplyr::select(
              .data = data.with.info,
              dplyr::one_of(c("STRATA", "INDIVIDUALS", "MARKERS", "GT", "GL"))
            ))))
  } else {
    res <- suppressWarnings(
      dplyr::arrange(.data = data.to.change, INDIVIDUALS) %>%
        dplyr::rename(INDIVIDUALS_N = INDIVIDUALS, GT_N = GT) %>%
        dplyr::inner_join(clean.id, by = "INDIVIDUALS_N") %>%
        dplyr::select(-INDIVIDUALS_N) %>%
        dplyr::left_join(clean.gt, by = c("MARKERS", "GT_N")) %>%
        dplyr::select(INDIVIDUALS, MARKERS, GT) %>%
        dplyr::bind_rows(
          tidyr::drop_na(
            data = dplyr::select(
              .data = data.with.info,
              dplyr::one_of(c("INDIVIDUALS", "MARKERS", "GT", "GL"))
            ))))
  }
  return(res)
}#End defactorize_gt


#' @title split_imp
#' @description Split markers into chunk for parallel imputations processing
#' @rdname split_imp
#' @keywords internal
#' @export
split_imp <- function(x, cpu.rounds, parallel.core = parallel::detectCores() - 1) {
  n.row <- nrow(x)
  split.vec <- as.integer(floor((parallel.core * cpu.rounds * (1:n.row - 1) / n.row) + 1))
  return(split.vec)
}#End split_imp
