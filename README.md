
<!-- README.md is generated from README.Rmd. Please edit that file -->

# grur<img src="README_grur_logo.png" align="right"/>

[![Travis-CI Build
Status](https://travis-ci.org/thierrygosselin/grur.svg?branch=master)](https://travis-ci.org/thierrygosselin/grur)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/thierrygosselin/grur?branch=master&svg=true)](https://ci.appveyor.com/project/thierrygosselin/grur)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/grur)](http://cran.r-project.org/package=grur)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/87596763.svg)](https://zenodo.org/badge/latestdoi/87596763)

[![packageversion](https://img.shields.io/badge/Package%20version-0.0.11-orange.svg)](commits/master)
[![Last-changedate](https://img.shields.io/badge/last%20change-2018--07--09-brightgreen.svg)](/commits/master)

-----

<https://thierrygosselin.github.io/grur/>

The name **grur** |ɡro͞oˈr| was chosen because the missing genotypes
dilemma with RADseq data reminds me of the cheese paradox.

Here, I don’t want to sustain a
[war](http://www.lefigaro.fr/flash-eco/2012/12/07/97002-20121207FILWWW00487-le-gruyere-francais-doit-avoir-des-trous.php)
or the controversy of cheese with holes, so choose as you like, the
*French Gruyère* or the *Swiss Emmental*. The paradox is that the more
cheese you have the more holes you’ll get. But, the more holes you have
means the less cheese you have… So, someone could conclude, the more
cheese = the less cheese ? I’ll leave that up to you, back to genomics…

Numerous genomic analysis are vulnerable to missing values, don’t get
trapped by missing genotypes in your RADseq dataset.

Use **grur** to **visualize patterns of missingness** and **perform
map-independent imputations of missing genotypes** (see
[features](https://github.com/thierrygosselin/grur#features) below).

## Installation

To try out the dev version of **grur**, copy/paste the code below:

``` r
if (!require("devtools")) install.packages("devtools") # to install
devtools::install_github("thierrygosselin/grur")
library(grur)
```

## Options and required packages

Please follow instructions in the [Notebook
vignette](https://www.dropbox.com/s/5npumwdo0cxtxi4/rad_genomics_computer_setup.nb.html?dl=0)
to install required packages for the selected imputation options
below:

| imputation options                 |      package      | install instructions                                                                                    |
| :--------------------------------- | :---------------: | ------------------------------------------------------------------------------------------------------- |
| **imputation.method = “lightgbm”** |    `lightgbm`     | [Notebook vignette](https://www.dropbox.com/s/5npumwdo0cxtxi4/rad_genomics_computer_setup.nb.html?dl=0) |
| **imputation.method = “xgboost”**  |     `xgboost`     | [Notebook vignette](https://www.dropbox.com/s/5npumwdo0cxtxi4/rad_genomics_computer_setup.nb.html?dl=0) |
| **imputation.method = “rf”**       | `randomForestSRC` | [Notebook vignette](https://www.dropbox.com/s/5npumwdo0cxtxi4/rad_genomics_computer_setup.nb.html?dl=0) |
| **imputation.method = “rf\_pred”** |     `ranger`      | `install.packages("ranger")`                                                                            |
| **imputation.method = “bpca”**     |   `pcaMethods`    | `source("https://bioconductor.org/biocLite.R")`<br>`biocLite("pcaMethods")`                             |
| **if using pmm \> 0**              |   `missRanger`    | `install.packages("missRanger")`                                                                        |

#### Troubleshooting

  - `rmetasim` needs to be modified in order to simulate more than 2000
    markers [notebook
    vignette](https://www.dropbox.com/s/5npumwdo0cxtxi4/rad_genomics_computer_setup.nb.html?dl=0)
  - **Parallel computing**: follow the steps in this [notebook
    vignette](https://www.dropbox.com/s/5npumwdo0cxtxi4/rad_genomics_computer_setup.nb.html?dl=0)
    to install the packages with OpenMP-enabled compiler and conduct
    imputations in parallel.
  - [Installation
    problems.](https://www.dropbox.com/s/5npumwdo0cxtxi4/rad_genomics_computer_setup.nb.html?dl=0)
  - **Windows users**: Install
    [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
  - The R GUI is unstable with functions using parallel ([more
    info](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/mclapply.html)),
    so I recommend using
    [RStudio](https://www.rstudio.com/products/rstudio/download/) for a
    better experience.
  - Running codes in chunks inside R Notebook might cause problem, run
    it outside in the
console.

## Features

| Caracteristics                                   | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| :----------------------------------------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Simulate RADseq data**                         | `simulate_rad`: simulate populations of RADseq data following island or stepping stone models. Inside the function, allele frequency can be created with [fastsimcoal2](http://cmpg.unibe.ch/software/fastsimcoal2/) and then used inside [rmetasim](https://github.com/stranda/rmetasim) simulation engine. *Vignette coming soon*.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| **Patterns of missingness**                      | `missing_visualization`: visualize patterns of missing data associated with different variables of your study (lanes, chips, sequencers, populations, sample sites, reads/samples, homozygosity, etc). Similar to PLINK’s identify-by-missingness analysis (IBM), **grur** is more powerful because it generates more analysis and automatically creates tables and figures. Vignette: [html](https://www.dropbox.com/s/4zf032g6yjatj0a/vignette_missing_data_analysis.nb.html?dl=0) and [Rmd](https://www.dropbox.com/s/5fxw2h9w1l1j391/vignette_missing_data_analysis.Rmd?dl=0)<br><br>`generate_missing`: allows to generate missing genotypes in dateset \[simulated\] based on a compound Dirichlet-multinomial distribution. *Vignette coming soon*.                                                                                                                                                                                                                                                                                                                                                                              |
| **Imputations**                                  | `grur_imputations`: **Map-independent** imputations of missing genotypes with several algorithms (including machine leaning):<br> \* **Random Forests** (on-the-fly-imputations with randomForestSRC or using predictive modelling using ranger and missRanger),<br>\* **Extreme Gradient Tree Boosting** (using XGBoost or LightGBM),<br>\* **Bayesian PCA** (using bpca in pcaMethods),<br>\* **Classic Strawman: ** the most frequently observed, non-missing, genotypes is used for imputation.<br><br>**Hierarchy: ** algorithm’s model can account for *strata* groupings, e.g. if patterns of missingness is found in the data.<br><br>**Haplotypes: ** automatically detect SNPs on the same LOCUS (read/haplotype) to impute the SNPs jointly, reducing imputation artifacts. *Vignette coming soon*.                                                                                                                                                                                                                                                                                                                          |
| **Input/Output**                                 | The imputations offered in **grur** are seamlesly integrated in **radiator** and **assigner**. Imputations are also integrated with usefull filters, blacklists and whitelists inside those 2 packages. Genetic formats recognized: [VCF, SNPs and haplotypes](https://samtools.github.io/hts-specs/), [PLINK tped/tfam](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr), [genind](https://github.com/thibautjombart/adegenet), [genlight](https://github.com/thibautjombart/adegenet), [strataG gtypes](https://github.com/EricArcher/strataG), [Genepop](http://genepop.curtin.edu.au), [STACKS haplotype file](http://catchenlab.life.illinois.edu/stacks/), [hierfstat](https://github.com/jgx65/hierfstat), [COLONY](https://www.zsl.org/science/software/colony), [betadiv](http://adn.biol.umontreal.ca/~numericalecology/Rcode/), [δaδi](http://gutengroup.mcb.arizona.edu/software/), [structure](http://pritchardlab.stanford.edu/structure.html), [Arlequin](http://cmpg.unibe.ch/software/arlequin35/), [SNPRelate](https://github.com/zhengxwen/SNPRelate), dataframes of genotypes in wide or long/tidy format. |
| **[ggplot2](http://ggplot2.org)-based plotting** | Visualization: publication-ready figures of important metrics and statistics.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| **Parallel**                                     | Codes designed and optimized for fast computations with progress bars. Works with all OS: Linux, Mac and yes PC\!                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |

## Vignettes and examples

Vignettes with real data for example in the form of R Notebooks take too
much space to be included in package, without CRAN complaining.
Consequently, vignettes will start to be distributed separately, follow
the links below.

  - Missing data visualization and analysis [(html
    vignette)](https://www.dropbox.com/s/4zf032g6yjatj0a/vignette_missing_data_analysis.nb.html?dl=0)
    and
    [(Rmd)](https://www.dropbox.com/s/5fxw2h9w1l1j391/vignette_missing_data_analysis.Rmd?dl=0)

## Citation

To get the citation, inside R:

``` r
citation("grur")
```

## New features

Change log, version, new features and bug history lives in the [NEWS.md
file](https://github.com/thierrygosselin/grur/blob/master/NEWS.md)

**grur v.0.0.10 2018-04-26**

`grur's` dependencies:

  - I transferred to `Suggests` section these packages: lightgbm,
    missRanger, randomForestSRC, ranger, rmarkdown, rmetasim, strataG,
    xgboost.
  - Functions thate requires specific package will now say so.
  - Reason: people only interested in `missing_visualization` don’t have
    to install all the required packages required for imputations or
    simulations.

`simulate_rad`: with the latest R release (3.5.0), Check now throw a new
note: **Note: next used in wrong context: no loop is visible at
simulate\_rad.R:189** I replaced `next` inside `sapply` with `while`.

**grur v.0.0.9 2017-10-27**

  - `lightGBM` option to conduct the imputations is fully functional

#### Roadmap of future developments

  - Integrate more imputation method.
  - Workflow tutorial to further explore some problems.
  - Use Shiny and ggvis (when subplots and/or facets becomes available
    for ggvis).
  - Until publication **grur** will change rapidly, stay updated with
    releases and contribute with bug reports.
  - Suggestions ?

#### Contributions

This package has been developed in the open, and it wouldn’t be nearly
as good without your contributions. There are a number of ways you can
help me make this package even better:

  - If you don’t understand something, please let me know and raise an
    [issue](https://github.com/thierrygosselin/grur/issues)
  - Your feedback on what is confusing or hard to understand is
    valuable.
  - If you spot a typo, feel free to edit the underlying page and send a
    pull request.
  - New to pull request on github ? [The process is very
    easy](http://r-pkgs.had.co.nz/git.html#git-pullreq).
