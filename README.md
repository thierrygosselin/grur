[![Travis-CI Build Status](https://travis-ci.org/thierrygosselin/grur.svg?branch=master)](https://travis-ci.org/thierrygosselin/grur) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/thierrygosselin/grur?branch=master&svg=true)](https://ci.appveyor.com/project/thierrygosselin/grur) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/grur)](http://cran.r-project.org/package=grur) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![DOI](https://zenodo.org/badge/14548/thierrygosselin/grur.svg)](https://zenodo.org/badge/latestdoi/14548/thierrygosselin/grur)

[![packageversion](https://img.shields.io/badge/Package%20version-0.0.1-orange.svg)](commits/master) [![Last-changedate](https://img.shields.io/badge/last%20change-2017--04--08-brightgreen.svg)](/commits/master)

------------------------------------------------------------------------

grur: an R package tailored for RADseq data imputation
======================================================

Don't get trapped by missing genotypes in your RADseq dataset. Numerous genomic analysis are vulnerable to missing values, so you can use **grur** to:

-   **Visualize patterns of missing genotypes:** find patterns associated with different variables of your study (lanes, chips, sequencers, populations, sample sites, reads/samples, homozygosity). The analysis is similar to PLINK's identify-by-missingness analysis (IBM), but more powerful (more output and figures) and as fast!

-   **Perform map-independent imputations of missing genotypes:** using **Random Forests** (on-the-fly-imputations or using predictive modeling), **Extreme Gradient Tree Boosting** and the classic Strawman imputation that uses the most frequently observed, non-missing genotypes. The algorithms can build the model for the imputations using the **population** groupings or not (i.e. **overall samples**).

Installation
------------

To try out the dev version of **grur**, copy/paste the code below:

``` r
if (!require("devtools")) install.packages("devtools") # to install
devtools::install_github("thierrygosselin/grur")
library(grur) # to load
```

<table style="width:100%;">
<colgroup>
<col width="26%" />
<col width="73%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Caracteristics</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><strong>Input/Output</strong></td>
<td align="left">The imputations offered in <strong>grur</strong> are seamlesly integrated in <strong>stackr</strong> and **assigner. Imputations are also integrated with usefull filters, blacklist and whitelist inside those 2 packages. Format includes: <a href="https://samtools.github.io/hts-specs/">VCF, SNPs and haplotypes</a>, <a href="http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr">PLINK tped/tfam</a>, <a href="https://github.com/thibautjombart/adegenet">genind</a>, <a href="https://github.com/thibautjombart/adegenet">genlight</a>, <a href="https://github.com/EricArcher/strataG">strataG gtypes</a>, <a href="http://genepop.curtin.edu.au">Genepop</a>, <a href="http://catchenlab.life.illinois.edu/stacks/">STACKS haplotype file</a>, <a href="https://github.com/jgx65/hierfstat">hierfstat</a>, <a href="https://www.zsl.org/science/software/colony">COLONY</a>, <a href="http://adn.biol.umontreal.ca/~numericalecology/Rcode/">betadiv</a>, <a href="http://gutengroup.mcb.arizona.edu/software/">δaδi</a>, <a href="http://pritchardlab.stanford.edu/structure.html">structure</a>, <a href="http://cmpg.unibe.ch/software/arlequin35/">Arlequin</a>, <a href="https://github.com/zhengxwen/SNPRelate">SNPRelate</a>, Dataframes of genotypes in wide or long/tidy format</td>
</tr>
<tr class="even">
<td align="left"><strong>Pattern of missingness</strong></td>
<td align="left"><code>missing_visualization</code>: Visualize patterns of missing data. Find patterns associated with different variables of your study (lanes, chips, sequencers, populations, sample sites, reads/samples, homozygosity, etc)</td>
</tr>
<tr class="odd">
<td align="left"><strong>Imputations</strong></td>
<td align="left"><strong>Map-independent</strong> imputations of missing genotypes.<br>Using <strong>Random Forests</strong> (on-the-fly-imputations or predictive modeling), <strong>Extreme Gradient Tree Boosting</strong> and Strawman imputations (~ max/mean/mode: the most frequently observed, non-missing genotypes is used).<br> Imputations can be conducted <strong>overall samples</strong> or <strong>by populations</strong>.<br><br>Imputations are integrated in several of <strong>stackr</strong> and <strong>assigner</strong> functions.</td>
</tr>
<tr class="even">
<td align="left"><strong><a href="http://ggplot2.org">ggplot2</a>-based plotting</strong></td>
<td align="left">Visualize distribution of important metric and statistics and create publication-ready figures.</td>
</tr>
<tr class="odd">
<td align="left"><strong>Parallel</strong></td>
<td align="left">Codes designed and optimized for fast computations running imputations, iterations, etc. in parallel. Works with all OS: Linux, Mac and now PC!</td>
</tr>
</tbody>
</table>

[More in grur workflow below](https://github.com/thierrygosselin/grur#grur-workflow)

Prerequisite - Suggestions - Troubleshooting
--------------------------------------------

-   **Parallel computing**: Follow the steps in this [vignette](https://github.com/thierrygosselin/grur/blob/master/vignettes/vignette_imputations_parallel.Rmd) to install [data.table](https://github.com/Rdatatable/data.table) and [XGBoost](https://github.com/dmlc/xgboost) packages (e.g. to do imputations in parallel).
-   **Installation problem:** see this [vignette](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_installation_problems.Rmd)
-   **Windows users**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
-   For a better experience in **stackr** and in R in general, I recommend using [RStudio](https://www.rstudio.com/products/rstudio/download/). The R GUI is unstable with functions using parallel ([more info](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/mclapply.html)).

Vignettes, R Notebooks and examples
-----------------------------------

**Vignettes (in development, check periodically for updates):**

-   Vignettes with real data for example in the form of R Notebooks take too much space to be included in package, without CRAN complaining. Consequently, vignettes are distributed separately, follow the links below.
-   [installation problems](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_installation_problems.Rmd)
-   [parallel computing during imputations](https://github.com/thierrygosselin/grur/blob/master/vignettes/vignette_imputations_parallel.Rmd)
-   Missing data visualization and analysis [(html vignette)](https://www.dropbox.com/s/4zf032g6yjatj0a/vignette_missing_data_analysis.nb.html?dl=0) and [(Rmd)](https://www.dropbox.com/s/5fxw2h9w1l1j391/vignette_missing_data_analysis.Rmd?dl=0)

Citation:
---------

To get the citation, inside R:

``` r
citation("grur")
```

New features
------------

Change log, version, new features and bug history lives in the [NEWS.md file](https://github.com/thierrygosselin/grur/blob/master/NEWS.md)

**v.0.0.1 2017-04-07**

-   **grur** package launch!

Roadmap of future developments:
-------------------------------

-   Integrate more imputation method.
-   Workflow tutorial and vignettes to further explore some problems: *in progress*
-   Use Shiny and ggvis (when subplots and/or facets becomes available for ggvis).
-   Until publication **grur** will change rapidly, stay updated with releases and contribute with bug reports.
-   Suggestions ?

Contributions:
--------------

This package has been developed in the open, and it wouldn’t be nearly as good without your contributions. There are a number of ways you can help me make this package even better:

-   If you don’t understand something, please let me know.
-   Your feedback on what is confusing or hard to understand is valuable.
-   If you spot a typo, feel free to edit the underlying page and send a pull request.

New to pull request on github ? The process is very easy:

-   Click the edit this page on the sidebar.
-   Make the changes using github’s in-page editor and save.
-   Submit a pull request and include a brief description of your changes.
-   “Fixing typos” is perfectly adequate.

GBS workflow
------------

The **grur** package should fit here in your RADseq/GBS workflow.

stackr workflow
---------------

Currently under construction. Come back soon!

**Pattern of missingness**

-   Use `missing_visualization` with/without your new blacklists (e.g. of genotypes, individuals) and with/without whitelist of markers to examine patterns of missingness in you dataset before more extensive filtering (there is a vignette for this step)
-   The trick here is to use the `strata` argument to find patterns associated with different variables of your study (lanes, chips, sequencers, populations, sample sites, reads/samples, etc).
-   Do you see a trend between your missing pattern and reads/samples ? Heterozygosity?
-   Do you need more sequencing? Do you have to re-run some lanes?
-   Use imputation methods provided inside some of **stackr** functions (e.g. `tidy_genomic_data` or `genomic_converter`), to assess the impact of lowering or increasing different filtering thresholds that impact missing data.
-   Use `missing_visualization` with your new blacklists (e.g. of genotypes, individuals) and with your whitelist of markers to examine patterns of missingness in your dataset after filtering (there is a vignette for this step)
-   The trick again here is to use the `strata` argument to find patterns associated with different variables of your study (lanes, chips, sequencers, populations, sample sites, reads/samples, etc).
-   Do you see a trend between your missing pattern and reads/samples ? Heterozygosity?
