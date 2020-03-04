# grur 0.1.2 2020-02-13

* lots of works on `simulate_rad` and `missing_distributions` from Eric Archer.


# grur 0.1.1 2019-05-01

* pkgdown website
* better function doc
* vignettes
* work on dependencies
* work on travis CI



# grur 0.1.0 2019-03-15

* several commits will come with this in the next 2 weeks
* `missing_visualization`: is getting an upgrade. Will work with GDS file/object.
* `grur_imputations`: will be fully reviewed and updated to work with GDS.
* `radiator::genomic_converter`: is now integrated at the end of the imputation 
module, instead of the later. This will reduce duplication codes for me, and help pass
CRAN...
* New documentations.
* Several arguments are now deprecated or moved to the advance section. Read documentation.
* new life cycle section in `grur_imputations`
* Requirement to pass some test before imputations will be generated so that
user with bad data don't waist time with imputations!


# grur 0.0.11 2018-07-09

* grur ready for R 3.5.1 "Feather Spray" released on 2018/07/05
* grur updated to work with ggplot2 v.3.0.0


# grur 0.0.10 2018-04-26

`grur's` dependencies:

  * I transferred to `Suggests` section these packages: 
  lightgbm, missRanger, randomForestSRC, ranger, rmarkdown, rmetasim, strataG,
  xgboost.
  * Functions thate requires specific package will now say so.
  * Reason: people only interested in `missing_visualization` don't have to install
  all the required packages required for imputations or simulations.

`simulate_rad`: with the latest R release (3.5.0), Check now throw a new note:
**Note: next used in wrong context: no loop is visible at simulate_rad.R:189**
I replaced `next` inside `sapply` with `while`.


# grur 0.0.9 2017-10-27

* `lightGBM` option to conduct the imputations is fully functional



# grur 0.0.8 2017-10-24

* worked on `missing_visualization` to get better output figures with large number of pop
* worked on imputation algorithm to prep for the use of `lightGBM` package


# grur 0.0.7 2017-10-02

worked on `missing_visualization`
* added Eric Archer's work on detecting patterns
* figures are now combined under themes (e.g. all PCoA for the same strata are together, etc)

# grur 0.0.6 2017-08-18

* several typos fixed
* work on `missing_visualization`

# grur 0.0.5 2017-08-15

* restored progress bar when using parallel computing by installing the new dev
version of `pbmcapply` package.



# grur 0.0.5 2017-08-15

* restored progress bar when using parallel computing by installing the new dev
version of `pbmcapply` package.


# grur 0.0.4 2017-08-15

* bug fix: removed the progress bar when using parallel computing.
This is temporary, while waiting for a fix with `pbmcapply` package.


# grur 0.0.3 2017-06-17

* `simulate_rad`: a new function developed by Eric Archer to simulate populations
of RADseq data following **island** or **stepping** stone models.
* **grur** works with `dplyr v.0.7.0`

# grur 0.0.2 2017-04-23

* Work on `generate_missing` function.
* Still not as I wish. But getting there.
* 2 new functions: `imputations_accuracy` and `memorize_missing`. Check functions
documentations for more info.


# grur 0.0.1 2017-04-07

New package is out!


