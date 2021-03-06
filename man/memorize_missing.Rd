% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/memorize_missing.R
\name{memorize_missing}
\alias{memorize_missing}
\title{Memorize missingness pattern and randomize attributes}
\usage{
memorize_missing(data, strata = NULL, randomize = NULL, filename = NULL)
}
\arguments{
\item{data}{14 options for input (\strong{diploid data only}): VCFs (SNPs or Haplotypes,
to make the vcf population ready),
plink (tped, bed), stacks haplotype file, genind (library(adegenet)),
genlight (library(adegenet)), gtypes (library(strataG)), genepop, DArT,
and a data frame in long/tidy or wide format. To verify that radiator detect
your file format use \code{\link[radiator]{detect_genomic_format}} (see example below).
Documented in \strong{Input genomic datasets} of \code{\link[radiator]{tidy_genomic_data}}.

\strong{DArT and VCF data}: \pkg{radiator} was not meant to generate alleles
and genotypes if you are using a VCF file with no genotype
(only genotype likelihood: GL or PL).
Neither is \pkg{radiator} able to magically generate a genind object
from a SilicoDArT dataset. Please look at the first few lines of your dataset
to understand it's limit before asking raditor to convert or filter your dataset.}

\item{strata}{(optional/required) Required for VCF and haplotypes files,
optional for the other formats supported.
See documentation of \code{\link[radiator]{tidy_genomic_data}} for more info.
Default: \code{strata = NULL}.}

\item{randomize}{(optional, string) To randomize the missingness of specific attributes.
Available options: \code{"markers", "populations", "individuals" and "overall"}.
All options can be selected in a string,
\code{randomize = c("markers", "populations", "individuals", "overall")}
Default: \code{randomize = NULL} will only keep the original missingness pattern.}

\item{filename}{(optional) The name of the file (extension not necessary)
written to the working directory and containing the missing info.
Default: \code{filename = NULL}, the missing info is in the global
environment only.

\code{grur} takes advantage of the lightweight and speedy file reading/writing
package \code{fst} (\emph{Lightning Fast Serialization of Data Frames for R})
to write the dataframe to the working directory. This file can be used inside
\code{generate_missing} (coming soon) function.}
}
\value{
A tidy dataframe in the global environment with columns: 
\code{POP_ID}, \code{INDIVIDUALS}, \code{MARKERS}, and in the subsequent
columns, the missingness info coded 0 for missing and 1 for genotyped.
Depending on the value chosen for the argument \code{randomize}, 
the columns are:
\itemize{
   \item \code{MISSING_ORIGINAL}: for the original missing pattern (always present)
   \item \code{MISSING_MARKERS_MIX}: for the missing pattern randomized by markers (optional) 
   \item \code{MISSING_POP_MIX}: for the missing pattern randomized by populations (optional) 
   \item \code{MISSING_INDIVIDUALS_MIX}: for the missing pattern randomized by individuals (optional) 
   \item \code{MISSING_OVERALL_MIX}: for the missing pattern randomized overall (optional)
}
}
\description{
Use this function to keep the pattern of missing data (0/1).
The pattern can be randomized based on dataset attributes/covariates.
This can be useful to generate missingness on simulated dataset with the
same number of individuals, populations and markers or to analyze the accuracy
of imputation algorithms. A vignette is under construction to leverage this
function.
}
\examples{
\dontrun{
missing.memory <- memorize_missing(
data = "batch_1.vcf",
strata = "population.map.strata.tsv", 
randomize = "populations", filename = "missing.memory.panda"
)
}
}
\author{
Thierry Gosselin \email{thierrygosselin@icloud.com}
}
