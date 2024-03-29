% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapper.R
\name{dive_phe2mash}
\alias{dive_phe2mash}
\title{Wrapper to run mash given a phenotype data frame}
\usage{
dive_phe2mash(
  df,
  snp,
  type = "linear",
  svd = NULL,
  suffix = "",
  outputdir = ".",
  min.phe = 200,
  ncores = NA,
  save.plots = TRUE,
  thr.r2 = 0.2,
  thr.m = c("max", "sum"),
  num.strong = 1000,
  num.random = NA,
  scale.phe = TRUE,
  roll.size = 50,
  U.ed = NA,
  U.hyp = NA,
  verbose = TRUE
)
}
\arguments{
\item{df}{Dataframe containing phenotypes for mash where the first column is
'sample.ID', which should match values in the snp$fam$sample.ID column.}

\item{snp}{A "bigSNP" object; load with \code{snp_attach()}.}

\item{type}{Character string, or a character vector the length of the number
of phenotypes. Type of univarate regression to run for GWAS.
Options are "linear" or "logistic".}

\item{svd}{A "big_SVD" object; Optional covariance matrix to use for
population structure correction.}

\item{suffix}{Optional character vector to give saved files a unique search string/name.}

\item{outputdir}{Optional file path to save output files.}

\item{min.phe}{Integer. Minimum number of individuals phenotyped in order to
include that phenotype in GWAS. Default is 200. Use lower values with
caution.}

\item{ncores}{Optional integer to specify the number of cores to be used
for parallelization. You can specify this with bigparallelr::nb_cores().}

\item{save.plots}{Logical. Should Manhattan and QQ-plots be generated and
saved to the working directory for univariate GWAS? Default is TRUE.}

\item{thr.r2}{Value between 0 and 1. Threshold of r2 measure of linkage
disequilibrium. Markers in higher LD than this will be subset using clumping.}

\item{thr.m}{"sum" or "max". Type of threshold to use to clump values for
mash inputs. "sum" sums the -log10pvalues for each phenotype and uses
the maximum of this value as the threshold. "max" uses the maximum
-log10pvalue for each SNP across all of the univariate GWAS.}

\item{num.strong}{Integer. Number of SNPs used to derive data-driven covariance
matrix patterns, using markers with strong effects on phenotypes.}

\item{num.random}{Integer. Number of SNPs used to derive the correlation structure
of the null tests, and the mash fit on the null tests.}

\item{scale.phe}{Logical. Should effects for each phenotype be scaled to fall
between -1 and 1? Default is TRUE.}

\item{roll.size}{Integer. Used to create the svd for GWAS.}

\item{U.ed}{Mash data-driven covariance matrices. Specify these as a list or a path
to a file saved as an .rds. Creating these can be time-consuming, and
generating these once and reusing them for multiple mash runs can save time.}

\item{U.hyp}{Other covariance matrices for mash. Specify these as a list. These
matrices must have dimensions that match the number of phenotypes where
univariate GWAS ran successfully.}

\item{verbose}{Output some information on the iterations? Default is \code{TRUE}.}
}
\value{
A mash object made up of all phenotypes where univariate GWAS ran
successfully.
}
\description{
Though step-by-step GWAS, preparation of mash inputs, and mash
allows you the most flexibility and opportunities to check your results
for errors, once those sanity checks are complete, this function allows
you to go from a phenotype data.frame of a few phenotypes you want to
compare to a mash result. Some exception handling has been built into
this function, but the user should stay cautious and skeptical of any
results that seem 'too good to be true'.
}
