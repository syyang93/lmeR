% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_lmer.R
\name{run_lmer}
\alias{run_lmer}
\title{A wrapper for deciding either to run a default lm or a two-stage lm.  I always use the default.  Adapted from Vamsee (github.com/vkp3/pillalamarRi)}
\usage{
run_lmer(expr, cov, rcov, SCORE, omit.outlier = T, outlier_sd = 3)
}
\arguments{
\item{expr}{Expression vector (numeric) with length = #samples (numeric vector)}

\item{cov}{Regression covariates (fixed effects) in form \link{cov x samples} (data frame)}

\item{rcov}{Regression covariates (random effects) in the form \link{rcov x samples} (data frame)}

\item{SCORE}{mtDNA-CN, not included in \code{cov} (numeric vector)}

\item{omit.outlier}{(no longer supported) Whether or not you want to omit gene expression outliers}
}
\value{
A \link{1 x 6} vector output from an lm() like below:
\link{'intercept', 'beta', 'SE', 't_value', 'pval', 'corr.rho'}

Author: Stephanie Yang
Adapted from: Vamsee Pillalamarri
}
\description{
formula: Expr ~ SCORE + Cov + (1|Rcov), where Cov = fixed covariates and Rcov = random covariates
}
\details{
Runs linear regression to estimate effect of SCORE on gene expression/metabolite/epigenetic data
}
