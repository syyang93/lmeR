#' Function that will run a regression for your variable of interest across all genes/metabolites/CpGs, in parallel.  Adapted from Vamsee (github.com/vkp3/pillalamarRi)
#' 
#' @param tx_expr Expression matrix in form: [genes x samples]. Will be converted to a list(!) of gene-expr vectors (if not input as list)
#' @param gene.ids Character vector of gene IDs, corresponding to rows in `tx_expr` [genes x samples]
#' @param cov Regression covariates (fixed effects) [cov x samples]
#' @param rcov Regression covariates to be treated as random effects
#' @param SCORE Main covariate to be permuted (not included in `cov`)
#' @param omit.outlier Whether or not you want to omit gene expression outliers
#' @param num.cores The number of cores you would like to use
#' @export 
#' 
#' @return All regression coefficients for lmer(gene expression ~ SCORE + cov + (1|rcov)) for all genes
#' 
#' @examples
#' tx_expr <- w.mt[,which(colnames(w.mt) %in% make.names(metablist))]
#' cov <- w.mt[,which(colnames(w.mt) %in% c('sex', 'age1c', 'PC1', 'PC2', 'PC3', 'bmi', 'Exam'))]
#' rcov <- w.mt[,which(colnames(w.mt) %in% c('idno'))]
#' SCORE <- w.mt[,which(colnames(w.mt) == "dpcrAdjMetric")]
#' lm_res.sort <- run.all.lmers(tx_expr, cov, rcov, colnames(tx_expr), SCORE, omit.outlier = T, num.cores = 10)
#' 
#' @note 
# Can't time test on local machine, since you can't run over multiple cores on your poor laptop.
# Benchmarked on data with 252 tested metabolites, 8 covariates, 1 random covariate, 1516 total samples on JHPCE with 2 GB. 
# For loop took 17.83 seconds, pbapply took 2.15 seconds

run.all.lmers <- function(tx_expr, cov, rcov, gene.ids, SCORE, omit.outlier = T, num.cores = 10)
{
  require(pbapply)
  require(lme4)
  start = Sys.time()
  lmer.res <-
    pblapply(tx_expr,            # Expression vector list for `pbapply::pblapply`
             run_lmer,             # This function
             cov = cov,           # Covariate matrix (fixed effects), as desribed above
             rcov = rcov,        # Covariates that are random varaibles
             SCORE = SCORE,       # PRS
             omit.outlier = omit.outlier,
             cl = num.cores)      # Number of cores to parallelize over
  end = Sys.time()
  total.time.confint = end-start
  total.time.confint
  # total.time.noconfint = end-start
  
  lmer.res <- simplify2array(lmer.res, higher=F)
  
  rownames(lmer.res) <-
    c(# 'intercept',
      'beta',
      'SE',
      't_value',
      'pval',
      # 'conf.low',
      # 'conf.high',
      'corr.rho')
  colnames(lmer.res) <- gene.ids
  lmer.res <- as.data.frame(t(lmer.res))
  return(lmer.res)
}
