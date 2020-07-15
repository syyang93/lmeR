#' A wrapper for deciding either to run a default lm or a two-stage lm.  I always use the default.  Adapted from Vamsee (github.com/vkp3/pillalamarRi)
#' 
#' Runs linear regression to estimate effect of SCORE on gene expression/metabolite/epigenetic data
#' @description formula: Expr ~ SCORE + Cov + (1|Rcov), where Cov = fixed covariates and Rcov = random covariates
#' 
#' @param expr Expression vector (numeric) with length = #samples (numeric vector)
#' @param cov Regression covariates (fixed effects) in form [cov x samples] (data frame)
#' @param rcov Regression covariates (random effects) in the form [rcov x samples] (data frame)
#' @param SCORE mtDNA-CN, not included in `cov` (numeric vector)
#' @param omit.outlier (no longer supported) Whether or not you want to omit gene expression outliers 
#' @return A [1 x 6] vector output from an lm() like below:
#'         ['intercept', 'beta', 'SE', 't_value', 'pval', 'corr.rho']
#' 
#' Author: Stephanie Yang
#' Adapted from: Vamsee Pillalamarri
#' @export

# for testing:
# w.mt = readRDS('/Volumes/JHPCE/dcs01/active/projects/mesa/syang/R_objects/metab.w.mt.vis1.rds')
# expr <- w.mt[,which(colnames(w.mt) == "D.Gluconic.acid")]
# cov <- dplyr::select(w.mt, gender1, age1c, PC1, PC2, PC3, bmi1c)
# rcov <- dplyr::select(w.mt, race1c)
# SCORE <- w.mt[,which(colnames(w.mt) == "dpcrAdjMetric")]

run_lmer <- function(expr, cov, rcov, SCORE, omit.outlier = T, outlier_sd = 3) {
    
  expr <- as.numeric(expr)
  
  # expr <- scale(expr) # uncomment if you would like to scale expression
  # expr <- inv.norm.transform(expr) # uncomment if you would like to inverse normal transform expression
  
  # Create dataframe
  expr_cov <- cbind(SCORE, expr, cov, rcov)

  # Create equation
  lmer.form = paste0('expr~SCORE+',paste0(colnames(cov), collapse = '+'), '+', paste0('(1|',colnames(rcov), ')', collapse = '+'))
  
  # Run mixed model
  lmer.fit <- lmer(lmer.form, data = expr_cov, na.action = na.exclude)
  
  # Get summary
  lmer.fit.summary <- summary(lmer.fit)
  
  # Get corr
  cor_expr_score <-
    cor(expr, SCORE)
  
  # Capture beta, tval, standard error.
  lmer.res_ <-
    as.data.frame(t(coef(lmer.fit.summary)['SCORE',]))
  
  # get pval christina's way:
  pval <- pchisq(lmer.res_$`t value`^2,1,lower.tail=F)
  
  lmer.res_ <-
    cbind(lmer.res_, pval,
    #      t(confint(lmer.fit)['SCORE', ]),
          cor_expr_score)
  
  return(as.matrix(lmer.res_))
}
