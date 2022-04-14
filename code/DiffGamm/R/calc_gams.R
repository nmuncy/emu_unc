#' Switch for determining desired family.
#'
#' @export
#' @param dist_fam Desired family/distribution to use with data (gamma, beta,
#' or gaus; str).
#' @return Family option for GAM/BAM function
switch_family <- function(dist_fam) {
  x_fam <- switch(dist_fam,
    "gamma" = "Gamma(link = \"logit\")",
    "beta" = "betar(link = \"logit\")",
    "gaus" = "gaussian()",
  )
  return(x_fam)
}

#' Conduct basic GAM analysis, G method.
#'
#' Model dti_fa for a tract, where sex has a
#' separate intercept, subject is random effect,
#' and using cubic spline for nodeID. K=50 determined
#' by testing.
#'
#' @export
#' @param df_tract Long-formatted dataframe for AFQ tract, containing
#' columns of dti_fa, sex, subjectID, and nodeID.
#' @param dist_fam Desired family/distribution to use with data (gamma, beta,
#' or gaus; str).
#' @return GAM object
#' @import mgcv
gam_model <- function(df_tract, dist_fam) {
  h_family <- switch_family(dist_fam)
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 50),
  data = df_tract,
  family = h_family,
  method = "fREML"
  )
  return(h_gam)
}

#' Conduct basic GAM analysis, G method with covariate.
#'
#' Similar to gam_model, but include covariate and control for
#' interaction of sex with covariate.
#'
#' @export
#' @param df_tract Long-formatted dataframe for AFQ tract, containing
#' columns of dti_fa, sex, subjectID, and nodeID.
#' @param dist_fam Desired family/distribution to use with data (gamma, beta,
#' or gaus; str).
#' @param cov Covariate column name (str).
#' @return GAM object
#' @import mgcv
gam_cov_model <- function(df_tract, dist_fam, cov) {
  h_family <- switch_family(dist_fam)
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 50) +
    s(get(cov), by = sex),
  data = df_tract,
  family = h_family,
  method = "fREML"
  )
  return(h_gam)
}

#' Model data with global + group smooths, GS method.
#'
#' Similar to gam_model, but each group factor has own smooth.
#'
#' @export
#' @param df_tract Long-formatted dataframe for AFQ tract, containing
#' columns of dti_fa, sex, subjectID, and nodeID.
#' @param dist_fam Desired family/distribution to use with data (gamma, beta,
#' or gaus; str).
#' @param col_group Column name of factored grouping value (str).
#' @return GAM object
#' @import mgcv
gam_GS_model <- function(df_tract, dist_fam, col_group) {
  h_family <- switch_family(dist_fam)
  names(df_tract)[names(df_tract) == col_group] <- "h_group"
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 50, m = 2) +
    s(nodeID, h_group, bs = "fs", k = 50, m = 2),
  data = df_tract,
  family = h_family,
  method = "fREML"
  )
  return(h_gam)
}

#' Model data with global + group smooths, GI method.
#'
#' Similar to gam_model, but each group factor has own smooth and wiggliness.
#'
#' @export
#' @param df_tract Long-formatted dataframe for AFQ tract, containing
#' columns of dti_fa, sex, subjectID, and nodeID.
#' @param dist_fam Desired family/distribution to use with data (gamma, beta,
#' or gaus; str).
#' @param col_group Column name of factored grouping value (str).
#' @return GAM object
#' @import mgcv
gam_GI_model <- function(df_tract, dist_fam, col_group) {
  h_family <- switch_family(dist_fam)
  names(df_tract)[names(df_tract) == col_group] <- "h_group"
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(h_group, bs = "re") +
    s(nodeID, bs = "cr", k = 50, m = 2) +
    s(nodeID, by = h_group, bs = "cr", k = 50, m = 1),
  data = df_tract,
  family = h_family,
  method = "fREML"
  )
  return(h_gam)
}

#' Model data with global + group smooths (ordered factors).
#'
#' Similar to gam_GS_model, but use ordered factors so a difference
#' smooth can be calculated.
#'
#' @export
#' @param df_tract Long-formatted dataframe for AFQ tract, containing
#' columns of dti_fa, sex, subjectID, and nodeID.
#' @param dist_fam Desired family/distribution to use with data (gamma, beta,
#' or gaus; str).
#' @param col_group Column name of ordered factored grouping value (str).
#' @return GAM object
#' @import mgcv
gam_GSOF_model <- function(df_tract, dist_fam, col_group) {
  h_family <- switch_family(dist_fam)
  names(df_tract)[names(df_tract) == col_group] <- "h_group"
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 50, m = 2) +
    s(nodeID, by = h_group, bs = "cr", k = 50, m = 2),
  data = df_tract,
  family = h_family,
  method = "fREML"
  )
  return(h_gam)
}

#' Model data with interaction by group.
#'
#' Model tract FA values, allow for interaction of nodeID, group, and
#' continuous variable (LGI, PPI, etc).
#'
#' @export
#' @param df_tract Long-formatted dataframe for AFQ tract, containing
#' columns of dti_fa, sex, subjectID, and nodeID.
#' @param dist_fam Desired family/distribution to use with data (gamma, beta,
#' or gaus; str).
#' @param col_group Column name of factored grouping value (str).
#' @param cont_var Continuous variable, which will interact with
#' predicted (str).
#' @return GAM object
#' @import mgcv
gam_intx_model <- function(df_tract, dist_fam, col_group, cont_var) {
  h_family <- switch_family(dist_fam)
  names(df_tract)[names(df_tract) == col_group] <- "h_group"
  names(df_tract)[names(df_tract) == cont_var] <- "h_var"
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 50, m = 1) +
    s(h_var, by = h_group, bs = "tp", k = 10, m = 2) +
    ti(
      nodeID, h_var, by = h_group, bs = c("cr", "tp"), k = c(50, 10), m = 2
    ),
  data = df_tract,
  family = h_family,
  method = "fREML",
  discrete = T
  )
  return(h_gam)
}

#' Model data with interaction by group (ordered factor).
#'
#' Model tract FA values, allow for interaction of nodeID, group, and
#' continuous variable (LGI, PPI, etc). Use ordered factor so interaction
#' of <cont_var> and nodeID difference from reference group can be
#' calculated.
#'
#' @export
#' @param df_tract Long-formatted dataframe for AFQ tract, containing
#' columns of dti_fa, sex, subjectID, and nodeID.
#' @param dist_fam Desired family/distribution to use with data (gamma, beta,
#' or gaus; str).
#' @param col_group Column name of ordered factored grouping value (str).
#' @param cont_var Continuous variable, which will interact with
#' predicted (str).
#' @return GAM object
#' @import mgcv
gam_intxOF_model <- function(df_tract, dist_fam, col_group, cont_var) {
  h_family <- switch_family(dist_fam)
  names(df_tract)[names(df_tract) == col_group] <- "h_group"
  names(df_tract)[names(df_tract) == cont_var] <- "h_var"
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 50, m = 1) +
    s(h_var, by = h_group, bs = "tp", k = 10, m = 2) +
    ti(nodeID, h_var, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    ti(
      nodeID, h_var, by = h_group, bs = c("cr", "tp"), k = c(50, 10), m = 2
    ),
  data = df_tract,
  family = h_family,
  method = "fREML",
  discrete = T
  )
  return(h_gam)
}
