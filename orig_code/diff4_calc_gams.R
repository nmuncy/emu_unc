library(mgcv)

switch_family <- function(dist_fam) {
  # Switch for determining desired family.
  #
  # Arguments:
  #   dist_fam (str) =  Desired family/distribution to use with data 
  #     (gamma, beta, or gaus)
  #
  # Returns:
  #   x_fam (str) = Family option for GAM/BAM function
  x_fam <- switch(dist_fam,
    "gamma" = "Gamma(link = \"logit\")",
    "beta" = "betar(link = \"logit\")",
    "gaus" = "gaussian()",
  )
  return(x_fam)
}


gam_G <- function(df, dist_fam) {
  # Conduct basic GAM analysis, G method.
  #
  # Model dti_fa for a tract, where sex has a
  # separate intercept, subject is random effect,
  # and using cubic spline for nodeID. K=50 determined
  # by testing.
  #
  # Arguments:
  #   df (dataframe) = Long-formatted dataframe for AFQ tract, containing
  #     columns of dti_fa, sex, subjectID, and nodeID
  #   dist_fam (str) = Desired family/distribution to use with data 
  #     (gamma, beta, or gaus)
  #
  # Returns:
  #   h_gam = GAM object
  
  h_family <- switch_family(dist_fam)
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 50),
  data = df,
  family = h_family,
  method = "fREML"
  )
  return(h_gam)
}


gam_Gcov <- function(df, dist_fam, cov) {
  # Conduct basic GAM analysis, G method with covariate.
  #
  # Similar to gam_G, but include covariate and control for
  # interaction of sex with covariate.
  #
  # Arguments:
  #   df (dataframe) = Long-formatted dataframe for AFQ tract, containing
  #     columns of dti_fa, sex, subjectID, and nodeID
  #   dist_fam (str) = Desired family/distribution to use with data 
  #     (gamma, beta, or gaus)
  #   cov (str) = Covariate column name
  #
  # Returns:
  #   h_gam = GAM object
  
  h_family <- switch_family(dist_fam)
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 50) +
    s(get(cov), by = sex),
  data = df,
  family = h_family,
  method = "fREML"
  )
  return(h_gam)
}


gam_GS <- function(df, dist_fam, col_group) {
  # Model data with global + group smooths, GS method.
  #
  # Similar to gam_G, but each group factor has own smooth.
  #
  # Arguments:
  #   df (dataframe) = Long-formatted dataframe for AFQ tract, containing
  #     columns of dti_fa, sex, subjectID, and nodeID
  #   dist_fam (str) = Desired family/distribution to use with data 
  #     (gamma, beta, or gaus)
  #   col_group (str) = Column name of factored grouping value
  #
  # Returns:
  #   h_gam = GAM object
  
  h_family <- switch_family(dist_fam)
  names(df)[names(df) == col_group] <- "h_group"
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 50, m = 2) +
    s(nodeID, h_group, bs = "fs", k = 50, m = 2),
  data = df,
  family = h_family,
  method = "fREML"
  )
  return(h_gam)
}


gam_GI <- function(df, dist_fam, col_group) {
  # Model data with global + group smooths, GI method.
  #
  # Similar to gam_G, but each group factor has own smooth and wiggliness.
  #
  # Arguments:
  #   df (dataframe) = Long-formatted dataframe for AFQ tract, containing
  #     columns of dti_fa, sex, subjectID, and nodeID
  #   dist_fam (str) = Desired family/distribution to use with data 
  #     (gamma, beta, or gaus)
  #   col_group (str) = Column name of factored grouping value
  #
  # Returns:
  #   h_gam = GAM object
  
  h_family <- switch_family(dist_fam)
  names(df)[names(df) == col_group] <- "h_group"
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(h_group, bs = "re") +
    s(nodeID, bs = "cr", k = 50, m = 2) +
    s(nodeID, by = h_group, bs = "cr", k = 50, m = 1),
  data = df,
  family = h_family,
  method = "fREML"
  )
  return(h_gam)
}


gam_GSOF <- function(df, dist_fam, col_group) {
  # Model data with global + group smooths (ordered factors).
  #
  # Similar to gam_GS_model, but use ordered factors so a difference
  # smooth can be calculated.
  #
  # Arguments:
  #   df (dataframe) = Long-formatted dataframe for AFQ tract, containing
  #     columns of dti_fa, sex, subjectID, and nodeID
  #   dist_fam (str) = Desired family/distribution to use with data 
  #     (gamma, beta, or gaus)
  #   col_group (str) = Column name of ordered factored grouping value
  #
  # Returns:
  #   h_gam = GAM object
  
  h_family <- switch_family(dist_fam)
  names(df)[names(df) == col_group] <- "h_group"
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 50, m = 2) +
    s(nodeID, by = h_group, bs = "cr", k = 50, m = 2),
  data = df,
  family = h_family,
  method = "fREML"
  )
  return(h_gam)
}


gam_Gintx <- function(df, dist_fam, col_group, cont_var) {
  # Plot the interaction of tract node, fa, group, covariate.
  #
  # Used for visualizing entire interaction of data.
  #
  # Arguments:
  #   df (dataframe) = Long-formatted dataframe for AFQ tract, containing
  #     columns of dti_fa, sex, subjectID, and nodeID
  #   dist_fam (str) = Desired family/distribution to use with data 
  #     (gamma, beta, or gaus)
  #   col_group (str) = Column name of factored grouping value
  #   cont_var (str) =  Continuous variable, which will interact with
  #     predicted
  #
  # Returns:
  #   h_gam = GAM object

  h_family <- switch_family(dist_fam)
  names(df)[names(df) == col_group] <- "h_group"
  names(df)[names(df) == cont_var] <- "h_var"
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, h_var, bs = c("cr", "tp"), k = c(50, 5), m = 2) +
    t2(
      nodeID, h_var, h_group,
      bs = c("cr", "tp", "re"), k = c(50, 5, 2), m = 2
    ),
  data = df,
  family = h_family,
  method = "fREML"
  )
  return(h_gam)
}


gam_GSintx <- function(df, dist_fam, col_group, cont_var) {
  # Model data with interaction by group.
  #
  # Model tract FA values, allow for interaction of nodeID, group, and
  # continuous variable (LGI, PPI, etc). Interaction term will not contain
  # main effect of nodeID i.e. no tract curvature.
  #
  # Arguments:
  #   df (dataframe) = Long-formatted dataframe for AFQ tract, containing
  #     columns of dti_fa, sex, subjectID, and nodeID
  #   dist_fam (str) = Desired family/distribution to use with data 
  #     (gamma, beta, or gaus)
  #   col_group (str) = Column name of factored grouping value
  #   cont_var (str) =  Continuous variable, which will interact with
  #     predicted
  #
  # Returns:
  #   h_gam = GAM object
  
  h_family <- switch_family(dist_fam)
  names(df)[names(df) == col_group] <- "h_group"
  names(df)[names(df) == cont_var] <- "h_var"
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 50, m = 2) +
    s(h_var, by = h_group, bs = "tp", k = 5, m = 2) +
    ti(
      nodeID, h_var,
      by = h_group, bs = c("cr", "tp"), k = c(50, 5), m = 2
    ),
  data = df,
  family = h_family,
  method = "fREML",
  discrete = T
  )
  return(h_gam)
}


gam_GSintxOF <- function(df, dist_fam, col_group, col_groupOF, cont_var) {
  # Model data with interaction by group (ordered factor).
  #
  # Model tract FA values, allow for interaction of nodeID, group, and
  # continuous variable (LGI, PPI, etc). Use ordered factor so interaction
  # of <cont_var> and nodeID difference from reference group can be
  # calculated.
  #
  # Arguments:
  #   df (dataframe) = Long-formatted dataframe for AFQ tract, containing
  #     columns of dti_fa, sex, subjectID, and nodeID
  #   dist_fam (str) = Desired family/distribution to use with data 
  #     (gamma, beta, or gaus)
  #   col_group (str) = Column name of factored grouping value
  #   col_groupOF (str) =  Column name of ordered grouping factor value,
  #     for make difference smooth for experimental group
  #   cont_var (str) =  Continuous variable, which will interact with
  #     predicted
  #
  # Returns:
  #   h_gam = GAM object
  
  h_family <- switch_family(dist_fam)
  names(df)[names(df) == col_group] <- "h_group"
  names(df)[names(df) == col_groupOF] <- "h_groupOF"
  names(df)[names(df) == cont_var] <- "h_var"
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 50, m = 2) +
    s(h_var, by = h_group, bs = "tp", k = 5, m = 2) +
    ti(nodeID, h_var, bs = c("cr", "tp"), k = c(50, 5), m = 2) +
    ti(
      nodeID, h_var,
      by = h_groupOF, bs = c("cr", "tp"), k = c(50, 5), m = 2
    ),
  data = df,
  family = h_family,
  method = "fREML",
  discrete = T
  )
  return(h_gam)
}
