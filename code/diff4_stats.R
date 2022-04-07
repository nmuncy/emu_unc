library("fitdistrplus")
library("itsadug")
library("tidymv")
library("dplyr")
library("mgcViz")
library("tools")
library("tidyr")
library("devtools")
install_local(path = ".")
library("DiffGamm")


# Functions ----
tract_fam <- function(tract){
  # Set family for each tract.
  #
  # These were determined through comparing various models
  # via itsadug::compareML.
  #
  # Arguments:
  #   tract (str) = AFQ tract name
  #
  # Returns:
  #   family (str) for use in DiffGamm::switch_family
  h_fam <- switch(tract,
                  "UNC_L" = "gamma",
                  "UNC_R" = "gaus",
                  "CGC_L" = "gaus",
                  "CGC_R" = "gaus",
  )
  return(h_fam)
}

write_gam_stats <- function(gam_obj, out_dir, gam_type, tract){
  # Write summary stats of GAM object.
  #
  # Arguments:
  #   gam_obj (obj) = GAM object returned by mgcv
  #   out_dir (str) = path to output dir
  #   gam_type (str) = diffentiate type of model used
  #   tract (str) = AFQ tract name
  capture.output(
    summary(gam_obj),
    file = paste0(
      out_dir, "/Stats_GAM_", tract, "_", gam_type, ".txt"
    )
  )
}

write_compare_stats <- function(model_a, model_b, tract, out_dir, out_str){
  # Write model comparison stats of two GAMs
  #
  # Arguments:
  #   model_a (obj) = GAM object returned by mgcv
  #   model_b (obj) = another GAM object returned by mgcv
  #   tract (str) = AFQ tract name
  #   out_dir (str) = path to output dir
  #   out_str (str) = extra string for specifying name of out file
  capture.output(
    compareML(model_a, model_b),
    file = paste0(
      out_dir, "/Stats_GAM_", tract, "_compare_", out_str, ".txt"
    )
  )
}


# Set Up ----

# set paths
# proj_dir <- file_path_as_absolute(paste0(getwd(), "/.."))
proj_dir <- "/Users/nmuncy/Projects/emu_unc"
data_dir <- paste0(proj_dir, "/data")
out_dir <- paste0(proj_dir, "/stats")

# capture session
capture.output(sessionInfo(), file = paste0(data_dir, "/R_session_info.txt"))

# import data, setup factors
df_afq <- read.csv(paste0(data_dir, "/AFQ_dataframe.csv"))
df_afq$sex <- factor(df_afq$sex)
ind_exp <- which(df_afq$dx_group == "Pat")
df_afq[ind_exp, ]$dx_group <- "Exp"
df_afq$dx_group <- factor(df_afq$dx_group)
df_afq$subjectID <- factor(df_afq$subjectID)
df_afq$dx_groupOF <- factor(df_afq$dx_group, ordered = T)

# clip tails
ind_keep <- which(
  df_afq$nodeID >= 10 & df_afq$nodeID <= 89
)
df_afq <- df_afq[ind_keep, ]
rm(ind_keep)
rm(ind_exp)

# incorporate PPI values of ses-S1 task-study
seed_list <- c("NSlacc", "NSldmpfc", "NSlsfs")
beh_list <- c("SPnegLF", "SPneuLF")
subj_list <- as.character(unique(df_afq$subjectID))
for (seed in seed_list) {
  df_ppi <- read.csv(
    paste0(data_dir, "/df_ses-S1_task-study_amgL-", seed, ".csv")
  )
  for (beh in beh_list) {
    h_col <- paste(seed, beh, sep = "_")
    df_afq[, h_col] <- NA
    for (subj in subj_list) {
      ind_afq <- which(df_afq$subjectID == subj)
      ind_ppi <- which(df_ppi$subj == paste0("sub-", subj))
      if (length(ind_ppi) == 0) {
        next
      }
      df_afq[ind_afq, h_col] <- df_ppi[ind_ppi, beh]
    }
  }
}
rm(df_ppi)



# Model Specification ----
#
# Determine the model that best fits the various tracts.
# Modeling individual tracts determined for that:
#   a) k=50 was a sufficient basis dimension for all tracts,
#   b) the required family arguments,
#   c) PDS did not increase model fit for any tract,
#   d) GS fit better than G, GI did not increase fit for all tracts.
#
# As each GAM for the tracts is very similar, only differing
# in the distribution, we can loop through the tracts.
#
# The dxGS models are saved.

tract_list <- c("UNC_L", "UNC_R", "CGC_L", "CGC_R")
for(tract in tract_list){

  # subset df_afq, keep people w/dx for group modeling
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract <- df_tract %>% drop_na(dx)

  # conduct basic model
  tract_dist <- tract_fam(tract)
  tract_G <- gam_model(df_tract, tract_dist)
  write_gam_stats(tract_G, out_dir, "G", tract)

  # does controlling for PDS help model fit
  tract_G_pds <- gam_cov_model(df_tract, tract_dist, "pds")
  write_gam_stats(tract_G_pds, out_dir, "G-PDS", tract)
  write_compare_stats(tract_G, tract_G_pds, tract, out_dir, "G-PDS")

  # test if group smooths increase fit, save model
  gam_file <- paste0(out_dir, "/Model_", tract, "_dxGS.Rda")
  if(!file.exists(gam_file)){
    h_gam <- gam_GS_model(df_tract, tract_dist, "dx_group")
    saveRDS(h_gam, file = gam_file)
    rm(h_gam)
  }
  tract_GS <- readRDS(gam_file)
  write_gam_stats(tract_GS, out_dir, "GS", tract)
  write_compare_stats(tract_GS, tract_G, tract, out_dir, "GS-G")

  # test if group wiggliness increases fit
  tract_GI <- gam_GI_model(df_tract, tract_dist, "dx_group")
  write_gam_stats(tract_GS, out_dir, "GI", tract)
  write_compare_stats(tract_GS, tract_GI, tract, out_dir, "GS-GI")

  # test if group smooths differ
  tract_GSOF <- gam_GSOF_model(df_tract, tract_dist, "dx_groupOF")
  write_gam_stats(tract_GSOF, out_dir, "GSOF", tract)

  # draw plots
  plot_tract_GS <- getViz(tract_GS)
  draw_global_smooth(plot_tract_GS, 2, tract, out_dir)
  draw_group_smooth(plot_tract_GS, 3, tract, out_dir)
  plot_tract_GSOF <- getViz(tract_GSOF)
  draw_group_smooth_diff(plot_tract_GSOF, 3, tract, out_dir)

  # clean env
  rm(tract_G)
  rm(tract_G_pds)
  rm(tract_GS)
  rm(tract_GSOF)
  rm(tract_GI)
  rm(plot_tract_GS)
  rm(plot_tract_GSOF)
}





# L. Unc LGI Interaction ----
#
# 1) Investigate tract-group-negLGI intx
# 2) Investigate tract-group-neuLGI intx
#
# Save interaction models as they take a while to generate

# 1) L. Unc group interaction with negLGI
# generate model if necessary
gam_file <- paste0(out_dir, "/Data_lunc_dxGS_neg.Rda")
if (!file.exists(gam_file)) {
  # h_gam <- bam(dti_fa ~ sex +
  #   s(subjectID, bs = "re") +
  #   te(nodeID, lgi_neg, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
  #   t2(
  #     nodeID, lgi_neg, dx_group,
  #     bs = c("cr", "tp", "re"),
  #     k = c(50, 10, 2),
  #     m = 2,
  #     full = TRUE
  #   ),
  # data = df_tract,
  # family = Gamma(link = "logit"),
  # method = "fREML"
  # )
  h_gam <- gam_intx_model(df_tract, "gamma", "dx_group", "lgi_neg")
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}

# read in model, get stats
lunc_dxGS_neg <- readRDS(gam_file)
compareML(lunc_dxGS, lunc_dxGS_neg) # lunc_dxGS_neg preferred
summary(lunc_dxGS_neg)

# draw, unpack tract-LGI intx by group
plot_lunc_dxGS_neg <- getViz(lunc_dxGS_neg)
draw_smooth_intx(plot_lunc_dxGS_neg, 2, "UNC_L", "lgi_neg", out_dir)
draw_group_intx(df_tract, lunc_dxGS_neg, "UNC_L", "lgi_neg", out_dir)

# test if experiment group differs from control (reference group)
gam_file <- paste0(out_dir, "/Data_lunc_dxGS_negOF.Rda")
if (!file.exists(gam_file)) {
  # h_gam <- bam(dti_fa ~ sex +
  #   s(subjectID, bs = "re") +
  #   te(nodeID, lgi_neg, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
  #   t2(
  #     nodeID, lgi_neg,
  #     by = dx_groupOF,
  #     bs = c("cr", "tp"),
  #     k = c(50, 10),
  #     m = 2,
  #     full = TRUE
  #   ),
  # data = df_tract,
  # family = Gamma(link = "logit"),
  # method = "fREML"
  # )
  h_gam <- gam_intxOF_model(df_tract, "gamma", "dx_groupOF", "lgi_neg")
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}

# read in model, get stats
lunc_dxGS_negOF <- readRDS(gam_file)
summary(lunc_dxGS_negOF)

# draw reference, difference interaction smooths
plot_lunc_dxGS_negOF <- getViz(lunc_dxGS_negOF)
draw_group_intx_ref(plot_lunc_dxGS_negOF, 2, "UNC_L", "lgi_neg", out_dir)
draw_group_intx_diff(plot_lunc_dxGS_negOF, 3, "UNC_L", "lgi_neg", out_dir)

# clean up
rm(lunc_dxGS_neg)
rm(plot_lunc_dxGS_neg)
rm(lunc_dxGS_negOF)
rm(plot_lunc_dxGS_negOF)


# 2) L. Unc group interaction with negLGI
#
# Same as above, but with neutral stimuli.
gam_file <- paste0(out_dir, "/Data_lunc_dxGS_neu.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, lgi_neu, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, lgi_neu, dx_group,
      bs = c("cr", "tp", "re"),
      k = c(50, 10, 2),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = Gamma(link = "logit"),
  method = "fREML"
  )
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}
lunc_dxGS_neu <- readRDS(gam_file)
compareML(lunc_dxGS, lunc_dxGS_neu) # lunc_dxGS_neu preferred
summary(lunc_dxGS_neu)

# draw, unpack tract-LGI intx by group
plot_lunc_dxGS_neu <- getViz(lunc_dxGS_neu)
draw_smooth_intx(plot_lunc_dxGS_neu, 2, "UNC_L", "lgi_neu", out_dir)
draw_group_intx(df_tract, lunc_dxGS_neu, "UNC_L", "lgi_neu", out_dir)

# test if experiment group differs from control (reference group)
gam_file <- paste0(out_dir, "/Data_lunc_dxGS_neuOF.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, lgi_neu, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, lgi_neu,
      by = dx_groupOF,
      bs = c("cr", "tp"),
      k = c(50, 10),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = Gamma(link = "logit"),
  method = "fREML"
  )
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}

# read in model, get stats, draw plots
lunc_dxGS_neuOF <- readRDS(gam_file)
summary(lunc_dxGS_neuOF)
plot_lunc_dxGS_neuOF <- getViz(lunc_dxGS_neuOF)
draw_group_intx_ref(plot_lunc_dxGS_neuOF, 2, "UNC_L", "lgi_neu", out_dir)
draw_group_intx_diff(plot_lunc_dxGS_neuOF, 3, "UNC_L", "lgi_neu", out_dir)

# clean up
rm(lunc_dxGS_neu)
rm(plot_lunc_dxGS_neu)
rm(lunc_dxGS_neuOF)
rm(plot_lunc_dxGS_neuOF)
rm(lunc_dxGS)

# # double-check OF method with plot_diff2 (which uses GI method)
# lunc_dxGS_neu2 <- bam(dti_fa ~ sex +
#   s(subjectID, bs = "re") +
#   te(nodeID, lgi_neu, bs = c("cr", "tp"), k = c(50, 10), by = dx_group),
# data = df_tract,
# family = gaussian(),
# method = "fREML"
# )
#
# plot(lunc_dxGS_neu2)
# plot_diff2(
#   lunc_dxGS_neu2,
#   view = c("nodeID", "lgi_neu"),
#   comp = list(dx_group = c("Pat", "Con")),
#   zlim = NULL,
#   se = 1.96,
#   color = "topo",
#   show.diff = T
# )


# L. Unc NSlacc PPI Interaction ----
#
# 1) Investigate tract-group-SPnegLFA intx for NSlacc
# 2) Investigate tract-group-SPneuLFA intx for NSlacc

# take complete behavior cases
df_tract <- df_tract %>% drop_na(NSlacc_SPnegLF)
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete = F)

adjust_outliers <- function(df_tract) {
  # Replace outliers with max/min values.
  #
  # Arguments:
  #
  # Returns:
  #

  col_list <- c(
    "NSlacc_SPnegLF",
    "NSlacc_SPneuLF",
    "NSldmpfc_SPnegLF",
    "NSldmpfc_SPneuLF",
    "NSlsfs_SPnegLF",
    "NSlsfs_SPneuLF"
  )

  df <- df_tract[which(df_tract$nodeID == 10), ]
  subj_list <- as.character(df$subjectID)

  for (col_name in col_list) {
    # find min/max
    h_iqr <- IQR(df[, col_name], na.rm = TRUE)
    h_quant <- quantile(df[, col_name], na.rm = TRUE, names = FALSE)
    h_min <- h_quant[2] - (1.5 * h_iqr)
    h_max <- h_quant[4] + (1.5 * h_iqr)

    # detect, replace outliers
    ind_out <- which(df[, col_name] < h_min | df[, col_name] > h_max)
    for (ind in ind_out) {
      if (df[ind, col_name] > h_max) {
        df[ind, col_name] <- h_max
      } else if (df[ind, col_name] < h_min) {
        df[ind, col_name] <- h_min
      }
    }

    # fill df_tract with orig/updated values
    for (subj in subj_list) {
      ind_df <- which(df$subjectID == subj)
      ind_tract <- which(df_tract$subjectID == subj)
      df_tract[ind_tract, col_name] <- df[ind_df, col_name]
    }
  }
  return(df_tract)
}
df_tract <- adjust_outliers(df_tract)


# 1) L. Unc group interaction with SPnegLFA
# generate model if necessary
gam_file <- paste0(out_dir, "/Data_lunc_dxGS_NSlacc-PPIneg.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, NSlacc_SPnegLF, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, NSlacc_SPnegLF, dx_group,
      bs = c("cr", "tp", "re"),
      k = c(50, 10, 2),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = Gamma(link = "logit"),
  method = "fREML"
  )
  # gam.check(h_gam, rep = 1000)
  # summary(h_gam)
  # plot(h_gam)

  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}

# read in model, get stats
lunc_dxGS_NSlacc_PPIneg <- readRDS(gam_file)
summary(lunc_dxGS_NSlacc_PPIneg)

# draw, unpack tract-LGI intx by group
plot_lunc_dxGS_NSlacc_PPIneg <- getViz(lunc_dxGS_NSlacc_PPIneg)
plot(sm(plot_lunc_dxGS_NSlacc_PPIneg, 1))
plot(sm(plot_lunc_dxGS_NSlacc_PPIneg, 2))

draw_smooth_intx(plot_lunc_dxGS_NSlacc_PPIneg, 2, "UNC_L", "NSlacc_SPnegLF", out_dir)
draw_group_intx(df_tract, lunc_dxGS_NSlacc_PPIneg, "UNC_L", "NSlacc_SPnegLF", out_dir)

# test if experiment group differs from control (reference group)
gam_file <- paste0(out_dir, "/Data_lunc_dxGS_NSlacc-PPIneg_OF.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, NSlacc_SPnegLF, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, NSlacc_SPnegLF,
      by = dx_groupOF,
      bs = c("cr", "tp"),
      k = c(50, 10),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = Gamma(link = "logit"),
  method = "fREML"
  )
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}

# read in model, get stats
lunc_dxGS_NSlacc_PPIneg_OF <- readRDS(gam_file)
summary(lunc_dxGS_NSlacc_PPIneg_OF)

# draw reference, difference interaction smooths
plot_lunc_dxGS_NSlacc_PPIneg_OF <- getViz(lunc_dxGS_NSlacc_PPIneg_OF)
draw_group_intx_ref(plot_lunc_dxGS_NSlacc_PPIneg_OF, 2, "UNC_L", "NSlacc_SPnegLF", out_dir)
draw_group_intx_diff(plot_lunc_dxGS_NSlacc_PPIneg_OF, 3, "UNC_L", "NSlacc_SPnegLF", out_dir)

# # clean up
# rm(lunc_dxGS_neg)
# rm(plot_lunc_dxGS_neg)
# rm(lunc_dxGS_negOF)
# rm(plot_lunc_dxGS_negOF)


# R. Unc LGI Interaction ----
#
# Same as L. UNC, but with runc.

# 1) R. Unc GS neg LGI dx intx
gam_file <- paste0(out_dir, "/Data_runc_dxGS_neg.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, lgi_neg, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, lgi_neg, dx_group,
      bs = c("cr", "tp", "re"),
      k = c(50, 10, 2),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
  )
  gam.check(h_gam, rep = 1000)
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}
runc_dxGS_neg <- readRDS(gam_file)
compareML(runc_dxGS, runc_dxGS_neg) # runc_dxGS_neg preferred
summary(runc_dxGS_neg)

# draw, unpack tract-LGI intx by group
plot_runc_dxGS_neg <- getViz(runc_dxGS_neg)
draw_smooth_intx(plot_runc_dxGS_neg, 2, "UNC_R", "lgi_neg", out_dir)
draw_group_intx(df_tract, runc_dxGS_neg, "UNC_R", "lgi_neg", out_dir)

# test if experiment group differs from control (reference group)
gam_file <- paste0(out_dir, "/Data_runc_dxGS_negOF.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, lgi_neg, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, lgi_neg,
      by = dx_groupOF,
      bs = c("cr", "tp"),
      k = c(50, 10),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
  )
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}

# read in model, get stats
runc_dxGS_negOF <- readRDS(gam_file)
summary(runc_dxGS_negOF)

# draw reference, difference interaction smooths
plot_runc_dxGS_negOF <- getViz(runc_dxGS_negOF)
draw_group_intx_ref(plot_runc_dxGS_negOF, 2, "UNC_R", "lgi_neg", out_dir)
draw_group_intx_diff(plot_runc_dxGS_negOF, 3, "UNC_R", "lgi_neg", out_dir)

# clean up
rm(runc_dxGS_neg)
rm(plot_runc_dxGS_neg)
rm(runc_dxGS_negOF)
rm(plot_runc_dxGS_negOF)


# 2) R. Unc GS neu LGI dx intx
gam_file <- paste0(out_dir, "/Data_runc_dxGS_neu.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, lgi_neu, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, lgi_neu, dx_group,
      bs = c("cr", "tp", "re"),
      k = c(50, 10, 2),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
  )
  gam.check(h_gam, rep = 1000)
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}
runc_dxGS_neu <- readRDS(gam_file)
compareML(runc_dxGS, runc_dxGS_neu) # runc_dxGS_neu preferred
summary(runc_dxGS_neu)

# draw, unpack tract-LGI intx by group
plot_runc_dxGS_neu <- getViz(runc_dxGS_neu)
draw_smooth_intx(plot_runc_dxGS_neu, 2, "UNC_R", "lgi_neu", out_dir)
draw_group_intx(df_tract, runc_dxGS_neu, "UNC_R", "lgi_neu", out_dir)

# test if experiment group differs from control (reference group)
gam_file <- paste0(out_dir, "/Data_runc_dxGS_neuOF.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, lgi_neu, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, lgi_neu,
      by = dx_groupOF,
      bs = c("cr", "tp"),
      k = c(50, 10),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
  )
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}

# read in model, get stats
runc_dxGS_neuOF <- readRDS(gam_file)
summary(runc_dxGS_neuOF)

# draw reference, difference interaction smooths
plot_runc_dxGS_neuOF <- getViz(runc_dxGS_neuOF)
draw_group_intx_ref(plot_runc_dxGS_neuOF, 2, "UNC_R", "lgi_neu", out_dir)
draw_group_intx_diff(plot_runc_dxGS_neuOF, 3, "UNC_R", "lgi_neu", out_dir)

# clean up
rm(runc_dxGS_neu)
rm(plot_runc_dxGS_neu)
rm(runc_dxGS_neuOF)
rm(plot_runc_dxGS_neuOF)
rm(runc_dxGS)


# L. Cing LGI Interaction ----
#
# Same as L. UNC, but with lcgc.

# 1) L. Cing GS neg LGI dx intx
gam_file <- paste0(out_dir, "/Data_lcgc_dxGS_neg.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, lgi_neg, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, lgi_neg, dx_group,
      bs = c("cr", "tp", "re"),
      k = c(50, 10, 2),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
  )
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}
lcgc_dxGS_neg <- readRDS(gam_file)
compareML(lcgc_dxGS, lcgc_dxGS_neg) # lcgc_dxGS_neg preferred
summary(lcgc_dxGS_neg)

# draw, unpack tract-LGI intx by group
plot_lcgc_dxGS_neg <- getViz(lcgc_dxGS_neg)
draw_smooth_intx(plot_lcgc_dxGS_neg, 2, "CGC_L", "lgi_neg", out_dir)
draw_group_intx(df_tract, lcgc_dxGS_neg, "CGC_L", "lgi_neg", out_dir)

# test if experiment group differs from control (reference group)
gam_file <- paste0(out_dir, "/Data_lcgc_dxGS_negOF.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, lgi_neg, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, lgi_neg,
      by = dx_groupOF,
      bs = c("cr", "tp"),
      k = c(50, 10),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
  )
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}

# read in model, get stats
lcgc_dxGS_negOF <- readRDS(gam_file)
summary(lcgc_dxGS_negOF)

# draw reference, difference interaction smooths
plot_lcgc_dxGS_negOF <- getViz(lcgc_dxGS_negOF)
draw_group_intx_ref(plot_lcgc_dxGS_negOF, 2, "CGC_L", "lgi_neg", out_dir)
draw_group_intx_diff(plot_lcgc_dxGS_negOF, 3, "CGC_L", "lgi_neg", out_dir)

# clean up
rm(lcgc_dxGS_neg)
rm(plot_lcgc_dxGS_neg)
rm(lcgc_dxGS_negOF)
rm(plot_lcgc_dxGS_negOF)

# 2) L. Cing GS neu LGI dx intx
gam_file <- paste0(out_dir, "/Data_lcgc_dxGS_neu.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, lgi_neu, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, lgi_neu, dx_group,
      bs = c("cr", "tp", "re"),
      k = c(50, 10, 2),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
  )
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}
lcgc_dxGS_neu <- readRDS(gam_file)
compareML(lcgc_dxGS, lcgc_dxGS_neu) # lcgc_dxGS_neu preferred
summary(lcgc_dxGS_neu)

# draw, unpack tract-LGI intx by group
plot_lcgc_dxGS_neu <- getViz(lcgc_dxGS_neu)
draw_smooth_intx(plot_lcgc_dxGS_neu, 2, "CGC_L", "lgi_neu", out_dir)
draw_group_intx(df_tract, lcgc_dxGS_neu, "CGC_L", "lgi_neu", out_dir)

# test if experiment group differs from control (reference group)
gam_file <- paste0(out_dir, "/Data_lcgc_dxGS_neuOF.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, lgi_neu, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, lgi_neu,
      by = dx_groupOF,
      bs = c("cr", "tp"),
      k = c(50, 10),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
  )
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}

# read in model, get stats
lcgc_dxGS_neuOF <- readRDS(gam_file)
summary(lcgc_dxGS_neuOF)

# draw reference, difference interaction smooths
plot_lcgc_dxGS_neuOF <- getViz(lcgc_dxGS_neuOF)
draw_group_intx_ref(plot_lcgc_dxGS_neuOF, 2, "CGC_L", "lgi_neu", out_dir)
draw_group_intx_diff(plot_lcgc_dxGS_neuOF, 3, "CGC_L", "lgi_neu", out_dir)

# clean up
rm(lcgc_dxGS_neu)
rm(plot_lcgc_dxGS_neu)
rm(lcgc_dxGS_neuOF)
rm(plot_lcgc_dxGS_neuOF)
rm(lcgc_dxGS)


# R. Cing Model Specification ----
#
# Same as L. Unc Model Specification, but for cgc_r.

# subset df_afq, take complete cases
df_tract <- df_afq[which(df_afq$tractID == "CGC_R"), ]
df_tract <- df_tract %>% drop_na(dx)

# 1) determine distribution
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete = F)

# build gam
rcgc_gaus <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(rcgc_gaus, rep = 1000)

# 2) pds effect
rcgc_pds <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50) +
  s(pds, by = sex),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(rcgc_pds, rep = 1000)
summary(rcgc_pds) # sex effect with females
plot(rcgc_pds)
compareML(rcgc_gaus, rcgc_pds) # gaus fits better
rm(rcgc_pds)

# 3) GS model by diagnosis
rcgc_dxGS <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(rcgc_dxGS, rep = 1000)
compareML(rcgc_gaus, rcgc_dxGS) # no real diff, cont w/dxGS method
summary(rcgc_dxGS)
rm(rcgc_gaus)

# GI model by diagnosis
rcgc_dxGI <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(dx_group, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_group, bs = "cr", k = 50, m = 1),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(rcgc_dxGI, rep = 1000)
compareML(rcgc_dxGI, rcgc_dxGS) # no real difference - continuing with dxGS
rm(rcgc_dxGI)

# draw
plot_rcgc_dxGS <- getViz(rcgc_dxGS)
draw_global_smooth(plot_rcgc_dxGS, 2, "CGC_R", out_dir)
draw_group_smooth(plot_rcgc_dxGS, 3, "CGC_R", out_dir)
rm(plot_rcgc_dxGS)


# 4) identify nodes that differ (make difference smooth)
rcgc_dxGS_OF <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_groupOF, bs = "cr", k = 50, m = 2),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(rcgc_dxGS_OF, rep = 1000)
summary(rcgc_dxGS_OF) # no real group diff
plot_rcgc_dxGS_OF <- getViz(rcgc_dxGS_OF)
draw_group_smooth_diff(plot_rcgc_dxGS_OF, 3, "CGC_R", out_dir)
rm(rcgc_dxGS_OF)
rm(plot_rcgc_dxGS_OF)


# R. Cing LGI Interaction ----
#
# Same as L. UNC, but with rcgc.

# 1) R. Cing GS neg LGI dx intx
gam_file <- paste0(out_dir, "/Data_rcgc_dxGS_neg.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, lgi_neg, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, lgi_neg, dx_group,
      bs = c("cr", "tp", "re"),
      k = c(50, 10, 2),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
  )
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}
rcgc_dxGS_neg <- readRDS(gam_file)
compareML(rcgc_dxGS, rcgc_dxGS_neg) # rcgc_dxGS_neg preferred
summary(rcgc_dxGS_neg)

# draw, unpack tract-LGI intx by group
plot_rcgc_dxGS_neg <- getViz(rcgc_dxGS_neg)
draw_smooth_intx(plot_rcgc_dxGS_neg, 2, "CGC_R", "lgi_neg", out_dir)
draw_group_intx(df_tract, rcgc_dxGS_neg, "CGC_R", "lgi_neg", out_dir)

# test if experiment group differs from control (reference group)
gam_file <- paste0(out_dir, "/Data_rcgc_dxGS_negOF.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, lgi_neg, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, lgi_neg,
      by = dx_groupOF,
      bs = c("cr", "tp"),
      k = c(50, 10),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
  )
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}

# read in model, get stats
rcgc_dxGS_negOF <- readRDS(gam_file)
summary(rcgc_dxGS_negOF)

# draw reference, difference interaction smooths
plot_rcgc_dxGS_negOF <- getViz(rcgc_dxGS_negOF)
draw_group_intx_ref(plot_rcgc_dxGS_negOF, 2, "CGC_R", "lgi_neg", out_dir)
draw_group_intx_diff(plot_rcgc_dxGS_negOF, 3, "CGC_R", "lgi_neg", out_dir)

# clean up
rm(rcgc_dxGS_neg)
rm(plot_rcgc_dxGS_neg)
rm(rcgc_dxGS_negOF)
rm(plot_rcgc_dxGS_negOF)

# 2) R. Cing GS neu LGI dx intx
gam_file <- paste0(out_dir, "/Data_rcgc_dxGS_neu.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, lgi_neu, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, lgi_neu, dx_group,
      bs = c("cr", "tp", "re"),
      k = c(50, 10, 2),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
  )
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}
rcgc_dxGS_neu <- readRDS(gam_file)
compareML(rcgc_dxGS, rcgc_dxGS_neu) # rcgc_dxGS_neu preferred
summary(rcgc_dxGS_neu)

# draw, unpack tract-LGI intx by group
plot_rcgc_dxGS_neu <- getViz(rcgc_dxGS_neu)
draw_smooth_intx(plot_rcgc_dxGS_neu, 2, "CGC_R", "lgi_neu", out_dir)
draw_group_intx(df_tract, rcgc_dxGS_neu, "CGC_R", "lgi_neu", out_dir)

# test if experiment group differs from control (reference group)
gam_file <- paste0(out_dir, "/Data_rcgc_dxGS_neuOF.Rda")
if (!file.exists(gam_file)) {
  h_gam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, lgi_neu, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    t2(
      nodeID, lgi_neu,
      by = dx_groupOF,
      bs = c("cr", "tp"),
      k = c(50, 10),
      m = 2,
      full = TRUE
    ),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
  )
  saveRDS(h_gam, file = gam_file)
  rm(h_gam)
}

# read in model, get stats
rcgc_dxGS_neuOF <- readRDS(gam_file)
summary(rcgc_dxGS_neuOF)

# draw reference, difference interaction smooths
plot_rcgc_dxGS_neuOF <- getViz(rcgc_dxGS_neuOF)
draw_group_intx_ref(plot_rcgc_dxGS_neuOF, 2, "CGC_R", "lgi_neu", out_dir)
draw_group_intx_diff(plot_rcgc_dxGS_neuOF, 3, "CGC_R", "lgi_neu", out_dir)

# clean up
rm(rcgc_dxGS_neu)
rm(plot_rcgc_dxGS_neu)
rm(rcgc_dxGS_neuOF)
rm(plot_rcgc_dxGS_neuOF)
rm(rcgc_dxGS)
