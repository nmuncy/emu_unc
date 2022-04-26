library("tidymv")
library("dplyr")
library("tools")
library("tidyr")
library("ez")
library("gridExtra")
library("ggpubr")

source("./diff4_calc_gams.R")
source("./diff4_plot_gams.R")
source("./diff4_pred_gams.R")


# General Functions ----
#
# Helper functions for organizing strings and writing files.
# Used by numerous sections below.

tract_fam <- function(tract) {
  # Set family for each tract.
  #
  # These were determined through comparing various models
  # via hist(), fitdistrplus::descdist, and itsadug::compareML.
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

switch_names <- function(name) {
  # Switch tract, Y-axis title names
  #
  # Convert tract names to long form for titles,
  # convert memory, ppi names to long forms for
  # y-axis titles.
  #
  # Arguments:
  #   name (str) = AFQ tract name, behavior name, PPI seed
  #
  # Returns:
  #   name_long (str) = converted name
  x_name <- switch(name,
    "UNC_L" = "L. Uncinate",
    "UNC_R" = "R. Uncinate",
    "CGC_L" = "L. Cingulum",
    "CGC_R" = "R. Cingulum",
    "lgi_neg" = "Negative LGI",
    "lgi_neu" = "Neutral LGI",
    "amgL" = "L. Amygdala",
    "amgR" = "R. Amygdala",
    "neg" = "Negative Scenes",
    "neu" = "Neutral Scenes",
  )
  return(x_name)
}

write_gam_stats <- function(gam_obj, out_dir, gam_type, tract) {
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

write_compare_stats <- function(model_a, model_b, tract, out_dir, out_str) {
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
#
# Set paths and tract list, get session info,
# then read-in data and covert certain columns to factors.
# Finally, clip off 10 nodes from beginning, end.

# set paths
proj_dir <- "/Users/nmuncy/Projects/emu_unc"
data_dir <- paste0(proj_dir, "/data")
out_dir <- paste0(proj_dir, "/stats")
tract_list <- c("UNC_L", "UNC_R", "CGC_L", "CGC_R")

# capture session
capture.output(
  sessionInfo(), file = paste0(proj_dir, "/env/R_session_info.txt")
)

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


# Model Specification ----
#
# Determine the model that best fits the various tracts. Modeling individual
# tracts determined that:
#   a) k=50 was a sufficient basis dimension for all tracts via
#       gam.check(model, rep = 1000)
#   b) the required family arguments via
#       hist(df$dti_fa), fitdistrplus::descdist(df$dti_fa, discrete = F),
#       and itsadug::compareML(model_A, model_B)
#   c) PDS did not increase model fit for any tract (compareML),
#   d) GS fit better than G, GI did not increase fit for all tracts (compareML).
#
# As each GAM for the tracts is very similar, only differing
# in the distribution, we can loop through the tracts.
#
# The dxGS models are saved.
for (tract in tract_list) {

  # subset df_afq, keep people w/dx for group modeling
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract <- df_tract %>% drop_na(dx)

  # conduct basic model
  tract_dist <- tract_fam(tract)
  tract_G <- gam_G(df_tract, tract_dist)
  write_gam_stats(tract_G, out_dir, "mG", tract)

  # does controlling for PDS help model fit
  tract_G_pds <- gam_Gcov(df_tract, tract_dist, "pds")
  write_gam_stats(tract_G_pds, out_dir, "mGPDS", tract)
  write_compare_stats(tract_G, tract_G_pds, tract, out_dir, "mG-mGPDS")

  # test if group smooths increase fit, save model
  gam_file <- paste0(out_dir, "/Model_", tract, "_mGS.Rda")
  if (!file.exists(gam_file)) {
    h_gam <- gam_GS(df_tract, tract_dist, "dx_group")
    saveRDS(h_gam, file = gam_file)
    rm(h_gam)
  }
  tract_GS <- readRDS(gam_file)
  write_gam_stats(tract_GS, out_dir, "mGS", tract)
  write_compare_stats(tract_GS, tract_G, tract, out_dir, "mGS-mG")

  # test if group wiggliness increases fit
  tract_GI <- gam_GI(df_tract, tract_dist, "dx_group")
  write_gam_stats(tract_GS, out_dir, "mGI", tract)
  write_compare_stats(tract_GS, tract_GI, tract, out_dir, "mGS-mGI")

  # test if group smooths differ
  tract_GSOF <- gam_GSOF(df_tract, tract_dist, "dx_groupOF")
  write_gam_stats(tract_GSOF, out_dir, "mGSOF", tract)

  # draw plots for global, group smooths
  plot_tract_GS <- getViz(tract_GS)
  draw_global_smooth(
    plot_tract_GS,
    2,
    tract,
    paste(switch_names(tract), "Tract Smooth"),
    out_dir
  )
  draw_group_smooth(
    plot_tract_GS,
    3,
    tract,
    paste(switch_names(tract), "Group Smooths"),
    out_dir
  )

  # draw diff smooth
  plot_tract_GSOF <- getViz(tract_GSOF)
  draw_group_smooth_diff(
    plot_tract_GSOF,
    3,
    tract,
    paste(switch_names(tract), "Exp-Con Difference Smooth"),
    out_dir
  )

  # clean env
  rm(tract_G)
  rm(tract_G_pds)
  rm(tract_GS)
  rm(tract_GSOF)
  rm(tract_GI)
  rm(plot_tract_GS)
  rm(plot_tract_GSOF)
  rm(df_tract)
}


# switch for selecting node for each tract
tract_node <- function(tract) {
  # Match tract to node, for interaction investiagtion
  #
  # Arguments:
  #   tract (str) = AFQ tract name
  #
  # Returns:
  #   id_node (int) = nodeID number
  id_node <- switch(tract,
   "UNC_L" = 37,
   "UNC_R" = 39,
   "CGC_L" = 57,
   "CGC_R" = 29
  )
  return(id_node)
}


# Interaction with LGI ----
#
# First, test for a group by valence interaction in LGI
# scores.
#
# Next, model the interaction of group, tract node, and
# a memory metric (negative/neutral LGI) in predicting
# tract FA values. Compare with tract_GS model to determine
# if including LGI improves model fit.
#
# Then, conduct interaction with ordered factors for group
# (ref = Con) to see if experimental group differs in
# interaction from reference group.

# set list of lgi behaviors
beh_list <- c("lgi_neg", "lgi_neu")

# subset dataframe, make long format
df_sub <- df_afq[which(df_afq$tractID == "UNC_L" & df_afq$nodeID == 10), ]
df_sub <- df_sub %>% drop_na(dx_group)

df_sub$dx_group <- as.character(df_sub$dx_group)
subj_list <- df_sub$subjectID
group_list <- unique(df_sub$dx_group)
num_beh <- length(beh_list)
num_subj <- length(subj_list)
num_group <- length(group_list)

df_long <- as.data.frame(matrix(NA, nrow = num_subj * num_beh, ncol = 4))
colnames(df_long) <- c("subj", "group", "mem", "value")
df_long$subj <- rep(subj_list, each = num_beh)
df_long$mem <- rep(beh_list, num_subj)

# mine data for e/subject
for (subj in subj_list) {

  # get group
  ind_long <- which(df_long$subj == subj)
  ind_sub <- which(df_sub$subjectID == subj)
  df_long[ind_long, ]$group <- df_sub[ind_sub, ]$dx_group

  # get lgi values
  for (beh in beh_list) {
    ind_long <- which(df_long$subj == subj & df_long$mem == beh)
    df_long[ind_long, ]$value <- df_sub[ind_sub, beh]
  }
}

# test for group x lgi diff
df_long$groupF <- factor(df_long$group)
df_long$memF <- factor(df_long$mem)
fit_anov <- ezANOVA(
  df_long, value,
  wid = subj, within = memF, between = groupF
)
fit_anov # MEs of group, valence; no group:valence

# rename vars for pretty plots
ind_neg <- which(df_long$mem == "lgi_neg")
ind_neu <- which(df_long$mem == "lgi_neu")
df_long[ind_neg, ]$mem <- "Negative"
df_long[ind_neu, ]$mem <- "Neutral"

ggplot(df_long, aes(x = mem, y = value, fill = group)) +
  geom_boxplot() +
  labs(x = "Stimulus Valence", y = "LGI") +
  scale_fill_discrete(name = "Group") +
  ggtitle("Memory Metric by Valence and Group") +
  theme(
    text = element_text(family = "Times New Roman")
  )
# ggsave(
#   "/Users/nmuncy/Desktop/group_lgi.png",
#   plot = last_plot(),
#   units = "in",
#   width = 4,
#   height = 3,
#   dpi = 600,
#   device = "png"
# )
rm(df_long)
rm(df_sub)


# conduct node-fa-lgi intx analyses via GAMs for e/tract
for (tract in tract_list) {

  # subset df_afq, keep people w/dx for group modeling
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract <- df_tract %>% drop_na(dx)

  # get distribution, dxGS model, and node
  tract_dist <- tract_fam(tract)
  tract_GS <- readRDS(paste0(out_dir, "/Model_", tract, "_mGS.Rda"))
  id_node <- tract_node(tract)

  # model intx w/e/behavior
  for (beh in beh_list) {

    # switch name for files
    beh_short <- switch(beh,
      "lgi_neg" = "Neg",
      "lgi_neu" = "Neu"
    )

    # model tract main effect with tract-fa-behavior intx
    gam_file <- paste0(
      out_dir, "/Model_", tract, "_mGIntx_LGI_", beh_short, ".Rda"
    )
    if (!file.exists(gam_file)) {
      h_gam <- gam_Gintx(df_tract, tract_dist, "dx_group", beh)
      saveRDS(h_gam, file = gam_file)
      rm(h_gam)
    }
    tract_Gintx <- readRDS(gam_file)

    # save stats, test against GS model
    write_gam_stats(
      tract_Gintx, out_dir, paste0("mGIntx_LGI_", beh_short), tract
    )
    write_compare_stats(
      tract_GS,
      tract_Gintx,
      tract,
      out_dir,
      paste0("mGS-mGIntx_LGI_", beh_short)
    )

    # model interaction of tract-behavior-fa to get interaction
    # smooths for each group
    gam_file <- paste0(
      out_dir, "/Model_", tract, "_mGSIntx_LGI_", beh_short, ".Rda"
    )
    if (!file.exists(gam_file)) {
      h_gam <- gam_GSintx(df_tract, tract_dist, "dx_group", beh)
      saveRDS(h_gam, file = gam_file)
      rm(h_gam)
    }
    tract_GSintx <- readRDS(gam_file)
    write_gam_stats(
      tract_GSintx, out_dir, paste0("mGSIntx_LGI_", beh_short), tract
    )
    
    plot_group_behs <- pred_group_covs(df_tract, id_node, tract_GSintx, beh)
    plot_group_intx <- pred_group_intx(df_tract, tract_GSintx, beh)
    
    # # sex
    # plot_group_behs_sex <- pred_group_sex_covs(df_tract, id_node, tract_GSintx, beh)
    # plot_group_intx_sex <- pred_group_sex_intx(df_tract, tract_GSintx, beh)
    

    # test if exp group node-fa-lgi intx term differs from control
    gam_file <- paste0(
      out_dir, "/Model_", tract, "_mGSOFIntx_LGI_", beh_short, ".Rda"
    )
    if (!file.exists(gam_file)) {
      h_gam <- gam_GSintxOF(
        df_tract, tract_dist, "dx_group", "dx_groupOF", beh
      )
      saveRDS(h_gam, file = gam_file)
      rm(h_gam)
    }
    tract_GSintxOF <- readRDS(gam_file)

    # get stats, plot
    write_gam_stats(
      tract_GSintxOF, out_dir, paste0("mGSOFIntx_LGI_", beh_short), tract
    )
    plot_group_intx_diff <- pred_group_intx_diff(df_tract, tract_GSintxOF, beh)
    
    # draw grid
    plot_list <- list(
      "beh" = plot_group_behs, 
      "intx" = plot_group_intx, 
      "intx_diff" = plot_group_intx_diff
      )
    name_list <- list(
      "col1" = paste("Node", id_node, "FA-Memory Smooth"),
      "col2" = "Node-FA-Memory Smooth",
      "rowL" = switch_names(beh),
      "rowR1" = "Control",
      "rowR2" = "Experimental",
      "rowR3" = "Difference",
      "bot1" = "Est. FA Fit",
      "bot2" = "Tract Node"
    )
    draw_two_three(plot_list, name_list, tract, beh_short, "LGI")

    # clean up
    rm(tract_Gintx)
    rm(tract_GSintx)
    rm(tract_GSintxOF)
    rm(plot_group_intx)
    rm(plot_group_behs)
    rm(plot_group_intx_diff)
  }
  rm(tract_GS)
}


# Interaction with ROI coefs ----
#
# Incorporate beta-coefficients from left and right amygdala during
# negative and neutral judgments into dataframe. Then test if a group x
# coefs x valence interaction exists for UNC_L-amgL and UNC_R-amgR.

# set seed and behavior lists
roi_list <- c("amgL", "amgR")
beh_list <- c("neg", "neu")

# incorporate ROI values of ses-S1 task-study decon-rVal in df_afq
subj_list <- as.character(unique(df_afq$subjectID))
for (roi in roi_list) {

  # get data
  df_roi <- read.csv(
    paste0(data_dir, "/df_ses-S1_task-study_decon-rVal_", roi, ".csv")
  )

  # work on e/behavior
  for (beh in beh_list) {
    h_col <- paste(roi, beh, sep = "_")
    df_afq[, h_col] <- NA

    # mine e/subj for coefs
    for (subj in subj_list) {
      ind_afq <- which(df_afq$subjectID == subj)
      ind_roi <- which(df_roi$subj == paste0("sub-", subj))
      if (length(ind_roi) == 0) {
        next
      }
      df_afq[ind_afq, h_col] <- df_roi[ind_roi, beh]
    }
  }
  rm(df_roi)
}

# conduct node-fa-roi intx analyses via GAMs for e/tract
for (tract in tract_list) {

  # only investigate UNC tracts, match hemisphere of tract to roi
  if (tract == "CGC_L" || tract == "CGC_R") {
    next
  }
  roi <- switch(tract,
    "UNC_L" = "amgL",
    "UNC_R" = "amgR"
  )

  # subset df_afq, keep people w/dx for group modeling, get dist
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract <- df_tract %>% drop_na(dx)
  tract_dist <- tract_fam(tract)
  id_node <- tract_node(tract)

  # get tract GS model for comparison
  tract_GS <- readRDS(paste0(out_dir, "/Model_", tract, "_mGS.Rda"))

  # model e/beh separately
  for (beh in beh_list) {

    # setup column name, output file name
    roi_beh <- paste(roi, beh, sep = "_")
    roi_beh_out <- switch(roi_beh,
      "amgL_neg" = "amgL_Neg",
      "amgL_neu" = "amgL_Neu",
      "amgR_neg" = "amgR_Neg",
      "amgR_neu" = "amgR_Neu"
    )

    # model tract main effect with tract-fa-ROI behavior coef intx
    gam_file <- paste0(
      out_dir, "/Model_", tract, "_mGIntx_ROI_", roi_beh_out, ".Rda"
    )
    if (!file.exists(gam_file)) {
      h_gam <- gam_Gintx(df_tract, tract_dist, "dx_group", roi_beh)
      saveRDS(h_gam, file = gam_file)
      rm(h_gam)
    }
    tract_Gintx <- readRDS(gam_file)

    # write stats, compare with GS model
    write_gam_stats(
      tract_Gintx, out_dir, paste0("mGIntx_ROI_", roi_beh_out), tract
    )
    write_compare_stats(
      tract_GS,
      tract_Gintx,
      tract,
      out_dir,
      paste0("mGS-mGIntx_ROI_", roi_beh_out)
    )

    # model tract node-fa-roi coef intx for e/group, write stats
    gam_file <- paste0(
      out_dir, "/Model_", tract, "_mGSIntx_ROI_", roi_beh_out, ".Rda"
    )
    if (!file.exists(gam_file)) {
      h_gam <- gam_GSintx(df_tract, tract_dist, "dx_group", roi_beh)
      saveRDS(h_gam, file = gam_file)
      rm(h_gam)
    }
    tract_GSintx <- readRDS(gam_file)
    write_gam_stats(
      tract_GSintx, out_dir, paste0("mGSIntx_ROI_", roi_beh_out), tract
    )
    
    plot_group_behs <- pred_group_covs(df_tract, id_node, tract_GSintx, roi_beh)
    plot_group_intx <- pred_group_intx(df_tract, tract_GSintx, roi_beh)

    # test if exp group node-fa-roi coef intx term differs from control
    gam_file <- paste0(
      out_dir, "/Model_", tract, "_mGSOFIntx_ROI_", roi_beh_out, ".Rda"
    )
    if (!file.exists(gam_file)) {
      h_gam <- gam_GSintxOF(
        df_tract, tract_dist, "dx_group", "dx_groupOF", roi_beh
      )
      saveRDS(h_gam, file = gam_file)
      rm(h_gam)
    }
    tract_GSintxOF <- readRDS(gam_file)
    write_gam_stats(
      tract_GSintxOF, out_dir, paste0("mGSOFIntx_ROI_", roi_beh_out), tract
    )
    
    plot_group_intx_diff <- 
      pred_group_intx_diff(df_tract, tract_GSintxOF, roi_beh)
    
    # draw grid
    plot_list <- list(
      "beh" = plot_group_behs, 
      "intx" = plot_group_intx, 
      "intx_diff" = plot_group_intx_diff
    )
    name_list <- list(
      "col1" = paste("Node", id_node, "FA-ROI Smooth"),
      "col2" = "Node-FA-ROI Smooth",
      "rowL" = paste(switch_names(roi), switch_names(beh)),
      "rowR1" = "Control",
      "rowR2" = "Experimental",
      "rowR3" = "Difference",
      "bot1" = "Est. FA Fit",
      "bot2" = "Tract Node"
    )
    draw_two_three(plot_list, name_list, tract, roi_beh_out, "ROI")
    
    # clean up
    rm(tract_Gintx)
    rm(tract_GSintx)
    rm(tract_GSintxOF)
    rm(plot_group_intx)
    rm(plot_group_behs)
    rm(plot_group_intx_diff)
  }
  rm(tract_GS)
}



# Interaction with PPI coefs ----
#
# Test for an interaction between group differences in the various tracts
# and the PPI term of the LAmg-ROI for Study trials preceding Test negative
# and neutral Lure FAs.

# set seed and behavior lists
seed_list <- c("NSlacc", "NSldmpfc", "NSlsfs")
beh_list <- c("neg", "neu")

# incorporate PPI values of ses-S1 task-study in
# df_afq - same as in ROI section
subj_list <- as.character(unique(df_afq$subjectID))
for (seed in seed_list) {
  df_ppi <- read.csv(
    paste0(data_dir, "/df_ses-S1_task-study_decon-rVal_amgL-", seed, ".csv")
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
  rm(df_ppi)
}

# conduct node-fa-ppi intx analyses via GAMs for e/tract
for (tract in tract_list) {

  # match tract to PPI region, only use L hemi tracts
  if (tract == "CGC_R" || tract == "UNC_R") {
    next
  }
  seed_list <- switch(tract,
    "UNC_L" = "NSlacc",
    "CGC_L" = c("NSlacc", "NSldmpfc", "NSlsfs")
  )

  # subset df_afq, keep people w/dx for group modeling, get dist
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract <- df_tract %>% drop_na(dx)
  tract_dist <- tract_fam(tract)
  id_node <- tract_node(tract)

  # get tract GS model for comparison
  tract_GS <- readRDS(paste0(out_dir, "/Model_", tract, "_mGS.Rda"))

  # models for each seed-beh combination
  for (seed in seed_list) {
    for (beh in beh_list) {

      # keep subjs who had sufficient number of behaviors to model
      h_seed_beh <- paste(seed, beh, sep = "_")
      df_seed <- df_tract %>% drop_na(h_seed_beh)

      # set up file output names
      h_name <- switch(h_seed_beh,
        "NSlacc_neg" = "LAmg-LACC_Neg",
        "NSldmpfc_neg" = "LAmg-LdmPFC_Neg",
        "NSlsfs_neg" = "LAmg-LSFS_Neg",
        "NSlacc_neu" = "LAmg-LACC_Neu",
        "NSldmpfc_neu" = "LAmg-LdmPFC_Neu",
        "NSlsfs_neu" = "LAmg-LSFS_Neu",
      )

      # set up y-axis titles
      y_name <- switch(h_seed_beh,
        "NSlacc_neg" = "LAmg-LACC Coef: Negative Ratings",
        "NSldmpfc_neg" = "LAmg-LdmPFC Coef: Negative Ratings",
        "NSlsfs_neg" = "LAmg-LSFS Coef: Negative Ratings",
        "NSlacc_neu" = "LAmg-LACC Coef: Neutral Ratings",
        "NSldmpfc_neu" = "LAmg-LdmPFC Coef: Neutral Ratings",
        "NSlsfs_neu" = "LAmg-LSFS Coef: Neutral Ratings",
      )

      # model tract main effect with tract-fa-ppi coef intx
      gam_file <- paste0(
        out_dir, "/Model_", tract, "_mGIntx_PPI_", h_name, ".Rda"
      )

      # model tract main effect with tract-fa-behavior intx
      gam_file <- 
        paste0(out_dir, "/Model_", tract, "_mGIntx_PPI_", h_name, ".Rda")
      if (!file.exists(gam_file)) {
        h_gam <- gam_Gintx(df_tract, tract_dist, "dx_group", h_seed_beh)
        saveRDS(h_gam, file = gam_file)
        rm(h_gam)
      }

      # write stats, compare with GS model
      tract_Gintx <- readRDS(gam_file)
      write_gam_stats(
        tract_Gintx, out_dir, paste0("mGSIntx_PPI_", h_name), tract
      )
      write_compare_stats(
        tract_GS,
        tract_Gintx,
        tract,
        out_dir,
        paste0("mGS-mGSIntx_PPI_", h_name)
      )

      # model tract node-fa-ppi coef intx for e/group, write stats
      gam_file <- paste0(
        out_dir, "/Model_", tract, "_mGSIntx_PPI_", h_name, ".Rda"
      )

      # model interaction of tract-group-behavior
      gam_file <- 
        paste0(out_dir, "/Model_", tract, "_mGSIntx_PPI_", h_name, ".Rda")
      if (!file.exists(gam_file)) {
        h_gam <- gam_GSintx(df_tract, tract_dist, "dx_group", h_seed_beh)
        saveRDS(h_gam, file = gam_file)
        rm(h_gam)
      }
      tract_GSintx <- readRDS(gam_file)
      write_gam_stats(
        tract_GSintx, out_dir, paste0("mGSIntx_PPI_", h_name), tract
      )
      
      plot_group_behs <- 
        pred_group_covs(df_tract, id_node, tract_GSintx, h_seed_beh)
      plot_group_intx <- pred_group_intx(df_tract, tract_GSintx, h_seed_beh)

      # test if exp group node-fa-ppi intx term differs from control
      gam_file <- paste0(
        out_dir, "/Model_", tract, "_mGSOFIntx_PPI_", h_name, ".Rda"
      )
      if (!file.exists(gam_file)) {
        h_gam <- gam_GSintxOF(
          df_tract, tract_dist, "dx_group", "dx_groupOF", h_seed_beh
        )
        saveRDS(h_gam, file = gam_file)
        rm(h_gam)
      }
      tract_GSintxOF <- readRDS(gam_file)
      write_gam_stats(
        tract_GSintxOF, out_dir, paste0("mGSOFIntx_PPI_", h_name), tract
      )
      
      plot_group_intx_diff <- 
        pred_group_intx_diff(df_tract, tract_GSintxOF, h_seed_beh)
      
      # draw grid
      plot_list <- list(
        "beh" = plot_group_behs, 
        "intx" = plot_group_intx, 
        "intx_diff" = plot_group_intx_diff
      )
      name_list <- list(
        "col1" = paste("Node", id_node, "FA-PPI Smooth"),
        "col2" = "Node-FA-PPI Smooth",
        "rowL" = y_name,
        "rowR1" = "Control",
        "rowR2" = "Experimental",
        "rowR3" = "Difference",
        "bot1" = "Est. FA Fit",
        "bot2" = "Tract Node"
      )
      draw_two_three(plot_list, name_list, tract, h_name, "PPI")
      
      # clean up
      rm(tract_Gintx)
      rm(tract_GSintx)
      rm(tract_GSintxOF)
      rm(plot_group_intx)
      rm(plot_group_behs)
      rm(plot_group_intx_diff)
    }
  }
  rm(tract_GS)
  rm(df_tract)
}
