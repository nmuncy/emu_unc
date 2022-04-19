library("fitdistrplus")
library("itsadug")
library("tidymv")
library("dplyr")
library("mgcViz")
library("tools")
library("tidyr")
library("devtools")

# install_local(path = ".", force = T)
library("DiffGamm")


# Functions ----
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
    "amgR" = "L. Amygdala",
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

  # draw plots
  plot_tract_GS <- getViz(tract_GS)
  draw_global_smooth(
    plot_tract_GS,
    2,
    tract,
    paste(switch_names(tract), "Global Smooth"),
    out_dir
  )
  draw_group_smooth(
    plot_tract_GS,
    3,
    tract,
    paste(switch_names(tract), "Group Smooths"),
    out_dir
  )

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


# Interaction with LGI ----
#
# First, model the interaction of group, tract node, and
# a memory metric (negative/neutral LGI) in predicting
# tract FA values. Compare with tract_GS model to determine
# if including LGI improves model fit.
#
# Then, conduct interaction with ordered factors for group
# (ref = Con) to see if experimental group differs in
# interaction from reference group.

beh_list <- c("lgi_neg", "lgi_neu")
for (tract in tract_list) {

  # subset df_afq, keep people w/dx for group modeling
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract <- df_tract %>% drop_na(dx)

  # get distribution, get dxGS model
  tract_dist <- tract_fam(tract)
  tract_GS <- readRDS(paste0(out_dir, "/Model_", tract, "_mGS.Rda"))

  # model intx w/e/behavior
  for (beh in beh_list) {

    # switch name for files
    beh_short <- switch(beh,
      "lgi_neg" = "Neg",
      "lgi_neu" = "Neu"
    )

    # model tract main effect with tract-fa-behavior intx
    gam_file <- paste0(out_dir, "/Model_", tract, "_mGIntx_LGI_", beh_short, ".Rda")
    if (!file.exists(gam_file)) {
      h_gam <- gam_Gintx(df_tract, tract_dist, "dx_group", beh)
      saveRDS(h_gam, file = gam_file)
      rm(h_gam)
    }
    tract_Gintx <- readRDS(gam_file)
    write_gam_stats(tract_Gintx, out_dir, paste0("mGIntx_LGI_", beh_short), tract)
    write_compare_stats(
      tract_GS, tract_Gintx, tract, out_dir, paste0("mGS-mGIntx_LGI_", beh_short)
    )

    # plot
    plot_Gintx <- getViz(tract_Gintx)
    draw_Gintx(
      plot_Gintx,
      2,
      tract,
      switch_names(beh),
      paste(
        switch_names(tract), "Node-FA-Memory Interaction"
      ),
      paste0(out_dir, "/Plot_", tract, "_mGIntx_LGI_", beh_short)
    )

    # model interaction of tract-group-behavior to get interaction
    # smooths for each group
    gam_file <- paste0(out_dir, "/Model_", tract, "_mGSIntx_LGI_", beh_short, ".Rda")
    if (!file.exists(gam_file)) {
      h_gam <- gam_GSintx(df_tract, tract_dist, "dx_group", beh)
      saveRDS(h_gam, file = gam_file)
      rm(h_gam)
    }
    tract_GSintx <- readRDS(gam_file)

    # draw
    plot_tract_GSintx <- getViz(tract_GSintx)

    draw_intx(
      plot_tract_GSintx,
      5,
      tract,
      switch_names(beh),
      paste(
        switch_names(tract), "Node-FA-Memory Interaction, Control"
      ),
      paste0(
        out_dir, "/Plot_", tract, "_mGSIntx_LGI_", beh_short, "_con"
      )
    )

    draw_intx(
      plot_tract_GSintx,
      6,
      tract,
      switch_names(beh),
      paste(
        switch_names(tract), "Node-FA-Memory Interaction, Experimental"
      ),
      paste0(
        out_dir, "/Plot_", tract, "_mGSIntx_LGI_", beh_short, "_exp"
      )
    )

    # test if experiment group interaction differs from control
    gam_file <- paste0(out_dir, "/Model_", tract, "_mGSOFIntx_LGI_", beh_short, ".Rda")
    if (!file.exists(gam_file)) {
      h_gam <- gam_GSintxOF(
        df_tract, tract_dist, "dx_group", "dx_groupOF", beh
      )
      saveRDS(h_gam, file = gam_file)
      rm(h_gam)
    }
    tract_GSintxOF <- readRDS(gam_file)

    # get stats
    write_gam_stats(
      tract_GSintxOF, out_dir, paste0("mGSOFIntx_LGI_", beh_short), tract
    )

    # draw
    plot_tract_GSintxOF <- getViz(tract_GSintxOF)
    plot(sm(plot_tract_GSintxOF, 1))
    plot(sm(plot_tract_GSintxOF, 2))
    plot(sm(plot_tract_GSintxOF, 3))
    plot(sm(plot_tract_GSintxOF, 4))
    plot(sm(plot_tract_GSintxOF, 5))
    plot(sm(plot_tract_GSintxOF, 6))
    draw_intx_diff(
      plot_tract_GSintxOF,
      6,
      tract,
      switch_names(beh),
      paste(
        switch_names(tract),
        "Node-FA-Memory Interaction, Experimental Difference"
      ),
      paste0(
        out_dir, "/Plot_", tract, "_mGSOFIntx_LGI_", beh_short, "_diff"
      )
    )

    # clean up
    rm(tract_GSintx)
    rm(tract_GSintxOF)
    rm(plot_tract_GSintx)
    rm(plot_tract_GSintxOF)
  }
  rm(df_tract)
  rm(tract_GS)
}


# Interaction with decon ROI coefs ----
#
# Desc.

# set seed and behavior lists
roi_list <- c("amgL", "amgR")
beh_list <- c("neg", "neu")

# incorporate ROI values of ses-S1 task-study decon-rVal in df_afq
subj_list <- as.character(unique(df_afq$subjectID))
for (roi in roi_list) {
  df_roi <- read.csv(
    paste0(data_dir, "/df_ses-S1_task-study_decon-rVal_", roi, ".csv")
  )
  for (beh in beh_list) {
    h_col <- paste(roi, beh, sep = "_")
    df_afq[, h_col] <- NA
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

# investigate tract interactions with decon ROI metrics
for (tract in tract_list) {

  # only investigate UNC tracts, match hemispheres
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

  # get tract GS model for comparison
  tract_GS <- readRDS(paste0(out_dir, "/Model_", tract, "_mGS.Rda"))

  # model e/beh separately
  for (beh in beh_list) {
    roi_beh <- paste(roi, beh, sep = "_")
    roi_beh_out <- switch(roi_beh,
      "amgL_neg" = "amgL_Neg",
      "amgL_neu" = "amgL_Neu",
      "amgR_neg" = "amgR_Neg",
      "amgR_neu" = "amgR_Neu"
    )

    # model tract main effect with tract-fa-ROI behavior coef intx
    gam_file <- paste0(out_dir, "/Model_", tract, "_mGIntx_ROI_", roi_beh_out, ".Rda")
    if (!file.exists(gam_file)) {
      h_gam <- gam_Gintx(df_tract, tract_dist, "dx_group", roi_beh)
      saveRDS(h_gam, file = gam_file)
      rm(h_gam)
    }
    tract_Gintx <- readRDS(gam_file)
    write_gam_stats(tract_Gintx, out_dir, paste0("mGIntx_ROI_", roi_beh_out), tract)
    write_compare_stats(
      tract_GS, tract_Gintx, tract, out_dir, paste0("mGS-mGIntx_ROI_", roi_beh_out)
    )

    # plot
    plot_Gintx <- getViz(tract_Gintx)
    draw_Gintx(
      plot_Gintx,
      2,
      tract,
      paste(switch_names(beh), switch_names(roi), "Coefs"),
      paste(
        switch_names(tract), "Node-FA-ROI Interaction"
      ),
      paste0(out_dir, "/Plot_", tract, "_mGIntx_ROI_", roi_beh_out)
    )

    # model interaction of tract-group-behavior
    gam_file <- paste0(out_dir, "/Model_", tract, "_mGSIntx_ROI_", roi_beh_out, ".Rda")
    if (!file.exists(gam_file)) {
      h_gam <- gam_GSintx(df_tract, tract_dist, "dx_group", roi_beh)
      saveRDS(h_gam, file = gam_file)
      rm(h_gam)
    }
    tract_GSintx <- readRDS(gam_file)
    write_gam_stats(tract_GSintx, out_dir, paste0("mGSIntx_ROI_", roi_beh_out), tract)

    # draw
    plot_tract_GSintx <- getViz(tract_GSintx)

    draw_intx(
      plot_tract_GSintx,
      attr_num = 5,
      tract,
      paste(switch_names(beh), switch_names(roi), "Coefs"),
      paste(
        switch_names(tract), "Node-FA-ROI Interaction, Control"
      ),
      paste0(
        out_dir, "/Plot_", tract, "_mGSIntx_ROI_", roi_beh_out, "_con"
      )
    )

    draw_intx(
      plot_obj = plot_tract_GSintx,
      attr_num = 6,
      tract,
      paste(switch_names(beh), switch_names(roi), "Coefs"),
      paste(
        switch_names(tract), "Node-FA-ROI Interaction, Experimental"
      ),
      paste0(
        out_dir, "/Plot_", tract, "_mGSIntx_ROI_", roi_beh_out, "_exp"
      )
    )

    # test if experiment group interaction differs from control
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
    write_gam_stats(tract_GSintxOF, out_dir, paste0("mGSOFIntx_ROI_", roi_beh_out), tract)

    # draw
    plot_tract_GSintxOF <- getViz(tract_GSintxOF)
    draw_intx_diff(
      plot_tract_GSintxOF,
      6,
      tract,
      paste(switch_names(beh), switch_names(roi), "Coefs"),
      paste(
        switch_names(tract), "Node-FA-ROI Interaction, Experimental Difference"
      ),
      paste0(
        out_dir, "/Plot_", tract, "_mGSOFIntx_ROI_", roi_beh_out, "_diff"
      )
    )

    rm(tract_GSintx)
    rm(tract_GSintxOF)
    rm(plot_tract_GSintx)
    rm(plot_tract_GSintxOF)
  }
}



# Interaction with PPI ----
#
# Test for an interaction between group differences in the various tracts
# and the PPI term of the LAmg-ROI for Study trials preceding Test negative
# and neutral Lure FAs.
#
# Distributions of tracts are luckily the same, even w/reduced data,
# so the same syntax as with the LGI section can be used.

# set seed and behavior (+ difference) lists
seed_list <- c("NSlacc", "NSldmpfc", "NSlsfs")
beh_list <- c("SPnegLF", "SPneuLF") # "SPnegLF.SPneuLF"

# incorporate PPI values of ses-S1 task-study in df_afq
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
  rm(df_ppi)
}

# investigate tract interactions with PPI metrics
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

  # get tract GS model for comparison
  tract_GS <- readRDS(paste0(out_dir, "/Model_", tract, "_mGS.Rda"))

  for (seed in seed_list) {
    for (beh in beh_list) {

      # keep subjs who had sufficient number of behaviors to model
      h_seed_beh <- paste(seed, beh, sep = "_")
      df_seed <- df_tract %>% drop_na(h_seed_beh)

      # set up file names
      h_name <- switch(h_seed_beh,
        "NSlacc_SPnegLF" = "LAmg-LACC_NegLF",
        "NSldmpfc_SPnegLF" = "LAmg-LdmPFC_NegLF",
        "NSlsfs_SPnegLF" = "LAmg-LSFS_NegLF",
        "NSlacc_SPneuLF" = "LAmg-LACC_NeuLF",
        "NSldmpfc_SPneuLF" = "LAmg-LdmPFC_NeuLF",
        "NSlsfs_SPneuLF" = "LAmg-LSFS_NeuLF",
        "NSlacc_SPnegLF.SPneuLF" = "LAmg-LACC_NegLF-NeuLF",
        "NSldmpfc_SPnegLF.SPneuLF" = "LAmg-LdmPFC_NegLF-NeuLF",
        "NSlsfs_SPnegLF.SPneuLF" = "LAmg-LSFS_NegLF-NeuLF"
      )

      # model tract main effect with tract-fa-behavior intx
      gam_file <- paste0(out_dir, "/Model_", tract, "_mGIntx_", h_name, ".Rda")
      if (!file.exists(gam_file)) {
        h_gam <- gam_Gintx(df_tract, tract_dist, "dx_group", h_seed_beh)
        saveRDS(h_gam, file = gam_file)
        rm(h_gam)
      }
      tract_Gintx <- readRDS(gam_file)
      write_gam_stats(tract_Gintx, out_dir, paste0("mGSIntx_", h_name), tract)
      write_compare_stats(
        tract_GS, tract_Gintx, tract, out_dir, paste0("mGS-mGSIntx_", h_name)
      )

      # plot
      plot_Gintx <- getViz(tract_Gintx)
      draw_Gintx(
        plot_Gintx,
        2,
        tract,
        paste(h_name, beh, sep = "_"),
        switch_names(h_seed_beh),
        paste(
          switch_names(tract), "Node-FA-PPI Interaction"
        ),
        out_dir
      )

      # # model interaction of tract-group-behavior
      # gam_file <- paste0(
      #   out_dir, "/Model_", tract, "_mGIntx_", h_name, ".Rda"
      # )
      # if (!file.exists(gam_file)) {
      #   h_gam <- gam_intx(df_seed, tract_dist, "dx_group", h_seed_beh)
      #   saveRDS(h_gam, file = gam_file)
      #   rm(h_gam)
      # }
      # tract_intx <- readRDS(gam_file)
      # write_gam_stats(tract_intx, out_dir, paste0("GS-", h_name), tract)
      #
      # # test if experiment group interaction differs from control
      # gam_file <- paste0(
      #   out_dir, "/Model_", tract, "_", h_name, "OF.Rda"
      # )
      # if (!file.exists(gam_file)) {
      #   h_gam <- gam_intxOF(
      #     df_seed, tract_dist, "dx_group", "dx_groupOF", h_seed_beh
      #   )
      #   saveRDS(h_gam, file = gam_file)
      #   rm(h_gam)
      # }
      # tract_intxOF <- readRDS(gam_file)
      # write_gam_stats(tract_intxOF, out_dir, paste0("GSOF-", h_name), tract)
      #
      # # draw
      # plot_tract_intx <- getViz(tract_intx)
      # draw_intx(
      #   plot_obj = plot_tract_intx,
      #   attr_num = 5,
      #   tract,
      #   y_var = paste(h_name, beh, "con", sep = "_"),
      #   y_name = switch_names(h_seed_beh),
      #   plot_title = paste(
      #     switch_names(tract), "Node-FA-PPI Interaction, Control"
      #   ),
      #   out_dir
      # )
      # draw_intx(
      #   plot_obj = plot_tract_intx,
      #   attr_num = 6,
      #   tract,
      #   y_var = paste(h_name, beh, "exp", sep = "_"),
      #   y_name = switch_names(h_seed_beh),
      #   plot_title = paste(
      #     switch_names(tract), "Node-FA-PPI Interaction, Experimental"
      #   ),
      #   out_dir
      # )
      #
      # plot_tract_intxOF <- getViz(tract_intxOF)
      # draw_intx_diff(
      #   plot_obj = plot_tract_intxOF,
      #   attr_num = 5,
      #   tract,
      #   y_var = paste(h_name, beh, sep = "_"),
      #   y_name = switch_names(h_seed_beh),
      #   plot_title = paste(
      #     switch_names(tract),
      #     "Node-FA-PPI Interaction, Experimental Difference"
      #   ),
      #   out_dir
      # )
      #
      # rm(df_seed)
      # rm(tract_intx)
      # rm(tract_intxOF)
      # rm(plot_tract_intx)
      # rm(plot_tract_intxOF)
    }
  }
}
