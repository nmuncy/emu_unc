library("fitdistrplus")
library("itsadug")
library("tidymv")
library("dplyr")
library("mgcViz")
library("tools")
library("tidyr")
library("devtools")
install_local(path = ".", force = T)
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
  x_name <- switch(name,
    "UNC_L" = "L. Uncinate",
    "UNC_R" = "R. Uncinate",
    "CGC_L" = "L. Cingulum",
    "CGC_R" = "R. Cingulum",
    "lgi_neg" = "Negative LGI",
    "lgi_neu" = "Neutral LGI",
    "NSlacc" = "LAmg-LACC: Study prec. Neg-Neu Lure FA",
    "NSldmpfc" = "LAmg-LdmPFC: Study prec. Neg-Neu Lure FA",
    "NSlsfs" = "LAmg-LSFS: Study prec. Neg-Neu Lure FA",
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


# Functions ----
#
# Switches and plotting functions repeatedly used.

switch_names <- function(name) {
  # Switch for decoding AFQ tract and behavior names
  #
  # Arguments:
  #   name (str) = AFQ tract, behavior name
  #
  # Returns:
  #   x_name (str) = reformatted name

  x_name <- switch(name,
    "UNC_L" = "L. Uncinate",
    "UNC_R" = "R. Uncinate",
    "CGC_L" = "L. Cingulum",
    "CGC_R" = "R. Cingulum",
    "lgi_neg" = "Negative LGI",
    "lgi_neu" = "Neutral LGI",
  )
  return(x_name)
}

draw_global_smooth <- function(plot_obj, attr_num, tract, out_dir) {
  # Draw tract global smooth
  #
  # Arguments:
  #   plot_obj_of (object) = plotable object returned by getViz
  #   attr_num (int) = list/attribute number of plot obj that contains
  #                     group smooths
  #   tract (str) = AFQ tract string
  #   out_dir (str) = path to output location
  #
  # Writes:
  #   <out_dir>/Plot_GAM_Global_<tract>.jpg

  # use plot to extract attribute of interest
  p <- plot(sm(plot_obj, attr_num))
  p_data <- as.data.frame(p$data$fit)
  colnames(p_data) <- c("nodeID", "est", "ty", "se")
  p_data$lb <- as.numeric(p_data$est - (2 * p_data$se))
  p_data$ub <- as.numeric(p_data$est + (2 * p_data$se))

  # draw
  tract_long <- switch_names(tract)
  ggplot(data = p_data, aes(x = nodeID, y = est)) +
    geom_line() +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2) +
    scale_x_continuous(breaks = c(seq(0, 99, by = 10), 99)) +
    ggtitle(paste(tract_long, "Smooth")) +
    ylab("Fit Est.") +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))

  ggsave(
    paste0(out_dir, "/Plot_GAM_Global_", tract, ".png"),
    plot = last_plot(),
    units = "in",
    width = 6,
    height = 6,
    dpi = 600,
    device = "png"
  )
}

draw_group_smooth <- function(plot_obj, attr_num, tract, out_dir) {
  # Draw group smooths.
  #
  # Plot group smooths separate from global smooth.
  #
  # Arguments:
  #   plot_obj_of (object) = plotable object returned by getViz
  #   attr_num (int) = list/attribute number of plot obj that contains
  #                     group smooths
  #   tract (str) = AFQ tract string
  #   out_dir (str) = path to output location
  #
  # Writes:
  #   <out_dir>/Plot_GAM_GS_<tract>.jpg

  # use plot to extract attribute of interest
  p <- plot(sm(plot_obj, attr_num))
  p_data <- as.data.frame(p$data$fit)
  colnames(p_data) <- c("nodeID", "est", "ty", "Group")

  # draw
  tract_long <- switch_names(tract)
  ggplot(data = p_data, aes(x = nodeID, y = est, group = Group)) +
    geom_line(aes(color = Group)) +
    scale_y_continuous(limits = c(-0.2, 0.2)) +
    scale_x_continuous(breaks = c(seq(0, 99, by = 10), 99)) +
    ggtitle(paste(tract_long, "Group Smooths")) +
    ylab("Fit Est.") +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))

  ggsave(
    paste0(out_dir, "/Plot_GAM_GS_", tract, ".png"),
    plot = last_plot(),
    units = "in",
    width = 6,
    height = 6,
    dpi = 600,
    device = "png"
  )
}

draw_group_smooth_diff <- function(plot_obj, attr_num, tract, out_dir) {
  # Draw group difference smooth.
  #
  # Plot an A-B difference smooth, identify nodes
  # which sig differ from 0, draw polygons to ID.
  #
  # Arguments:
  #   plot_obj (object) = plotable object returned by getViz
  #   attr_num (int) = list/attribute number of plot obj that contains
  #                     group difference smooth
  #   tract (str) = AFQ tract string
  #   out_dir (str) = path to output location
  #
  # Writes:
  #   <out_dir>/Plot_GAM_Diff_<tract>.jpg

  # unpack difference smooth data
  p <- plot(sm(plot_obj, attr_num)) +
    geom_hline(yintercept = 0)
  p_data <- as.data.frame(p$data$fit)
  colnames(p_data) <- c("nodeID", "est", "ty", "se")

  # find sig nodes
  p_data$lb <- as.numeric(p_data$est - (2 * p_data$se))
  p_data$ub <- as.numeric(p_data$est + (2 * p_data$se))
  sig_nodes <- which(
    (p_data$est < 0 & p_data$ub < 0) |
      (p_data$est > 0 & p_data$lb > 0)
  )

  # find start, end points of sig regions
  vec_start <- sig_nodes[1]
  vec_end <- vector()
  y_min <- min(p_data$lb)
  num_nodes <- length(sig_nodes)
  c <- 2
  while (c < num_nodes) {
    cc <- c + 1
    if (sig_nodes[cc] > sig_nodes[c] + 1) {
      vec_end <- append(vec_end, sig_nodes[c])
      vec_start <- append(vec_start, sig_nodes[cc])
    }
    c <- cc
  }
  vec_end <- append(vec_end, sig_nodes[num_nodes])

  # make df for drawing rectangles, adjust for 0-index nodeID
  d_rect <- data.frame(
    x_start = vec_start,
    x_end = vec_end,
    y_start = rep(y_min, length(vec_start)),
    y_end = rep(0, length(vec_start))
  )
  d_rect$x_start <- d_rect$x_start - 1
  d_rect$x_end <- d_rect$x_end - 1

  # get tract name
  tract_long <- switch_names(tract)

  # draw
  ggplot(data = p_data, aes(x = nodeID, y = est)) +
    geom_hline(yintercept = 0) +
    geom_line() +
    geom_ribbon(
      aes(ymin = lb, ymax = ub),
      alpha = 0.2
    ) +
    annotate(
      "rect",
      xmin = c(d_rect$x_start),
      xmax = c(d_rect$x_end),
      ymin = c(d_rect$y_start),
      ymax = c(d_rect$y_end),
      alpha = 0.2,
      fill = "red"
    ) +
    scale_x_continuous(breaks = c(seq(0, 99, by = 10))) +
    ggtitle(paste(tract_long, "Exp-Con Difference Smooth")) +
    ylab("Est. Difference") +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))
  # print(p)

  ggsave(
    paste0(out_dir, "/Plot_GAM_Diff_", tract, ".png"),
    plot = last_plot(),
    units = "in",
    width = 6,
    height = 6,
    dpi = 600,
    device = "png"
  )
}

draw_smooth_intx <- function(plot_obj, attr_num, tract, y_var, out_dir) {
  # Draw beavior-nodeID interactions.
  #
  # Arguments:
  #   plot_obj (object) = plotable object returned by getViz
  #   attr_num (int) = list/attribute number of plot obj that contains
  #                     group difference smooth
  #   tract (str) = AFQ tract string
  #   y_var (str) = behavior of interest, used for Y-axis
  #   out_dir (str) = path to output location
  #
  # Writes:
  #   <out_dir>/Plot_GAM_Intx_<tract>_<beh>.png

  beh_long <- switch_names(y_var)
  tract_long <- switch_names(tract)

  p <- plot(sm(plot_obj, attr_num)) +
    scale_x_continuous(breaks = c(seq(0, 99, by = 10), 99)) +
    ggtitle(paste0(tract_long, "-Memory Metric Intx")) +
    ylab(beh_long) +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))
  print(p)

  ggsave(
    paste0(out_dir, "/Plot_GAM_Intx_", tract, "_", y_var, ".png"),
    units = "in",
    width = 6,
    height = 6,
    device = "png"
  )
}

draw_group_intx <- function(df_tract, gam_obj, tract, y_var, out_dir) {
  # Draw behavior-nodeID interactions by group.
  #
  # Plot factorial 3D interaction between nodeID, dti_fa, and behavior
  # as a function of diagnosis group.
  #
  # Arguments:
  #   df_tract (dataframe) = tract dataframe supplied to GAM
  #   gam_obj (object) = returned object from GAM/BAM tool
  #   tract (str) = AFQ tract name
  #   y_var (str) = behavior of interest, used for Y-axis
  #   out_dir (str) = path to output location
  #
  # Writes:
  #   <out_dir>/Plot_GAM_Group-Intx_<tract>_<beh>.png

  df_pred <- transform(
    df_tract,
    h_pred = predict(gam_obj, type = "response")
  )

  beh_long <- switch_names(y_var)
  tract_long <- switch_names(tract)

  ggplot(
    data = df_pred,
    aes(
      x = nodeID,
      y = get(y_var),
      fill = h_pred,
      color = h_pred,
      height = get(y_var)
    )
  ) +
    geom_tile() +
    facet_wrap(~dx_group, ncol = 2) +
    scale_fill_viridis("dti_fa") +
    scale_color_viridis("dti_fa") +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 50, 99)) +
    labs(x = "Tract Node", y = beh_long) +
    ggtitle(paste0(tract_long, "-Memory Metric Intx by Group")) +
    theme(
      legend.position = "right",
      text = element_text(family = "Times New Roman")
    )

  ggsave(
    paste0(out_dir, "/Plot_GAM_Group-Intx_", tract, "_", y_var, ".png"),
    units = "in",
    width = 6,
    height = 6,
    device = "png"
  )
}

draw_group_intx_diff <- function(plot_obj, attr_num, tract, y_var, out_dir) {
  # Draw group interaction difference 3D smooth.
  #
  # Using the output of an ordered-factor group interaction
  # model, draw how group B differs in their nodeID-FA-continuous
  # interaction from the reference group (group A).
  #
  # Arguments:
  #   plot_obj (object) = plotable object returned by getViz
  #   attr_num (int) = list/attribute number of plot obj that contains
  #                     group difference smooth
  #   tract (str) = AFQ tract name
  #   y_var (str) = behavior of interest, used for Y-axis
  #   out_dir (str) = path to output location
  #
  # Writes:
  #   <out_dir>/Plot_GAM_Group-Intx-Diff_<tract>_<beh>.png

  beh_long <- switch_names(y_var)
  tract_long <- switch_names(tract)

  p <- plot(sm(plot_obj, attr_num)) +
    scale_x_continuous(breaks = c(seq(0, 99, by = 10), 99)) +
    ggtitle(paste(tract_long, "Exp-Con Intx Difference Smooth")) +
    ylab(beh_long) +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))
  print(p)

  ggsave(
    paste0(out_dir, "/Plot_GAM_Group-Intx-Diff_", tract, "_", y_var, ".png"),
    units = "in",
    width = 6,
    height = 6,
    device = "png"
  )
}



# Set Up ----

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
  tract_G <- gam_model(df_tract, tract_dist)
  write_gam_stats(tract_G, out_dir, "G", tract)

  # does controlling for PDS help model fit
  tract_G_pds <- gam_cov_model(df_tract, tract_dist, "pds")
  write_gam_stats(tract_G_pds, out_dir, "G-PDS", tract)
  write_compare_stats(tract_G, tract_G_pds, tract, out_dir, "G-PDS")

  # test if group smooths increase fit, save model
  gam_file <- paste0(out_dir, "/Model_", tract, "_dxGS.Rda")
  if (!file.exists(gam_file)) {
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
  draw_global_smooth(
    plot_obj = plot_tract_GS,
    attr_num = 2,
    tract,
    plot_title = paste(switch_names(tract), "Global Smooth"),
    out_dir
  )
  draw_group_smooth(
    plot_obj = plot_tract_GS,
    attr_num = 3,
    tract,
    plot_title = paste(switch_names(tract), "Group Smooths"),
    out_dir
  )

  plot_tract_GSOF <- getViz(tract_GSOF)
  draw_group_smooth_diff(
    plot_obj = plot_tract_GSOF,
    attr_num = 3,
    tract,
    plot_title = paste(switch_names(tract), "Exp-Con Difference Smooth"),
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
#
# Interaction models take a while to generate ...
beh_list <- c("lgi_neg", "lgi_neu")
for (tract in tract_list) {

  # subset df_afq, keep people w/dx for group modeling
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract <- df_tract %>% drop_na(dx)

  # get distribution, get dxGS model
  tract_dist <- tract_fam(tract)
  tract_GS <- readRDS(paste0(out_dir, "/Model_", tract, "_dxGS.Rda"))

  for (beh in beh_list) {

    # model interaction of tract-group-behavior
    gam_file <- paste0(out_dir, "/Model_", tract, "_", beh, ".Rda")
    if (!file.exists(gam_file)) {
      h_gam <- gam_intx_model(df_tract, tract_dist, "dx_group", beh)
      saveRDS(h_gam, file = gam_file)
      rm(h_gam)
    }
    tract_intx <- readRDS(gam_file)
    beh_short <- switch(beh,
      "lgi_neg" = "Neg",
      "lgi_neu" = "Neu"
    )
    write_gam_stats(tract_intx, out_dir, paste0("GS-", beh_short), tract)
    write_compare_stats(
      tract_GS, tract_intx, tract, out_dir, paste0("GS-Intx", beh_short)
    )

    # test if experiment group interaction differs from control
    gam_file <- paste0(out_dir, "/Model_", tract, "_", beh, "_OF.Rda")
    if (!file.exists(gam_file)) {
      h_gam <- gam_intxOF_model(df_tract, tract_dist, "dx_groupOF", beh)
      saveRDS(h_gam, file = gam_file)
      rm(h_gam)
    }
    tract_intxOF <- readRDS(gam_file)
    write_gam_stats(tract_intxOF, out_dir, paste0("GSOF-", beh_short), tract)

    # draw
    plot_tract_intx <- getViz(tract_intx)
    draw_smooth_intx(
      plot_obj = plot_tract_intx,
      attr_num = 2,
      tract,
      y_var = beh,
      y_name = switch_names(beh),
      plot_title = paste(
        switch_names(tract), "Node-FA-Memory Interaction"
      ),
      out_dir
    )
    draw_group_intx(
      df = df_tract,
      gam_obj = tract_intx,
      tract,
      y_var = beh,
      y_name = switch_names(beh),
      plot_title = paste(
        switch_names(tract), "Node-FA-Memory Interaction, by Group"
      ),
      out_dir
    )

    plot_tract_intxOF <- getViz(tract_intxOF)
    draw_group_intx_ref(
      plot_obj = plot_tract_intxOF,
      attr_num = 2,
      tract,
      y_var = beh,
      y_name = switch_names(beh),
      plot_title = paste(
        switch_names(tract), "Node-FA-Memory Interaction, Control"
      ),
      out_dir
    )
    draw_group_intx_diff(
      plot_obj = plot_tract_intxOF,
      attr_num = 3,
      tract,
      y_var = beh,
      y_name = switch_names(beh),
      plot_title = paste(
        switch_names(tract),
        "Node-FA-Memory Interaction, Experimental Difference"
      ),
      out_dir
    )

    # clean up
    rm(tract_intx)
    rm(tract_intxOF)
    rm(plot_tract_intx)
    rm(plot_tract_intxOF)
  }
  rm(df_tract)
  rm(tract_GS)
}


# Interaction with PPI ----
#
# Dist are same, even w/reduced data

# incorporate PPI values of ses-S1 task-study in df_afq
seed_list <- c("NSlacc", "NSldmpfc", "NSlsfs")
beh <- "SPnegLF.SPneuLF"
subj_list <- as.character(unique(df_afq$subjectID))
for (seed in seed_list) {
  df_ppi <- read.csv(
    paste0(data_dir, "/df_ses-S1_task-study_amgL-", seed, ".csv")
  )
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
  rm(df_ppi)
}

for (tract in tract_list) {

  # match tract to PPI region
  if (tract == "CGC_R" || tract == "UNC_R") {
    next
  }
  seed_list <- switch(tract,
    "UNC_L" = "NSlacc",
    "CGC_L" = c("NSldmpfc", "NSlsfs")
  )

  # subset df_afq, keep people w/dx for group modeling, get dist
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract <- df_tract %>% drop_na(dx)
  tract_dist <- tract_fam(tract)

  for (seed in seed_list) {

    # keep subjs who had sufficient number of behaviors to model
    h_seed_beh <- paste(seed, beh, sep = "_")
    df_seed <- df_tract %>% drop_na(h_seed_beh)

    # set up file names
    h_name <- switch(h_seed_beh,
      "NSlacc_SPnegLF.SPneuLF" = "LAmg-LACC_NegLF-NeuLF",
      "NSldmpfc_SPnegLF.SPneuLF" = "LAmg-LdmPFC_NegLF-NeuLF",
      "NSlsfs_SPnegLF.SPneuLF" = "LAmg-LSFS_NegLF-NeuLF"
    )

    # model interaction of tract-group-behavior
    gam_file <- paste0(
      out_dir, "/Model_", tract, "_", h_name, ".Rda"
    )
    if (!file.exists(gam_file)) {
      h_gam <- gam_intx_model(df_seed, tract_dist, "dx_group", h_seed_beh)
      saveRDS(h_gam, file = gam_file)
      rm(h_gam)
    }
    tract_intx <- readRDS(gam_file)
    write_gam_stats(tract_intx, out_dir, paste0("GS-", h_name), tract)

    # test if experiment group interaction differs from control
    gam_file <- paste0(
      out_dir, "/Model_", tract, "_", h_name, "OF.Rda"
    )
    if (!file.exists(gam_file)) {
      h_gam <- gam_intxOF_model(df_seed, tract_dist, "dx_groupOF", h_seed_beh)
      saveRDS(h_gam, file = gam_file)
      rm(h_gam)
    }
    tract_intxOF <- readRDS(gam_file)
    write_gam_stats(tract_intxOF, out_dir, paste0("GSOF-", h_name), tract)

    # draw
    plot_tract_intx <- getViz(tract_intx)
    draw_smooth_intx(
      plot_obj = plot_tract_intx,
      attr_num = 2,
      tract,
      y_var = h_seed_beh,
      y_name = switch_names(seed),
      plot_title = paste(
        switch_names(tract), "Node-FA-PPI Interaction"
      ),
      out_dir
    )
    draw_group_intx(
      df = df_seed,
      gam_obj = tract_intx,
      tract,
      y_var = h_seed_beh,
      y_name = switch_names(seed),
      plot_title = paste(
        switch_names(tract), "Node-FA-PPI Interaction, by Group"
      ),
      out_dir
    )

    plot_tract_intxOF <- getViz(tract_intxOF)
    draw_group_intx_ref(
      plot_obj = plot_tract_intxOF,
      attr_num = 2,
      tract,
      y_var = h_seed_beh,
      y_name = switch_names(seed),
      plot_title = paste(
        switch_names(tract), "Node-FA-PPI Interaction, Control"
      ),
      out_dir
    )
    draw_group_intx_diff(
      plot_obj = plot_tract_intxOF,
      attr_num = 3,
      tract,
      y_var = h_seed_beh,
      y_name = switch_names(seed),
      plot_title = paste(
        switch_names(tract), "Node-FA-PPI Interaction, Experimental Difference"
      ),
      out_dir
    )

    rm(df_seed)
    rm(tract_intx)
    rm(tract_intxOF)
    rm(plot_tract_intx)
    rm(plot_tract_intxOF)
  }
}
