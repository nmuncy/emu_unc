library("ggplot2")
library("fitdistrplus")
library("mgcv")
library("itsadug")
library("tidymv")
library("dplyr")
library("mgcViz")
library("tools")
library("viridis")
library("cowplot")


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
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    ggtitle(paste(tract_long, "Smooth")) +
    ylab("Fit Est.") +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))

  ggsave(
    paste0(out_dir, "/Plot_GAM_", tract, "_Global.png"),
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
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    ggtitle(paste(tract_long, "Group Smooths")) +
    ylab("Fit Est.") +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))

  ggsave(
    paste0(out_dir, "/Plot_GAM_", tract, "_GS.png"),
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
  sig_rows <- which(
    (p_data$est < 0 & p_data$ub < 0) |
      (p_data$est > 0 & p_data$lb > 0)
  )
  sig_nodes <- p_data[sig_rows, ]$nodeID

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

  # make df for drawing rectangles
  d_rect <- data.frame(
    x_start = vec_start,
    x_end = vec_end,
    y_start = rep(y_min, length(vec_start)),
    y_end = rep(0, length(vec_start))
  )
  d_rect$x_start <- d_rect$x_start
  d_rect$x_end <- d_rect$x_end

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
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    ggtitle(paste(tract_long, "Exp-Con Difference Smooth")) +
    ylab("Est. Difference") +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))
  # print(p)

  ggsave(
    paste0(out_dir, "/Plot_GAM_", tract, "_Diff.png"),
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
  h_title <- ifelse(
    (y_var == "lgi_neg" | y_var == "lgi_neu"), "Memory Metric", "PPI Term"
  )

  p <- plot(sm(plot_obj, attr_num)) +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    ggtitle(paste0(tract_long, " Node-FA-", h_title, " Smooth")) +
    ylab(beh_long) +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))
  print(p)

  ggsave(
    paste0(out_dir, "/Plot_GAM_", tract, "_Intx_", y_var, ".png"),
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
  h_title <- ifelse(
    (y_var == "lgi_neg" | y_var == "lgi_neu"), "Memory Metric", "PPI Term"
  )

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
    scale_x_continuous(expand = c(0, 0), breaks = c(10, 50, 89)) +
    labs(x = "Tract Node", y = beh_long) +
    ggtitle(paste0(tract_long, " Node-FA-", h_title, " Smooth by Group")) +
    theme(
      legend.position = "right",
      text = element_text(family = "Times New Roman")
    )

  ggsave(
    paste0(out_dir, "/Plot_GAM_", tract, "_Group-Intx_", y_var, ".png"),
    units = "in",
    width = 6,
    height = 6,
    device = "png"
  )
}

draw_group_intx_ref <- function(plot_obj, attr_num, tract, y_var, out_dir) {
  # Draw group interaction reference (control) 3D smooth.
  #
  # Using the output of an ordered-factor group interaction
  # model, draw how the reference group A interacts with continuous
  # metric, nodeID, and predicted FA (s(x)).
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
  #   <out_dir>/Plot_GAM_Group-Intx-Ref_<tract>_<beh>.png

  beh_long <- switch_names(y_var)
  tract_long <- switch_names(tract)
  h_title <- ifelse(
    (y_var == "lgi_neg" | y_var == "lgi_neu"), "Memory Metric", "PPI Term"
  )

  p <- plot(sm(plot_obj, attr_num)) +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    ggtitle(paste0(tract_long, " Control Node-FA-", h_title, " Smooth")) +
    ylab(beh_long) +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))
  print(p)

  ggsave(
    paste0(out_dir, "/Plot_GAM_", tract, "_Group-Intx-Ref_", y_var, ".png"),
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
  h_title <- ifelse(
    (y_var == "lgi_neg" | y_var == "lgi_neu"), "Memory Metric", "PPI Term"
  )

  p <- plot(sm(plot_obj, attr_num)) +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    ggtitle(paste0(
      tract_long, " Experiment Node-FA-", h_title, " Difference Smooth"
    )) +
    ylab(beh_long) +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))
  print(p)

  ggsave(
    paste0(out_dir, "/Plot_GAM_", tract, "_Group-Intx-Diff_", y_var, ".png"),
    units = "in",
    width = 6,
    height = 6,
    device = "png"
  )
}



# Set Up ----

# set paths
proj_dir <- file_path_as_absolute(paste0(getwd(), "/.."))
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


# L. Unc Model Specification ----
#
# 1) Identify distribution
# 2) Determine role of PDS
# 3) Assess GS vs GI fit
# 4) Identify nodes which differ
#
# Steps will be commented for L. Unc models
# but then repeated with basic comments for
# other tracts.

# subset df_afq, take complete cases
df_tract <- df_afq[which(df_afq$tractID == "UNC_L"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

# 1) determine distribution
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete = F)

# build gam with Gamma dist
lunc_gamma <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
gam.check(lunc_gamma, rep = 1000)


# 2) pds effect
lunc_pds <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50) +
  s(pds, by = sex),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)

gam.check(lunc_pds, rep = 1000)
summary(lunc_pds) # no pds effect
compareML(lunc_gamma, lunc_pds) # adding pds does not help
rm(lunc_pds) # keep env trim by removing losing  model


# 3) GS model by diagnosis (dxGS)
lunc_dxGS <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
gam.check(lunc_dxGS, rep = 1000)
compareML(lunc_gamma, lunc_dxGS) # lunc_dxGS preferred
rm(lunc_gamma)

# GI model by diagnosis (dxGI)
lunc_dxGI <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(dx_group, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_group, bs = "cr", k = 50, m = 1),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
gam.check(lunc_dxGI, rep = 1000)
compareML(lunc_dxGI, lunc_dxGS) # no real diff, continue with dxGS
rm(lunc_dxGI)

# get plot object via getViz
plot_lunc_dxGS <- getViz(lunc_dxGS)
plot(sm(plot_lunc_dxGS, 2))
plot(sm(plot_lunc_dxGS, 3))

# make nice plots
draw_global_smooth(plot_lunc_dxGS, 2, "UNC_L", out_dir)
draw_group_smooth(plot_lunc_dxGS, 3, "UNC_L", out_dir)
rm(plot_lunc_dxGS)


# 4) GS dx - identify nodes that differ (make difference smooth)
lunc_dxGS_OF <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_groupOF, bs = "cr", k = 50, m = 2),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
gam.check(lunc_dxGS_OF, rep = 1000)

# get stat of group diff
summary(lunc_dxGS_OF)

# get plot object via getViz
plot_lunc_dxGS_OF <- getViz(lunc_dxGS_OF)
plot(sm(plot_lunc_dxGS_OF, 2))
plot(sm(plot_lunc_dxGS_OF, 3))

# make nice plot, clean up
draw_group_smooth_diff(plot_lunc_dxGS_OF, 3, "UNC_L", out_dir)
rm(lunc_dxGS_OF)
rm(plot_lunc_dxGS_OF)


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
  family = Gamma(link = "logit"),
  method = "fREML"
  )
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
  family = Gamma(link = "logit"),
  method = "fREML"
  )
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


# R. Unc Model Specification ----
#
# Same as L. Unc Model Specification, but for runc.

# subset df_afq, take complete cases
df_tract <- df_afq[which(df_afq$tractID == "UNC_R"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

# 1) determine distribution
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete = F)

# build gam
runc_gaus <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(runc_gaus, rep = 1000)

runc_gamma <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
gam.check(runc_gamma, rep = 1000)
compareML(runc_gaus, runc_gamma) # runc_gaus preferred
rm(runc_gamma)


# 2) pds effect
runc_pds <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50) +
  s(pds, by = sex),
data = df_tract,
family = gaussian(),
method = "fREML"
)

gam.check(runc_pds, rep = 1000)
summary(runc_pds) # no pds by sex effect
compareML(runc_gaus, runc_pds) # gaus preferred
rm(runc_pds)


# 3) GS model by diagnosis
runc_dxGS <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(runc_dxGS, rep = 1000)
compareML(runc_gaus, runc_dxGS) # runc_dxGS preferred
summary(runc_dxGS)
rm(runc_gaus)

# GI model by diagnosis
runc_dxGI <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(dx_group, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_group, bs = "cr", k = 50, m = 1),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(runc_dxGI, rep = 1000)
compareML(runc_dxGI, runc_dxGS) # no real difference, using dxGS
rm(runc_dxGI)

# draw
plot_runc_dxGS <- getViz(runc_dxGS)
draw_global_smooth(plot_runc_dxGS, 2, "UNC_R", out_dir)
draw_group_smooth(plot_runc_dxGS, 3, "UNC_R", out_dir)
rm(plot_runc_dxGS)


# 4) GS dx - identify nodes that differ
runc_dxGS_OF <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_groupOF, bs = "cr", k = 50, m = 2),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(runc_dxGS_OF, rep = 1000)
summary(runc_dxGS_OF) # group diff

# draw
plot_runc_dxGS_OF <- getViz(runc_dxGS_OF)
draw_group_smooth_diff(plot_runc_dxGS_OF, 3, "UNC_R", out_dir)
rm(plot_runc_dxGS_OF)
rm(runc_dxGS_OF)


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


# L. Cing Model Specification ----
#
# Same as L. Unc Model Specification, but for cgc_l.

# subset df_afq, take complete cases
df_tract <- df_afq[which(df_afq$tractID == "CGC_L"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

# 1) determine distribution
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete = F)

# build gam
lcgc_gaus <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lcgc_gaus, rep = 1000)

lcgc_beta <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = betar(link = "logit"),
method = "fREML"
)
gam.check(lcgc_beta, rep = 1000)
compareML(lcgc_beta, lcgc_gaus) # going with gaus
rm(lcgc_beta)

lcgc_gamma <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
gam.check(lcgc_gamma, rep = 1000)
compareML(lcgc_gaus, lcgc_gamma) # going with gaus
rm(lcgc_gamma)


# 2) pds effect
lcgc_pds <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50) +
  s(pds, by = sex),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lcgc_pds, rep = 1000)
summary(lcgc_pds) # no effect of pds
compareML(lcgc_gaus, lcgc_pds) # lcgc_gaus preferred
rm(lcgc_pds)


# 3) GS model by diagnosis
lcgc_dxGS <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lcgc_dxGS, rep = 1000)
compareML(lcgc_gaus, lcgc_dxGS) # lcgc_dxGS preferred
summary(lcgc_dxGS)
rm(lcgc_gaus)

# GI model by diagnosis
lcgc_dxGI <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(dx_group, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_group, bs = "cr", k = 50, m = 1),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lcgc_dxGI, rep = 1000)
compareML(lcgc_dxGI, lcgc_dxGS) # no real diff, continuing with dxGS
rm(lcgc_dxGI)

# draw
plot_lcgc_dxGS <- getViz(lcgc_dxGS)
draw_global_smooth(plot_lcgc_dxGS, 2, "CGC_L", out_dir)
draw_group_smooth(plot_lcgc_dxGS, 3, "CGC_L", out_dir)
rm(plot_lcgc_dxGS)


# 4) GS dx - identify nodes that differ (make difference smooth)
lcgc_dxGS_OF <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_groupOF, bs = "cr", k = 50, m = 2),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lcgc_dxGS_OF, rep = 1000)
summary(lcgc_dxGS_OF) # group diff

# draw
plot_lcgc_dxGS_OF <- getViz(lcgc_dxGS_OF)
draw_group_smooth_diff(plot_lcgc_dxGS_OF, 3, "CGC_L", out_dir)
rm(lcgc_dxGS_OF)
rm(plot_lcgc_dxGS_OF)


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
#
#   1) Identify distribution
#   2) Determine role of PDS
#   3) Assess GS vs GI fit
#   4) Identify nodes which differ

# subset df_afq, take complete cases
df_tract <- df_afq[which(df_afq$tractID == "CGC_R"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

# 1) determine distribution
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete = F) # beta or gamma
ggplot(df_tract, aes(y = dti_fa, x = nodeID)) +
  geom_point()

# build gam
rcgc_beta <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = betar(link = "logit"),
method = "fREML"
)
# gam.check(rcgc_beta, rep = 1000)

rcgc_gamma <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
# gam.check(rcgc_gamma, rep = 1000)
compareML(rcgc_beta, rcgc_gamma) # going with gamma
summary(rcgc_gamma)


# 2) pds effect
rcgc_pds <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50) +
  s(pds, by = sex),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
# gam.check(rcgc_pds, rep = 1000)
summary(rcgc_pds)
compareML(rcgc_gamma, rcgc_pds) # no real diff, continuing with gamma


# 3) GS model by diagnosis
rcgc_dxGS <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
# gam.check(rcgc_dxGS, rep = 1000)
compareML(rcgc_gamma, rcgc_dxGS) # no difference
summary(rcgc_dxGS)

# draw
plot_rcgc_dxGS <- getViz(rcgc_dxGS)
draw_global_smooth(plot_rcgc_dxGS, 2, "CGC_R", out_dir)
draw_group_smooth(plot_rcgc_dxGS, 3, "CGC_R", out_dir)

# GI model by diagnosis
rcgc_dxGI <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(dx_group, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_group, bs = "cr", k = 50, m = 1),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
# gam.check(rcgc_dxGI, rep = 1000)
summary(rcgc_dxGI) # slight group diff, is this model preferred?
compareML(rcgc_dxGI, rcgc_gamma) # no difference - continuing with gamma model

# draw
plot_rcgc_gamma <- getViz(rcgc_gamma)
draw_global_smooth(plot_rcgc_gamma, 2, "CGC_R", out_dir)


# 4) Gamma - identify nodes that differ (make difference smooth)
rcgc_dxGS_OF <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_groupOF, bs = "cr", k = 50, m = 2),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
# gam.check(rcgc_dxGS_OF, rep = 1000)
summary(rcgc_dxGS_OF) # no real group diff
plot_rcgc_dxGS_OF <- getViz(rcgc_dxGS_OF)
draw_group_smooth_diff(plot_rcgc_dxGS_OF, 3, "CGC_R", out_dir)


# R. Cing LGI Interaction ----
#
# Same as L. UNC, but with rcgc.
#
#   1) Investigate tract-group-negLGI intx
#   2) Investigate tract-group-neuLGI intx

# 1) R. Cing GS neg LGI dx intx
rcgc_dxGS_neg <- bam(dti_fa ~ sex +
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
family = Gamma(link = "logit"),
method = "fREML"
)
# gam.check(rcgc_dxGS_neg, rep = 1000)
compareML(rcgc_gamma, rcgc_dxGS_neg) # rcgc_dxGS_neg preferred
summary(rcgc_dxGS_neg)

# draw, unpack tract-LGI intx by group
plot_rcgc_dxGS_neg <- getViz(rcgc_dxGS_neg)
draw_smooth_intx(plot_rcgc_dxGS_neg, 2, "CGC_R", "lgi_neg", out_dir)
draw_group_intx(df_tract, rcgc_dxGS_neg, "CGC_R", "lgi_neg", out_dir)


# 2) R. Cing GS neu LGI dx intx
rcgc_dxGS_neu <- bam(dti_fa ~ sex +
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
# gam.check(rcgc_dxGS_neu, rep = 1000)
compareML(rcgc_gamma, rcgc_dxGS_neu) # rcgc_dxGS_neu preferred
summary(rcgc_dxGS_neu)

# draw, unpack tract-LGI intx by group
plot_rcgc_dxGS_neu <- getViz(rcgc_dxGS_neu)
draw_smooth_intx(plot_rcgc_dxGS_neu, 2, "CGC_R", "lgi_neu", out_dir)
draw_group_intx(df_tract, rcgc_dxGS_neu, "CGC_R", "lgi_neu", out_dir)
