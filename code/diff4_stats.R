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
  p_data$lb <- p_data$est - p_data$se
  p_data$ub <- p_data$est + p_data$se

  # draw
  tract_long <- switch_names(tract)
  ggplot(data = p_data, aes(x = nodeID, y = est)) +
    geom_line() +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2) +
    scale_x_continuous(breaks = c(seq(0, 99, by = 10), 99)) +
    ggtitle(paste0(tract_long, " Smooth")) +
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
    ggtitle(paste0(tract_long, ", Group Smooths")) +
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

draw_group_diff <- function(plot_obj, attr_num, tract, out_dir) {
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

  # find sig nodes
  p_data$lb <- as.numeric(p_data$y - (2 * p_data$se))
  p_data$ub <- as.numeric(p_data$y + (2 * p_data$se))
  sig_nodes <- which(
    (p_data$y < 0 & p_data$ub < 0) |
      (p_data$y > 0 & p_data$lb > 0)
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
  p <- p +
    annotate(
      "rect",
      xmin = c(d_rect$x_start),
      xmax = c(d_rect$x_end),
      ymin = c(d_rect$y_start),
      ymax = c(d_rect$y_end),
      alpha = 0.2,
      fill = "red"
    ) +
    scale_x_continuous(breaks = c(seq(0, 99, by = 10), 99)) +
    ggtitle(paste0(tract_long, ", Difference Smooth, Experimental-Control")) +
    ylab("Est. Difference") +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))
  print(p)

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
df_afq$dx_group <- factor(df_afq$dx_group)
df_afq$subjectID <- factor(df_afq$subjectID)
df_afq$dx_groupOF <- factor(df_afq$dx_group, ordered = T)


# L. Unc Model Specification ----
#
# 1) Identify distribution
# 2) Determine role of PDS
# 3) Assess GS vs GI fit
# 4) Identify nodes which differ

# subset df_afq, take complete cases
df_tract <- df_afq[which(df_afq$tractID == "UNC_L"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

# 1) determine distribution
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete = F)

# build gam
lunc_gaus <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = gaussian(),
method = "fREML"
)
# gam.check(lunc_gaus, rep = 1000)
summary(lunc_gaus)
# plot(lunc_gaus)


# 2) pds effect
lunc_pds <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50) +
  s(pds, by = sex),
data = df_tract,
family = gaussian(),
method = "fREML"
)
# gam.check(lunc_pds, rep = 1000)
summary(lunc_pds)
# plot(lunc_pds)
compareML(lunc_gaus, lunc_pds) # gaus preferred


# 3) GS model by diagnosis
lunc_dxGS <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2),
data = df_tract,
family = gaussian(),
method = "fREML"
)
# gam.check(lunc_dxGS, rep = 1000)
compareML(lunc_gaus, lunc_dxGS) # lunc_dxGS preferred
summary(lunc_dxGS)
# plot(lunc_dxGS)

# GI model by diagnosis
lunc_dxGI <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(dx_group, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_group, bs = "cr", k = 50, m = 1),
data = df_tract,
family = gaussian(),
method = "fREML"
)
# gam.check(lunc_dxGI, rep = 1000)
summary(lunc_dxGI)
# plot(lunc_dxGI)
compareML(lunc_dxGI, lunc_dxGS) # lunc_dxGS preferred

# draw
plot_lunc_dxGS <- getViz(lunc_dxGS)
# plot(sm(plot_lunc_dxGS, 2))
# plot(sm(plot_lunc_dxGS, 3))
draw_global_smooth(plot_lunc_dxGS, 2, "UNC_L", out_dir)
draw_group_smooth(plot_lunc_dxGS, 3, "UNC_L", out_dir)


# 4) GS dx - identify nodes that differ (make difference smooth)
lunc_dxGS_OF <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_groupOF, bs = "cr", k = 50, m = 2),
data = df_tract,
family = gaussian(),
method = "fREML"
)
# gam.check(lunc_dxGS_OF, rep = 1000)
# plot(lunc_dxGS_OF)
summary(lunc_dxGS_OF) # group diff

# draw
plot_lunc_dxGS_OF <- getViz(lunc_dxGS_OF)
# plot(sm(plot_lunc_dxGS_OF, 2))
# plot(sm(plot_lunc_dxGS_OF, 3))
draw_group_diff(plot_lunc_dxGS_OF, 3, "UNC_L", out_dir)


# L. Unc LGI Interaction ----
#
# 1) Investigate tract-group-negLGI intx
# 2) Investigate tract-group-neuLGI intx

# 1) L. Unc GS neg LGI dx intx
lunc_dxGS_neg <- bam(dti_fa ~ sex +
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
# gam.check(lunc_dxGS_neg, rep = 1000)
compareML(lunc_dxGS, lunc_dxGS_neg) # lunc_dxGS_neg preferred
summary(lunc_dxGS_neg)
# plot(lunc_dxGS_neg)

# draw, unpack tract-LGI intx by group
plot_lunc_dxGS_neg <- getViz(lunc_dxGS_neg)
# plot(sm(plot_lunc_dxGS_neg, 2))
draw_smooth_intx(plot_lunc_dxGS_neg, 2, "UNC_L", "lgi_neg", out_dir)
draw_group_intx(df_tract, lunc_dxGS_neg, "UNC_L", "lgi_neg", out_dir)

# L. Unc GS neg LGI intx, decompose tract curve
lunc_dxGS_neg_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2) +
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
# gam.check(lunc_dxGS_neg_decomp, rep = 1000)
summary(lunc_dxGS_neg_decomp) # group diff in curvature driving intx above?
# plot(lunc_dxGS_neg_decomp)


# 2) L. Unc GS neu LGI dx intx
lunc_dxGS_neu <- bam(dti_fa ~ sex +
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
# gam.check(lunc_dxGS_neu, rep = 1000)
compareML(lunc_dxGS, lunc_dxGS_neu) # lunc_dxGS_neu preferred
summary(lunc_dxGS_neu)
# plot(lunc_dxGS_neu)

# draw, unpack tract-LGI intx by group
plot_lunc_dxGS_neu <- getViz(lunc_dxGS_neu)
# plot(sm(plot_lunc_dxGS_neu, 2)) # lgi_neu linear intx shows nicely
draw_smooth_intx(plot_lunc_dxGS_neu, 2, "UNC_L", "lgi_neu", out_dir)
draw_group_intx(df_tract, lunc_dxGS_neu, "UNC_L", "lgi_neu", out_dir)

# L. Unc GS neu LGI intx, decompose tract curve
lunc_dxGS_neu_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2) +
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
# gam.check(lunc_dxGS_neu_decomp, rep = 1000)
summary(lunc_dxGS_neu_decomp) # strong intx despite controlling for group smooth diff
# plot(lunc_dxGS_neu_decomp)
# plot_lunc_dxGS_neu_decomp <- getViz(lunc_dxGS_neu_decomp)
# plot(sm(plot_lunc_dxGS_neu_decomp, 4))


# R. Unc Model Specification ----
#
# Same as L. Unc Model Specification, but for runc.
#
#   1) Identify distribution
#   2) Determine role of PDS
#   3) Assess GS vs GI fit
#   4) Identify nodes which differ

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
# gam.check(runc_gaus, rep = 1000)
summary(runc_gaus)


# 2) pds effect
runc_pds <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50) +
  s(pds, by = sex),
data = df_tract,
family = gaussian(),
method = "fREML"
)
# gam.check(runc_pds, rep = 1000)
summary(runc_pds)
compareML(runc_gaus, runc_pds) # gaus preferred, again


# 3) GS model by diagnosis
runc_dxGS <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2),
data = df_tract,
family = gaussian(),
method = "fREML"
)
# gam.check(runc_dxGS, rep = 1000)
compareML(runc_gaus, runc_dxGS) # runc_dxGS preferred, again
summary(runc_dxGS)

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
# gam.check(runc_dxGI, rep = 1000)
summary(runc_dxGI)
compareML(runc_dxGI, runc_dxGS) # no real difference, using dxGS

# draw
plot_runc_dxGS <- getViz(runc_dxGS)
draw_global_smooth(plot_runc_dxGS, 2, "UNC_R", out_dir)
draw_group_smooth(plot_runc_dxGS, 3, "UNC_R", out_dir)


# 4) GS dx - identify nodes that differ (make difference smooth)
runc_dxGS_OF <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_groupOF, bs = "cr", k = 50, m = 2),
data = df_tract,
family = gaussian(),
method = "fREML"
)
# gam.check(runc_dxGS_OF, rep = 1000)
summary(runc_dxGS_OF) # group diff

# draw
plot_runc_dxGS_OF <- getViz(runc_dxGS_OF)
draw_group_diff(plot_runc_dxGS_OF, 3, "UNC_R", out_dir)


# R. Unc LGI Interaction ----
#
# Same as L. UNC, but with runc.
#
#   1) Investigate tract-group-negLGI intx
#   2) Investigate tract-group-neuLGI intx

# 1) R. Unc GS neg LGI dx intx
runc_dxGS_neg <- bam(dti_fa ~ sex +
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
# gam.check(runc_dxGS_neg, rep = 1000)
compareML(runc_dxGS, runc_dxGS_neg) # runc_dxGS_neg preferred, again
summary(runc_dxGS_neg)

# draw, unpack tract-LGI intx by group
plot_runc_dxGS_neg <- getViz(runc_dxGS_neg)
draw_smooth_intx(plot_runc_dxGS_neg, 2, "UNC_R", "lgi_neg", out_dir)
draw_group_intx(df_tract, runc_dxGS_neg, "UNC_R", "lgi_neg", out_dir)

# R. Unc GS neg LGI intx, decompose tract curve
runc_dxGS_neg_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2) +
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
# gam.check(runc_dxGS_neg_decomp, rep = 1000)
summary(runc_dxGS_neg_decomp) # still strong lgi_neg intx with nodeID by group


# 2) R. Unc GS neu LGI dx intx
runc_dxGS_neu <- bam(dti_fa ~ sex +
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
# gam.check(runc_dxGS_neu, rep = 1000)
compareML(runc_dxGS, runc_dxGS_neu) # runc_dxGS_neu preferred, again
summary(runc_dxGS_neu)

# draw, unpack tract-LGI intx by group
plot_runc_dxGS_neu <- getViz(runc_dxGS_neu)
draw_smooth_intx(plot_runc_dxGS_neu, 2, "UNC_R", "lgi_neu", out_dir)
draw_group_intx(df_tract, runc_dxGS_neu, "UNC_R", "lgi_neu", out_dir)

# R. Unc GS neu LGI intx, decompose tract curve
runc_dxGS_neu_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2) +
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
# gam.check(runc_dxGS_neu_decomp, rep = 1000)
summary(runc_dxGS_neu_decomp) # no group intx with lgi_neu, main effect of nodeID, lgi_neu


# L. Cing Model Specification ----
#
# Same as L. Unc Model Specification, but for cgc_l.
#
#   1) Identify distribution
#   2) Determine role of PDS
#   3) Assess GS vs GI fit
#   4) Identify nodes which differ

# subset df_afq, take complete cases
df_tract <- df_afq[which(df_afq$tractID == "CGC_L"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

# 1) determine distribution
hist(df_tract$dti_fa) # ugh
descdist(df_tract$dti_fa, discrete = F) # beta or gamma
ggplot(df_tract, aes(y = dti_fa, x = nodeID)) +
  geom_point()

# build gam
lcgc_beta <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = betar(link = "logit"),
method = "fREML"
)
# gam.check(lcgc_beta, rep = 1000)

lcgc_gamma <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
# gam.check(lcgc_gamma, rep = 1000)
compareML(lcgc_beta, lcgc_gamma) # going with gamma
# plot(lcgc_gamma)
summary(lcgc_gamma)


# 2) pds effect
lcgc_pds <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50) +
  s(pds, by = sex),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
# gam.check(lcgc_pds, rep = 1000)
summary(lcgc_pds)
# plot(lcgc_pds)
compareML(lcgc_gamma, lcgc_pds) # lcgc_gamma preferred


# 3) GS model by diagnosis
lcgc_dxGS <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
# gam.check(lcgc_dxGS, rep = 1000)
compareML(lcgc_gamma, lcgc_dxGS) # lcgc_dxGS preferred
summary(lcgc_dxGS)
# plot(lcgc_dxGS)

# GI model by diagnosis
lcgc_dxGI <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(dx_group, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_group, bs = "cr", k = 50, m = 1),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
# gam.check(lcgc_dxGI, rep = 1000)
summary(lcgc_dxGI)
# plot(lcgc_dxGI)
compareML(lcgc_dxGI, lcgc_dxGS) # lcgc_dxGS preferred

# draw
plot_lcgc_dxGS <- getViz(lcgc_dxGS)
draw_global_smooth(plot_lcgc_dxGS, 2, "CGC_L", out_dir)
draw_group_smooth(plot_lcgc_dxGS, 3, "CGC_L", out_dir)

# 4) GS dx - identify nodes that differ (make difference smooth)
lcgc_dxGS_OF <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_groupOF, bs = "cr", k = 50, m = 2),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
# gam.check(lcgc_dxGS_OF, rep = 1000)
plot(lcgc_dxGS_OF)
summary(lcgc_dxGS_OF) # group diff
plot_lcgc_dxGS_OF <- getViz(lcgc_dxGS_OF)
plot(sm(plot_lcgc_dxGS_OF, 2))
plot(sm(plot_lcgc_dxGS_OF, 3))

# draw
draw_group_diff(plot_lcgc_dxGS_OF, 3, "CGC_L", out_dir)


# L. Cing LGI Interaction ----
#
# Same as L. UNC, but with lcgc.
#
#   1) Investigate tract-group-negLGI intx
#   2) Investigate tract-group-neuLGI intx

# 1) L. Cing GS neg LGI dx intx
lcgc_dxGS_neg <- bam(dti_fa ~ sex +
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
# gam.check(lcgc_dxGS_neg, rep = 1000)
compareML(lcgc_dxGS, lcgc_dxGS_neg) # lcgc_dxGS_neg preferred, again
summary(lcgc_dxGS_neg)
plot(lcgc_dxGS_neg)
plot_lcgc_dxGS_neg <- getViz(lcgc_dxGS_neg)
plot(sm(plot_lcgc_dxGS_neg, 2))



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
draw_group_diff(plot_rcgc_dxGS_OF, 3, "CGC_R", out_dir)


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

# R. Cing GS neg LGI intx, decompose tract curve
rcgc_dxGS_neg_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2) +
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
# gam.check(rcgc_dxGS_neg_decomp, rep = 1000)
summary(rcgc_dxGS_neg_decomp) # still strong group intx


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

# R. Cing GS neu LGI intx, decompose tract curve
rcgc_dxGS_neu_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2) +
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
# gam.check(rcgc_dxGS_neu_decomp, rep = 1000)
summary(rcgc_dxGS_neu_decomp) # strong intx still
