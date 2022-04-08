
#' Draw global smooth of AFQ tract.
#'
#' Draw global smooth resulting from GS model. Use plot(sm(obj, int))
#' to extract values, calculate confidence intervals, and ggplot
#' on a scaled X-axis.
#'
#' @export
#' @param plot_obj Plotable object returned by getViz (object)
#' @param attr_num List/attribute number of plot_obj that contains
#' global smooth (int)
#' @param tract AFQ tract name (str)
#' @param plot_title Title of plot (str)
#' @param out_dir Path to output location
#' @details Writes <out_dir>/Plot_GAM_<tract>_Global.png
#' @import ggplot2
#' @import mgcViz
draw_global_smooth <- function(plot_obj, attr_num, tract, plot_title, out_dir) {

  # use plot to extract attribute of interest
  p <- plot(sm(plot_obj, attr_num))
  p_data <- as.data.frame(p$data$fit)
  colnames(p_data) <- c("nodeID", "est", "ty", "se")
  p_data$lb <- as.numeric(p_data$est - (2 * p_data$se))
  p_data$ub <- as.numeric(p_data$est + (2 * p_data$se))

  # make, save ggplot
  pp <- ggplot(data = p_data, aes(x = .data$nodeID, y = .data$est)) +
    geom_line() +
    geom_ribbon(aes(ymin = .data$lb, ymax = .data$ub), alpha = 0.2) +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    ggtitle(plot_title) +
    ylab("Fit Est.") +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))
  print(pp)

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

#' Draw group smooths of AFQ tract.
#'
#' Draw group smooths resulting from a GS model. Use plot(sm(obj, int))
#' to extract values and ggplot on a scaled Y- and X-axis.
#'
#' @export
#' @param plot_obj Plotable object returned by getViz (object)
#' @param attr_num List/attribute number of plot_obj that contains
#' group smooths (int)
#' @param tract AFQ tract name (str)
#' @param plot_title Title of plot (str)
#' @param out_dir Path to output location
#' @details Writes <out_dir>/Plot_GAM_<tract>_GS.png
#' @import ggplot2
#' @import mgcViz
draw_group_smooth <- function(plot_obj, attr_num, tract, plot_title, out_dir) {

  # use plot to extract attribute of interest
  p <- plot(sm(plot_obj, attr_num))
  p_data <- as.data.frame(p$data$fit)
  colnames(p_data) <- c("nodeID", "est", "ty", "Group")

  # make, save ggplot
  pp <- ggplot(
    data = p_data,
    aes(x = .data$nodeID, y = .data$est, group = .data$Group)
  ) +
    geom_line(aes(color = .data$Group)) +
    scale_y_continuous(limits = c(-0.2, 0.2)) +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    ggtitle(plot_title) +
    ylab("Fit Est.") +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))
  print(pp)

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

#' Draw difference of group smooths of AFQ tract.
#'
#' Plot an A-B difference smooth, identify nodes which sig differ from 0,
#' draw polygons to ID.
#'
#' @export
#' @param plot_obj Plotable object returned by getViz (object)
#' @param attr_num List/attribute number of plot_obj that contains group
#' difference smooth (int)
#' @param tract AFQ tract name (str)
#' @param plot_title Title of plot (str)
#' @param out_dir Path to output location
#' @details Writes <out_dir>/Plot_GAM_<tract>_GS-Diff.png
#' @import ggplot2
#' @import mgcViz
draw_group_smooth_diff <- function(plot_obj, attr_num, tract, plot_title, out_dir) {

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

  # draw
  pp <- ggplot(data = p_data, aes(x = .data$nodeID, y = .data$est)) +
    geom_hline(yintercept = 0) +
    geom_line() +
    geom_ribbon(
      aes(ymin = .data$lb, ymax = .data$ub),
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
    ggtitle(plot_title) +
    ylab("Est. Difference") +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))
  print(pp)

  ggsave(
    paste0(out_dir, "/Plot_GAM_", tract, "_GS-Diff.png"),
    plot = last_plot(),
    units = "in",
    width = 6,
    height = 6,
    dpi = 600,
    device = "png"
  )
}

#' Draw behavior-nodeID interactions.
#'
#' @export
#' @param plot_obj Plotable object returned by getViz (object)
#' @param attr_num List/attribute number of plot_obj that contains interaction
#'  smooth (int)
#' @param tract AFQ tract name (str)
#' @param y_var Behavior of interest, used for Y-axis (str)
#' @param y_name Y-axis title (str)
#' @param plot_title Title of plot (str)
#' @param out_dir Path to output location
#' @details Writes <out_dir>/Plot_GAM_<tract>_Intx_<y_var>.png
#' @import ggplot2
#' @import mgcViz
draw_smooth_intx <- function(plot_obj, attr_num, tract, y_var, y_name, plot_title, out_dir) {
  p <- plot(sm(plot_obj, attr_num)) +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    ggtitle(plot_title) +
    ylab(y_name) +
    xlab("Tract Node") +
    guides(fill = guide_legend(title = "Est. FA Fit", reverse = TRUE)) +
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

#' Draw behavior-nodeID interactions by group.
#'
#' Plot factorial 3D interaction between nodeID, dti_fa, and behavior
#' as a function of diagnosis group.
#'
#' @export
#' @param df Tract dataframe supplied to GAM
#' @param gam_obj Returned object from GAM/BAM tool
#' @param tract AFQ tract name (str)
#' @param y_var Column name of behavior of interest (str)
#' @param y_name Y-axis title (str)
#' @param plot_title Title for plot (str)
#' @param out_dir Path to output location (str)
#' @details Writes <out_dir>/Plot_GAM_<tract>_Group-Intx_<y_var>.png
#' @import mgcViz
#' @import ggplot2
#' @import viridis
#' @importFrom stats predict
draw_group_intx <- function(df, gam_obj, tract, y_var, y_name, plot_title, out_dir) {

  # get model predictions
  df_pred <- transform(
    df,
    h_pred = predict(gam_obj, type = "response")
  )
  colnames(df_pred)[which(names(df_pred) == "dti_fa")] <- "Est. FA Fit"

  # draw plot and save
  p <- ggplot(
    data = df_pred,
    aes(
      x = .data$nodeID,
      y = get(y_var),
      fill = .data$h_pred,
      color = .data$h_pred,
      height = get(y_var)
    )
  ) +
    geom_tile() +
    facet_wrap(~ .data$dx_group, ncol = 2) +
    scale_fill_viridis("Est. FA Fit") +
    scale_color_viridis("Est. FA Fit") +
    scale_x_continuous(expand = c(0, 0), breaks = c(10, 50, 89)) +
    labs(x = "Tract Node", y = y_name) +
    ggtitle(plot_title) +
    theme(
      legend.position = "right",
      text = element_text(family = "Times New Roman")
    )
  print(p)

  ggsave(
    paste0(out_dir, "/Plot_GAM_", tract, "_Group-Intx_", y_var, ".png"),
    units = "in",
    width = 6,
    height = 6,
    device = "png"
  )
}

#' Draw 3D smooth for reference (control) group interaction.
#'
#' Using the output of an ordered-factor group interaction
#' model, draw how the reference group A interacts with continuous
#' metric, nodeID, and predicted FA (s(x)).
#'
#' @export
#' @param plot_obj Plotable object returned by getViz (object)
#' @param attr_num List/attribute number of plot_obj that contains
#' reference group interaction smooth of reference group (int)
#' @param tract AFQ tract name (str)
#' @param y_var Behavior of interest, used for Y-axis (str)
#' @param y_name Y-axis title (str)
#' @param plot_title Title of plot (str)
#' @param out_dir Path to output location
#' @details Writes <out_dir>/Plot_GAM_<tract>_Group-Intx-Ref_<y_var>.png
#' @import ggplot2
#' @import mgcViz
draw_group_intx_ref <- function(plot_obj, attr_num, tract, y_var, y_name, plot_title, out_dir) {
  p <- plot(sm(plot_obj, attr_num)) +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    ggtitle(plot_title) +
    ylab(y_name) +
    xlab("Tract Node") +
    guides(fill = guide_legend(title = "Est. FA Fit", reverse = TRUE)) +
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

#' Draw 3D interaction difference smooth for experimental group.
#'
#' Using the output of an ordered-factor group interaction
#' model, draw how group B differs in their nodeID-FA-continuous
#' interaction from the reference group (group A).
#'
#' @export
#' @param plot_obj Plotable object returned by getViz (object)
#' @param attr_num List/attribute number of plot_obj that contains
#' reference group difference interaction smooth (int)
#' @param tract AFQ tract name (str)
#' @param y_var Behavior of interest, used for Y-axis (str)
#' @param y_name Y-axis title (str)
#' @param plot_title Title of plot (str)
#' @param out_dir Path to output location
#' @details Writes <out_dir>/Plot_GAM_<tract>_Group-Intx-Diff_<y_var>.png
#' @import ggplot2
#' @import mgcViz
draw_group_intx_diff <- function(plot_obj, attr_num, tract, y_var, y_name, plot_title, out_dir) {
  p <- plot(sm(plot_obj, attr_num)) +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    ggtitle(plot_title) +
    ylab(y_name) +
    xlab("Tract Node") +
    guides(fill = guide_legend(title = "Est. FA Fit", reverse = TRUE)) +
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
