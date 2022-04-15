
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
    paste0(out_dir, "/Plot_GAM_", tract, "_mGS-Global.png"),
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
    paste0(out_dir, "/Plot_GAM_", tract, "_mGS-Group.png"),
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
#' Plot an A-B difference smooth from an GAM using an ordered factor for group,
#' identify nodes which sig differ from 0, draw polygons to ID.
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
    paste0(out_dir, "/Plot_GAM_", tract, "_mGSOF-Diff.png"),
    plot = last_plot(),
    units = "in",
    width = 6,
    height = 6,
    dpi = 600,
    device = "png"
  )
}

#' Plot interaction of global and covariate smooth.
#'
#' Draw the overall interaction between the tract smooth and a covariate.
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
#' @details Writes <out_dir>/Plot_GAM_<tract>_Intx_<x_var>.png
#' @import ggplot2
#' @import mgcViz
draw_Gintx <- function(plot_obj, attr_num, tract, y_var, y_name, plot_title, out_dir) {

  p <- plot(sm(plot_obj, attr_num))
  p_data <- p$data$fit
  colnames(p_data) <- c("z", "tz", "node", "cov", "se")

  pp <- ggplot(
    data = p_data, aes(x = .data$node, y = .data$cov, z = .data$z)
  ) +
    geom_tile(aes(fill = .data$z)) +
    geom_contour(colour = "black") +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    scale_fill_viridis(option = "D", name = "Est. FA Fit") +
    labs(y = y_name, x = "Tract Node") +
    ggtitle(plot_title) +
    theme(text = element_text(family = "Times New Roman"))
  print(pp)

  ggsave(
    paste0(out_dir, "/Plot_GAM_", tract, "_mGSIntx_", y_var, ".png"),
    units = "in",
    width = 6,
    height = 6,
    device = "png"
  )
}

#' Draw interaction smooths of AFQ tract node, FA, and behavioral covariate.
#'
#' Plot using the gamViz plottable object and corresponding attribute. Draws
#' tract node as the X-axis.
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
#' @details Writes <out_dir>/Plot_GAM_<tract>_Intx_<x_var>.png
#' @import ggplot2
#' @import mgcViz
draw_intx <- function(plot_obj, attr_num, tract, y_var, y_name, plot_title, out_dir) {

  # switch x-y so node is X-axis
  p <- plot(sm(plot_obj, attr_num))
  p_data <- p$data$fit
  colnames(p_data) <- c("z", "tz", "cov", "node", "se")

  pp <- ggplot(
    data = p_data, aes(x = .data$node, y = .data$cov, z = .data$z)
  ) +
    geom_tile(aes(fill = .data$z)) +
    geom_contour(colour = "black") +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    scale_fill_viridis(option = "D", name = "Est. FA Fit") +
    labs(y = y_name, x = "Tract Node") +
    ggtitle(plot_title) +
    theme(text = element_text(family = "Times New Roman"))
  print(pp)

  ggsave(
    paste0(out_dir, "/Plot_GAM_", tract, "_mGSIntx_", y_var, ".png"),
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
#' interaction from the reference group (group A). The plot is inverted
#' to illustrate B-A.
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
#' @import viridis
draw_intx_diff <- function(plot_obj, attr_num, tract, y_var, y_name, plot_title, out_dir) {

  # invert direction for ease of interpretation, switch x-y so node is X-axis
  p <- plot(sm(plot_obj, attr_num))
  p_data <- p$data$fit
  p_data$zI <- -1 * p_data$z
  colnames(p_data) <- c("z", "tz", "cov", "node", "se", "zI")

  pp <- ggplot(
    data = p_data, aes(x = .data$node, y = .data$cov, z = .data$zI)
  ) +
    geom_tile(aes(fill = .data$zI)) +
    geom_contour(colour = "black") +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    scale_fill_viridis(option = "D", name = "Est. FA Fit") +
    labs(y = y_name, x = "Tract Node") +
    ggtitle(plot_title) +
    theme(text = element_text(family = "Times New Roman"))
  print(pp)

  ggsave(
    paste0(out_dir, "/Plot_GAM_", tract, "_mGSOFIntx-Diff_", y_var, ".png"),
    units = "in",
    width = 6,
    height = 6,
    device = "png"
  )
}
