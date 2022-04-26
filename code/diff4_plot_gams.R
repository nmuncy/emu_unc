library(ggplot2)
library(mgcViz)
library(viridis)
library(itsadug)


draw_global_smooth <- function(plot_obj, attr_num, tract, plot_title, out_dir) {
  # Draw global smooth of AFQ tract.
  #
  # Draw global smooth resulting from GS model. Use plot(sm(obj, int))
  # to extract values, calculate confidence intervals, and ggplot
  # on a scaled X-axis.
  #
  # Arguments:
  #   plot_obj (list) = Plotable object returned by getViz
  #   attr_num (int) = List/attribute number of plot_obj that contains
  #     global smooth
  #   tract (str) = AFQ tract name
  #   plot_title (str) = Plot title for ggplot
  #   out_dir (str) = Path to output location
  #
  # Writes:
  #   <out_dir>/Plot_<tract>_mGS-Global.png

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
    ylab("Est. FA Fit") +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))
  print(pp)

  ggsave(
    paste0(out_dir, "/Plot_", tract, "_mGS-Global.png"),
    plot = last_plot(),
    units = "in",
    width = 4,
    height = 3,
    dpi = 600,
    device = "png"
  )
}


draw_group_smooth <- function(plot_obj, attr_num, tract, plot_title, out_dir) {
  # Draw group smooths of AFQ tract.
  #
  # Draw group smooths resulting from a GS model. Use plot(sm(obj, int))
  # to extract values and ggplot on a scaled Y- and X-axis.
  #
  # Arguments:
  #   plot_obj (list) = Plotable object returned by getViz
  #   attr_num (int) = List/attribute number of plot_obj that contains
  #     group smooths
  #   tract (str) = AFQ tract name
  #   plot_title (str) = Title of plot
  #   out_dir (str) = Path to output location
  #
  # Writes:
  #   <out_dir>/Plot_<tract>_mGS-Group.png

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
    ylab("Est. FA Fit") +
    xlab("Tract Node") +
    theme(text = element_text(family = "Times New Roman"))
  print(pp)

  ggsave(
    paste0(out_dir, "/Plot_", tract, "_mGS-Group.png"),
    plot = last_plot(),
    units = "in",
    width = 4,
    height = 3,
    dpi = 600,
    device = "png"
  )
}


draw_group_smooth_diff <- function(plot_obj, attr_num, tract, plot_title, out_dir) {
  # Draw difference of group smooths of AFQ tract.
  #
  # Plot an A-B difference smooth from an GAM using an ordered factor for group,
  # identify nodes which sig differ from 0, draw polygons to ID.
  #
  # Arguments:
  #   plot_obj (list) = Plotable object returned by getViz
  #   attr_num (int) = List/attribute number of plot_obj that contains group
  #     difference smooth
  #   tract (str) = AFQ tract name
  #   plot_title (str) = Title of plot
  #   out_dir (str) = Path to output location
  #
  # Writes:
  #   <out_dir>/Plot_<tract>_mGS-Diff.png

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
    paste0(out_dir, "/Plot_", tract, "_mGSOF-Diff.png"),
    plot = last_plot(),
    units = "in",
    width = 4,
    height = 3,
    dpi = 600,
    device = "png"
  )
}


draw_two_three <- function(plot_list, name_list, tract, beh_short, out_name){
  # Title
  #
  # Desc
  
  # unpack, organize plots
  r1A <- plot_list$beh$con
  r1B <- plot_list$intx$con
  r2A <- plot_list$beh$exp
  r2B <- plot_list$intx$exp
  r3A <- plot_list$beh$diff
  r3B <- plot_list$intx_diff$diff
  
  # make col titles, y axis, x axis, and row names
  col1_name <- text_grob(name_list$col1, size = 12, family = "Times New Roman")
  col2_name <- text_grob(name_list$col2, size = 12, family = "Times New Roman")
  bot1_name <- text_grob(name_list$bot1, size = 10, family = "Times New Roman")
  bot2_name <- text_grob(name_list$bot2, size = 10, family = "Times New Roman")
  
  l1_name <- l3_name <- ""
  l2_name <-
    text_grob(name_list$rowL, size = 10, family = "Times New Roman", rot = 90)
  
  r1_name <- text_grob(
    name_list$rowR1, size = 12, family = "Times New Roman", rot = 270
  )
  r2_name <- text_grob(
    name_list$rowR2, size = 12, family = "Times New Roman", rot = 270
  )
  r3_name <- text_grob(
    name_list$rowR3, size = 12, family = "Times New Roman", rot = 270
  )
  
  pOut <- grid.arrange(
    arrangeGrob(r1A, top = col1_name, left = l1_name),
    arrangeGrob(r1B, top = col2_name, right = r1_name),
    arrangeGrob(r2A, left = l2_name),
    arrangeGrob(r2B, right = r2_name),
    arrangeGrob(r3A, bottom = bot1_name, left = l3_name),
    arrangeGrob(r3B, bottom = bot2_name, right = r3_name),
    nrow = 3,
    ncol= 2,
    widths = c(0.75, 1), 
    heights = c(1, 1, 1)
  )
  print(pOut)
  
  ggsave(
    paste0(out_dir, "/Plot_", tract, "_", out_name, "_", beh_short, ".png"),
    plot = pOut,
    units = "in",
    height = 6,
    width = 6,
    dpi = 600,
    device = "png"
  )
}
