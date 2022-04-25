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


# draw_Gintx <- function(plot_obj, attr_num, tract, y_name, plot_title, out_file) {
#   # Plot interaction of global and covariate smooth.
#   #
#   # Draw the overall interaction between the tract smooth and a covariate.
#   #
#   # Arguments:
#   #   plot_obj (list) = Plotable object returned by getViz
#   #   attr_num (int) = List/attribute number of plot_obj that contains
#   #     interaction smooth
#   #   tract (str) = AFQ tract name
#   #   y_name (str) = Y-axis title
#   #   plot_title (str) = Title of plot
#   #   out_file (str) = Path to output location, file name
#   #
#   # Writes:
#   #   <out_file>.png
# 
#   p <- plot(sm(plot_obj, attr_num))
#   p_data <- p$data$fit
#   colnames(p_data) <- c("z", "tz", "node", "cov", "se")
# 
#   pp <- ggplot(
#     data = p_data, aes(x = .data$node, y = .data$cov, z = .data$z)
#   ) +
#     geom_tile(aes(fill = .data$z)) +
#     geom_contour(colour = "black") +
#     scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
#     scale_fill_viridis(option = "D", name = "Est. FA Fit") +
#     labs(y = y_name, x = "Tract Node") +
#     ggtitle(plot_title) +
#     theme(text = element_text(family = "Times New Roman"))
#   print(pp)
# 
#   ggsave(
#     paste0(out_file, ".png"),
#     units = "in",
#     width = 4,
#     height = 3,
#     dpi = 600,
#     device = "png"
#   )
# }


# pred_group_intx <- function(){
#   # Title.
#   #
#   # Desc.
#   
#   
# }


draw_intx <- function(plot_obj, attr_num, tract, y_name, plot_title, out_file) {
  # Draw interaction smooths of AFQ tract node, FA, and behavioral covariate.
  #
  # Plot using the gamViz plottable object and corresponding attribute. Draws
  # tract node as the X-axis.
  #
  # Arguments:
  #   plot_obj (list) = Plotable object returned by getViz
  #   attr_num (int) = List/attribute number of plot_obj that contains
  #     interaction smooth
  #   tract (str) = AFQ tract name
  #   y_name (str) = Y-axis title
  #   plot_title (str) = Title of plot
  #   out_file (str) = Path to output location, file name
  #
  # Writes:
  #   <out_file>.png

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
    paste0(out_file, ".png"),
    units = "in",
    width = 4,
    height = 3,
    dpi = 600,
    device = "png"
  )
}

draw_covS_diff <- function(gam_obj, x_name, plot_title, out_file) {
  # Draw difference smooth of covariate.
  #
  # Use plot_diff to calculate group difference smooth for covariate
  # from gam_GSintxOF model. Unpack, determine regions which differ,
  # and plot the smooth with sig boxes.
  #
  # Arguments:
  #   gam_obj (gam) = GAM object returned by mgcv, specifically the
  #     function gam_GSintxOF
  #   x_name (str) = name of X-axis label
  #   plot_title (str) = Title of plot
  #   out_file (str) = Path to output location, file name
  #
  # Writes:
  #   <out_file>.png

  # calc cov group diff smooth
  p_data <- plot_diff(
    gam_obj,
    view = c("h_var"),
    comp = list(h_group = c("Con", "Exp")),
    rm.ranef = T,
    plot = F
  )

  # determine regions that differ from zero
  p_data$lb <- p_data$est - p_data$CI
  p_data$ub <- p_data$est + p_data$CI
  sig_rows <- which(
    (p_data$est < 0 & p_data$ub < 0) |
      (p_data$est > 0 & p_data$lb > 0)
  )
  sig_nodes <- p_data[sig_rows, ]$h_var

  # find start, end points of sig regions, deal w/no sig regions
  vec_start <- sig_nodes[1]
  if(!is.na(vec_start)){
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
  }

  # draw smooth, shade diff regions if they exist
  if(!is.na(vec_start)){
    pp <- ggplot(data = p_data, aes(x = .data$h_var, y = .data$est)) +
      geom_hline(yintercept = 0) +
      geom_line() +
      geom_ribbon(aes(ymin = .data$lb, ymax = .data$ub), alpha = 0.2) +
      annotate(
        "rect",
        xmin = c(d_rect$x_start),
        xmax = c(d_rect$x_end),
        ymin = c(d_rect$y_start),
        ymax = c(d_rect$y_end),
        alpha = 0.2,
        fill = "red"
      ) +
      labs(x = x_name, y = "Est. FA Fit") +
      ggtitle(plot_title) +
      theme(
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(size = 12)
      )
  } else {
    pp <- ggplot(data = p_data, aes(x = .data$h_var, y = .data$est)) +
      geom_hline(yintercept = 0) +
      geom_line() +
      geom_ribbon(aes(ymin = .data$lb, ymax = .data$ub), alpha = 0.2) +
      labs(x = x_name, y = "Est. FA Fit") +
      ggtitle(plot_title) +
      theme(
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(size = 12)
      )
  }
  print(pp)
  ggsave(
    filename = paste0(out_file, ".png"),
    plot = last_plot(),
    units = "in",
    width = 4,
    height = 3,
    dpi = 600,
    device = "png"
  )
}

draw_intx_diff <- function(plot_obj, attr_num, tract, y_name, plot_title, out_file) {
  # Draw 3D interaction difference smooth for experimental group.
  #
  # Using the output of an ordered-factor group interaction
  # model, draw how group B differs in their nodeID-FA-continuous
  # interaction from the reference group (group A). The plot is inverted
  # to illustrate B-A.
  #'
  # Arguments:
  #   plot_obj (list) = Plotable object returned by getViz
  #   attr_num  (int) = List/attribute number of plot_obj that contains
  #     reference group difference interaction smooth (int)
  #   tract (str) = AFQ tract name
  #   y_name (str) = Y-axis title
  #   plot_title (str) = Title of plot
  #   out_file (str) = Path to output location, file name
  #
  # Writes:
  #   <out_file>.png

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
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 12)
    )
  print(pp)

  ggsave(
    paste0(out_file, ".png"),
    units = "in",
    width = 4,
    height = 3,
    dpi = 600,
    device = "png"
  )
}

pred_group_intx <- function(df_tract, gam_obj, y_name, beh){
  # Predict Node-FA-Behavior interactions.
  #
  # Make predicted interaction plots for Exp and Con groups.
  #
  # Arguments:
  #   df_tract (dataframe) = input data for GAM
  #   gam_obj (gam) = GAM object returned by mgcv, specifically the
  #     function gam_GSintx
  #   y_name (str) = name of Y-axis label
  #   beh (str) = column name for covariate of interest
  #
  # Returns:
  #   named list of control (con) and experimental (exp) plots
  
  # setup for predicting con, exp smooths
  node_list <- unique(df_tract$nodeID)
  num_node <- length(node_list)
  
  df_con <- df_tract[which(df_tract$dx_group == "Con"), ]
  df_con_cov <- df_con[which(df_con$nodeID == node_list[1]), ]
  subj_con <- as.character(df_con_cov$subjectID)
  num_con <- length(subj_con)
  seq_con_cov <- seq(
    min(df_con_cov[, beh]), 
    max(df_con_cov[, beh]), 
    length = num_con
  )
  
  df_exp <- df_tract[which(df_tract$dx_group == "Exp"), ]
  df_exp_cov <- df_exp[which(df_exp$nodeID == node_list[1]), ]
  subj_exp <- as.character(df_exp_cov$subjectID)
  num_exp <- length(subj_exp)
  seq_exp_cov <- seq(
    min(df_exp_cov[, beh]), 
    max(df_exp_cov[, beh]), 
    length = num_exp
  )
  
  # predict node-fa-beh intx, Con group - use term names from GAM obj
  df_pred_con_intx <- data.frame(
    subjectID = rep(subj_con[1], each = num_node, num_con),
    sex = df_con$sex,
    h_group = df_con$dx_group,
    nodeID = df_con$nodeID,
    h_var = rep(seq_con_cov, each = num_node)
  )
  pred_con_fit <- predict(gam_obj, df_pred_con_intx)
  
  # limit extrapolation
  ind_excl <- exclude.too.far(
    df_pred_con_intx$nodeID, df_pred_con_intx$h_var,
    df_con$nodeID, df_con[, beh],
    dist = 0.1
  )
  pred_con_fit[ind_excl] <- NA
  df_pred_con_intx <- cbind(df_pred_con_intx, fit = pred_con_fit)
  
  # predict node-fa-beh intx, Exp group, limit extrapolation
  df_pred_exp_intx <- data.frame(
    subjectID = rep(subj_exp[1], each = num_node, num_exp),
    sex = df_exp$sex,
    h_group = df_exp$dx_group,
    nodeID = df_exp$nodeID,
    h_var = rep(seq_exp_cov, each = num_node)
  )
  pred_exp_fit <- predict(tract_GSintx, df_pred_exp_intx)
  
  ind_excl <- exclude.too.far(
    df_pred_exp_intx$nodeID, df_pred_exp_intx$h_var,
    df_exp$nodeID, df_exp[, beh],
    dist = 0.1
  )
  pred_exp_fit[ind_excl] <- NA
  df_pred_exp_intx <- cbind(df_pred_exp_intx, fit = pred_exp_fit)
  
  # get min/max fit values, for consistint plot scales
  z_min <- min(c(df_pred_con_intx$fit, df_pred_exp_intx$fit), na.rm = T)
  z_max <- max(c(df_pred_con_intx$fit, df_pred_exp_intx$fit), na.rm = T)
  
  # plot
  pC <- ggplot(df_pred_con_intx, aes(x = nodeID, y = h_var, z = fit)) +
    geom_tile(aes(fill = fit)) +
    geom_contour(colour = "black") +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    scale_fill_viridis(
      option = "D", 
      name = "Est. FA Fit",
      limits = c(z_min, z_max)
    ) +
    labs(y = y_name, x = "Tract Node") +
    ggtitle("Tract Node-FA-Covariate Interaction, Con") +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 12)
    )
  print(pC)
  
  pE <- ggplot(df_pred_exp_intx, aes(x = nodeID, y = h_var, z = fit)) +
    geom_tile(aes(fill = fit)) +
    geom_contour(colour = "black") +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    scale_fill_viridis(
      option = "D", 
      name = "Est. FA Fit",
      limits = c(z_min, z_max)
    ) +
    labs(y = y_name, x = "Tract Node") +
    ggtitle(paste("Tract Node-FA-Covariate Interaction, Exp")) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 12)
    )
  print(pE)
  return(list("con" = pC, "exp" = pE))
}


pred_group_intx_diff <- function(df_tract, gam_obj, y_name, beh){
  # Title.
  #
  # Desc.
  
  # # For testing
  # gam_obj <- tract_GSintxOF
  # y_name <- "Neutral LGI"
  # id_node <- 37
  
  # set up for predicting node-fa-beh group difference intx
  node_list <- unique(df_tract$nodeID)
  num_node <- length(node_list)
  
  df_exp <- df_tract[which(df_tract$dx_group == "Exp"), ]
  df_exp_cov <- df_exp[which(df_exp$nodeID == node_list[1]), ]
  subj_exp <- as.character(df_exp_cov$subjectID)
  num_exp <- length(subj_exp)
  seq_exp_cov <- seq(
    min(df_exp_cov[, beh]), 
    max(df_exp_cov[, beh]), 
    length = num_exp
  )
  
  df_pred_diff <- data.frame(
    subjectID = rep(subj_exp[1], each = num_node, num_exp),
    sex = df_exp$sex,
    h_group = df_exp$dx_group,
    h_groupOF = df_exp$dx_groupOF,
    nodeID = df_exp$nodeID,
    h_var = rep(seq_exp_cov, each = num_node)
  )
  pred_intx_diff <- as.data.frame(predict.gam(
    gam_obj, df_pred_diff, type = "terms"
  ))
  
  colnames(pred_intx_diff) <- c(
    "sex", "subjectID", "nodeID", "behRef", "behExp", "intxRef", "fit"
  )
  
  df_pred_diff <- cbind(df_pred_diff, fit = pred_intx_diff$fit)
  p <- ggplot(df_pred_diff, aes(x = nodeID, y = h_var, z = fit)) +
    geom_tile(aes(fill = fit)) +
    geom_contour(colour = "black") +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    scale_fill_viridis(option = "D", name = "Est. FA Fit") +
    labs(y = y_name, x = "Tract Node") +
    ggtitle("Tract Node-FA-Covariate Interaction, Diff") +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 12)
    )
  print(p)
  
  return(list("diff" = p))
}


pred_group_behs <- function(df_tract, id_node, gam_obj, x_name, beh, tract_name){
  # Predict FA-Behavior interactions.
  #
  # First, make prediction smooths of covariate and FA for each group at
  # a specific node. Then find the difference between group covariate smooths.
  #
  # Arguments:
  #   df_tract (dataframe) = input data for GAM
  #   id_node (int) = Node ID to hold constant
  #   gam_obj (gam) = GAM object returned by mgcv, specifically the
  #     function gam_GSintx
  #   y_name (str) = name of Y-axis label
  #   beh (str) = column name for covariate of interest
  #   tract_name (str) = long-formatted tract name for plot title
  #
  # Returns:
  #   named list of control (con) and experimental (exp) plots
  
  # # For testing
  # gam_obj <- tract_GSintx
  # x_name <- "Neutral LGI"
  # id_node <- 37
  # tract_name <- "L. Uncinate"
  
  # set up for predicting
  num_node <- length(unique(df_tract$nodeID))
  
  df_con <- df_tract[which(df_tract$dx_group == "Con"), ]
  df_con_cov <- df_con[which(df_con$nodeID == id_node), ]
  subj_con <- as.character(df_con_cov$subjectID)
  num_con <- length(subj_con)
  seq_con_cov <- seq(
    min(df_con_cov[, beh]), 
    max(df_con_cov[, beh]), 
    length = num_con
  )
  
  df_exp <- df_tract[which(df_tract$dx_group == "Exp"), ]
  df_exp_cov <- df_exp[which(df_exp$nodeID == id_node), ]
  subj_exp <- as.character(df_exp_cov$subjectID)
  num_exp <- length(subj_exp)
  seq_exp_cov <- seq(
    min(df_exp_cov[, beh]), 
    max(df_exp_cov[, beh]), 
    length = num_exp
  )
  
  # predict beh smooths, Con
  df_pred_con_cov <- data.frame(
    subjectID = rep(subj_con[1], num_con),
    sex = rep("F", num_con),
    h_group = df_con_cov$dx_group,
    nodeID = df_con_cov$nodeID,
    h_var = seq_con_cov
  )
  pred_con_cov <- predict(gam_obj, df_pred_con_cov, se.fit = T)
  pred_con_cov <- transform(
    pred_con_cov, ub = fit + (2*se.fit), lb = fit - (2*se.fit)
  )
  df_pred_con_cov <- cbind(df_pred_con_cov, pred_con_cov)
  
  # predict beh smooths, Exp
  df_pred_exp_cov <- data.frame(
    subjectID = rep(subj_exp[1], num_exp),
    sex = rep("F", num_exp),
    h_group = df_exp_cov$dx_group,
    nodeID = df_exp_cov$nodeID,
    h_var = seq_exp_cov
  )
  pred_exp_cov <- predict(gam_obj, df_pred_exp_cov, se.fit = T)
  pred_exp_cov <- transform(
    pred_exp_cov, ub = fit + (2*se.fit), lb = fit - (2*se.fit)
  )
  df_pred_exp_cov <- cbind(df_pred_exp_cov, pred_exp_cov)

  # Plot
  y_min <- min(c(df_pred_con_cov$lb, df_pred_exp_cov$lb), na.rm = T)
  y_max <- max(c(df_pred_con_cov$ub, df_pred_exp_cov$ub), na.rm = T)
  
  pC <- ggplot(data = df_pred_con_cov, aes(x = h_var, y = fit)) +
    geom_line() +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2) +
    labs(x = x_name, y = "Est. FA Fit") +
    ggtitle(paste(tract_name, "Node", id_node, "FA-Covariate Smooth, Con")) +
    theme(
      text = element_text(family = "Times New Roman"),
      legend.position = "none",
      plot.title = element_text(size = 12)
    ) +
    coord_cartesian(ylim = c(y_min, y_max))
  print(pC)
  
  pE <- ggplot(data = df_pred_exp_cov, aes(x = h_var, y = fit)) +
    geom_line() +
    ylim(y_min, y_max) +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2) +
    labs(x = x_name, y = "Est. FA Fit") +
    ggtitle(paste(tract_name, "Node", id_node, "FA-Covariate Smooth, Exp")) +
    theme(
      text = element_text(family = "Times New Roman"),
      legend.position = "none",
      plot.title = element_text(size = 12)
    ) +
    coord_cartesian(ylim = c(y_min, y_max))
  print(pE)
  
  # plot group difference smooth of covariate-fa intx as node X
  p_data <- plot_diff(
    gam_obj, 
    view = "h_var", 
    comp = list(h_group = c("Exp", "Con")),
    cond = list(nodeID = id_node),
    rm.ranef = T,
    plot = F
  )
  
  # determine regions that differ from zero
  p_data$lb <- p_data$est - p_data$CI
  p_data$ub <- p_data$est + p_data$CI
  sig_rows <- which(
    (p_data$est < 0 & p_data$ub < 0) |
      (p_data$est > 0 & p_data$lb > 0)
  )
  sig_nodes <- p_data[sig_rows, ]$h_var
  
  # find start, end points of sig regions
  vec_start <- sig_nodes[1]
  vec_end <- vector()
  y_min <- min(p_data$lb)
  num_rows <- length(sig_rows)
  c <- 1
  while (c < num_rows) {
    cc <- c + 1
    if (sig_rows[cc] > sig_rows[c] + 1) {
      vec_end <- append(vec_end, sig_nodes[c])
      vec_start <- append(vec_start, sig_nodes[cc])
    }
    c <- cc
  }
  vec_end <- append(vec_end, sig_nodes[num_rows])
  
  # make df for drawing rectangles
  d_rect <- data.frame(
    x_start = vec_start,
    x_end = vec_end,
    y_start = rep(y_min, length(vec_start)),
    y_end = rep(0, length(vec_start))
  )
  d_rect$x_start <- d_rect$x_start
  d_rect$x_end <- d_rect$x_end
  
  # draw smooth, shade diff regions
  pD <- ggplot(data = p_data, aes(x = h_var, y = est)) +
    geom_hline(yintercept = 0) +
    geom_line() +
    geom_ribbon(aes(ymin = .data$lb, ymax = .data$ub), alpha = 0.2) +
    annotate(
      "rect",
      xmin = c(d_rect$x_start),
      xmax = c(d_rect$x_end),
      ymin = c(d_rect$y_start),
      ymax = c(d_rect$y_end),
      alpha = 0.2,
      fill = "red"
    ) +
    labs(x = x_name, y = "Est. FA Fit") +
    ggtitle(paste(tract_name, "Node ", id_node,"FA-Covariate Smooth, Diff")) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 12)
    )
  print(pD)
  
  
  return(list("con" = pC, "exp" = pE, "diff" = pD))
}


