pred_group_covs <- function(df_tract, id_node, gam_obj, beh){
  # Predict FA-Behavior interactions.
  #
  # First, make prediction smooths of covariate and FA for each group at
  # a specific node. Then find the difference between group covariate smooths.
  # Excludes Sex term.
  #
  # Arguments:
  #   df_tract (dataframe) = input data for GAM
  #   id_node (int) = Node ID to hold constant
  #   gam_obj (gam) = GAM object returned by mgcv, specifically the
  #     function gam_GSintx
  #   beh (str) = column name for covariate of interest
  #
  # Returns:
  #   Named list of control (con) and experimental (exp) plots
  
  # set up for predicting
  num_node <- length(unique(df_tract$nodeID))
  
  df_con <- df_tract[which(df_tract$dx_group == "Con"), ]
  df_con <- na.omit(df_con)
  df_con_cov <- df_con[which(df_con$nodeID == id_node), ]
  subj_con <- as.character(df_con_cov$subjectID)
  num_con <- length(subj_con)
  seq_con_cov <- seq(
    min(df_con_cov[, beh]), 
    max(df_con_cov[, beh]), 
    length = num_con
  )
  
  df_exp <- df_tract[which(df_tract$dx_group == "Exp"), ]
  df_exp <- na.omit(df_exp)
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
    sex = df_con_cov$sex,
    h_group = df_con_cov$dx_group,
    nodeID = df_con_cov$nodeID,
    h_var = seq_con_cov
  )
  pred_con_cov <- predict(gam_obj, df_pred_con_cov, se.fit = T, exclude = "sex")
  pred_con_cov <- transform(
    pred_con_cov, ub = fit + (2*se.fit), lb = fit - (2*se.fit)
  )
  df_pred_con_cov <- cbind(df_pred_con_cov, pred_con_cov)
  
  # predict beh smooths, Exp
  df_pred_exp_cov <- data.frame(
    subjectID = rep(subj_exp[1], num_exp),
    sex = df_exp_cov$sex,
    h_group = df_exp_cov$dx_group,
    nodeID = df_exp_cov$nodeID,
    h_var = seq_exp_cov
  )
  pred_exp_cov <- predict(gam_obj, df_pred_exp_cov, se.fit = T, exclude = "sex")
  pred_exp_cov <- transform(
    pred_exp_cov, ub = fit + (2*se.fit), lb = fit - (2*se.fit)
  )
  df_pred_exp_cov <- cbind(df_pred_exp_cov, pred_exp_cov)
  
  # Plot
  x_min <- min(c(df_pred_con_cov$lb, df_pred_exp_cov$lb), na.rm = T)
  x_max <- max(c(df_pred_con_cov$ub, df_pred_exp_cov$ub), na.rm = T)
  
  # ggtitle(paste("Node", id_node, "FA-Behavior Smooth, Con"))
  pC <- ggplot(data = df_pred_con_cov, aes(x = h_var, y = fit)) +
    geom_line() +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2) +
    theme(
      text = element_text(family = "Times New Roman"),
      legend.position = "none",
      plot.title = element_text(size = 12),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    ) +
    coord_flip(ylim = c(x_min, x_max)) 
  # print(pC)
  
  # ggtitle(paste("Node", id_node, "FA-Behavior Smooth, Exp"))
  pE <- ggplot(data = df_pred_exp_cov, aes(x = h_var, y = fit)) +
    geom_line() +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2) +
    theme(
      text = element_text(family = "Times New Roman"),
      legend.position = "none",
      plot.title = element_text(size = 12),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    ) +
    coord_flip(ylim = c(x_min, x_max)) 
  # print(pE)
  
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
  
  # deal with whether sig nodes exist
  sig_nodes <- p_data[sig_rows, ]$h_var
  if(length(sig_nodes) > 0){
    
    # find start, end points of sig regions
    vec_start <- sig_nodes[1]
    vec_end <- vector()
    y_min <- min(p_data$lb)
    y_max <- max(p_data$ub)
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
      y_end = rep(y_max, length(vec_start))
    )
    d_rect$x_start <- d_rect$x_start
    d_rect$x_end <- d_rect$x_end
    
    # draw smooth, shade diff regions
    # ggtitle(paste("Node ", id_node,"FA-Behavior Smooth, Diff"))
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
      theme(
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()
      ) +
      coord_flip()
    # print(pD)
  } else {
    pD <- ggplot(data = p_data, aes(x = h_var, y = est)) +
      geom_hline(yintercept = 0) +
      geom_line() +
      geom_ribbon(aes(ymin = .data$lb, ymax = .data$ub), alpha = 0.2) +
      theme(
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()
      ) +
      coord_flip()
  }
  
  return(list("con" = pC, "exp" = pE, "diff" = pD))
}

pred_group_intx <- function(df_tract, gam_obj, beh){
  # Predict Node-FA-Behavior interactions.
  #
  # Make predicted interaction plots for Exp and Con groups. Excludes
  # Sex term.
  #
  # Arguments:
  #   df_tract (dataframe) = input data for GAM
  #   gam_obj (gam) = GAM object returned by mgcv, specifically the
  #     function gam_GSintx
  #   beh (str) = column name for covariate of interest
  #
  # Returns:
  #   named list of control (con) and experimental (exp) plots
  
  # setup for predicting con, exp smooths
  node_list <- unique(df_tract$nodeID)
  num_node <- length(node_list)
  
  df_con <- df_tract[which(df_tract$dx_group == "Con"), ]
  df_con <- na.omit(df_con)
  df_con_cov <- df_con[which(df_con$nodeID == node_list[1]), ]
  subj_con <- as.character(df_con_cov$subjectID)
  num_con <- length(subj_con)
  seq_con_cov <- seq(
    min(df_con_cov[, beh]), 
    max(df_con_cov[, beh]), 
    length = num_con
  )
  
  df_exp <- df_tract[which(df_tract$dx_group == "Exp"), ]
  df_exp <- na.omit(df_exp)
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
  pred_con_fit <- predict(gam_obj, df_pred_con_intx, exclude = "sex")
  
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
  pred_exp_fit <- predict(tract_GSintx, df_pred_exp_intx, exclude = "sex")
  
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
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 12),
      legend.text=element_text(size = 10),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    )
  # print(pC)
  
  pE <- ggplot(df_pred_exp_intx, aes(x = nodeID, y = h_var, z = fit)) +
    geom_tile(aes(fill = fit)) +
    geom_contour(colour = "black") +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    scale_fill_viridis(
      option = "D", 
      name = "Est. FA Fit",
      limits = c(z_min, z_max)
    ) +
    labs(x = "Tract Node") +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 12),
      legend.text=element_text(size = 10),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    )
  # print(pE)
  return(list("con" = pC, "exp" = pE))
}

pred_group_intx_diff <- function(df_tract, gam_obj,  beh, var_range){
  # Draw difference interaction smooth for node-fa-beh.
  #
  # Predict and plot how the Exp group differs in their interaction term
  # of node-fa-beh from the reference group. Excludes Sex term.
  #
  # Arguments:
  #   df_tract (dataframe) = input data for GAM
  #   gam_obj (gam) = GAM object returned by mgcv, specifically the
  #     function gam_GSintx
  #   beh (str) = column name for covariate of interest
  #   var_range (list) = named list containing $min and $max covariate
  #     values derived from plot_diff, for aligning covariate axes in
  #     2D, 3D difference smooths
  #
  # Returns:
  #   named list of difference (diff) plot
  
  # set up for predicting node-fa-beh group difference intx
  node_list <- unique(df_tract$nodeID)
  num_node <- length(node_list)
  
  df_exp <- df_tract[which(df_tract$dx_group == "Exp"), ]
  df_exp <- na.omit(df_exp)
  df_exp_cov <- df_exp[which(df_exp$nodeID == node_list[1]), ]
  subj_exp <- as.character(df_exp_cov$subjectID)
  num_exp <- length(subj_exp)
  seq_exp_cov <- seq(var_range$min, var_range$max, length = num_exp)
  
  df_pred_diff <- data.frame(
    subjectID = rep(subj_exp[1], each = num_node, num_exp),
    sex = df_exp$sex,
    h_group = df_exp$dx_group,
    h_groupOF = df_exp$dx_groupOF,
    nodeID = df_exp$nodeID,
    h_var = rep(seq_exp_cov, each = num_node)
  )
  pred_intx_diff <- as.data.frame(predict.gam(
    gam_obj, df_pred_diff, type = "terms", exclude = "sex"
  ))
  
  colnames(pred_intx_diff) <- c(
    "subjectID", "nodeID", "behRef", "behExp", "intxRef", "fit"
  )
  
  df_pred_diff <- cbind(df_pred_diff, fit = pred_intx_diff$fit)
  p <- ggplot(df_pred_diff, aes(x = nodeID, y = h_var, z = fit)) +
    geom_tile(aes(fill = fit)) +
    geom_contour(colour = "black") +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    scale_fill_viridis(option = "D", name = "Est. FA Fit") +
    labs(x = "Tract Node") +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 12),
      legend.text=element_text(size = 10),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    )
  # print(p)
  
  return(list("diff" = p))
}

pred_group_sex_covs <- function(df_tract, id_node, gam_obj, beh){
  # Predict FA-Behavior interactions, for each Sex.
  #
  # First, make prediction smooths of covariate and FA for each group at
  # a specific node. Then find the difference between group covariate smooths.
  #
  # Arguments:
  #   df_tract (dataframe) = input data for GAM
  #   id_node (int) = Node ID to hold constant
  #   gam_obj (gam) = GAM object returned by mgcv, specifically the
  #     function gam_GSintx
  #   beh (str) = column name for covariate of interest
  #
  # Returns:
  #   Named list of control (con) and experimental (exp) plots
  
  # set up for predicting
  num_node <- length(unique(df_tract$nodeID))
  sex_list <- unique(as.character(df_tract$sex))
  
  # get range across both sexes
  df_conNode <- df_tract[which(
    df_tract$dx_group == "Con" & 
      df_tract$nodeID == id_node
  ), ]
  
  df_expNode <- df_tract[which(
    df_tract$dx_group == "Exp" & 
      df_tract$nodeID == id_node
  ), ]
  
  cov_min = min(c(df_expNode[, beh], df_conNode[, beh]), na.rm = T)
  cov_max = max(c(df_expNode[, beh], df_conNode[, beh]), na.rm = T)
  
  h_pred_beh <- function(df_tract, sex, beh){
    # Title.
    #
    # Desc.
    
    df_con <- 
      df_tract[which(df_tract$dx_group == "Con" & df_tract$sex == sex), ]
    df_con_cov <- df_con[which(df_con$nodeID == id_node), ]
    subj_con <- as.character(df_con_cov$subjectID)
    num_con <- length(subj_con)
    seq_con_cov <- seq(
      min(df_con_cov[, beh]), 
      max(df_con_cov[, beh]), 
      length = num_con
    )
    
    df_exp <- 
      df_tract[which(df_tract$dx_group == "Exp" & df_tract$sex == sex), ]
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
      sex = rep(sex, num_con),
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
      sex = rep(sex, num_exp),
      h_group = df_exp_cov$dx_group,
      nodeID = df_exp_cov$nodeID,
      h_var = seq_exp_cov
    )
    pred_exp_cov <- predict(gam_obj, df_pred_exp_cov, se.fit = T)
    pred_exp_cov <- transform(
      pred_exp_cov, ub = fit + (2*se.fit), lb = fit - (2*se.fit)
    )
    df_pred_exp_cov <- cbind(df_pred_exp_cov, pred_exp_cov)
    
    return(list("df_con" = df_pred_con_cov, "df_exp" = df_pred_exp_cov))
  }
  
  df_all <- as.data.frame(matrix(NA, nrow = 1, ncol = 9))
  colnames(df_all) <- c(
    "subjectID", "sex", "h_group", "nodeID", 
    "h_var", "fit", "se.fit", "lb", "ub"
  )
  
  for(sex in sex_list){
    df_list <- h_pred_beh(df_tract, sex, beh)
    df_all <- rbind(df_all, df_list$df_con, df_list$df_exp)
  }
  
  df_con <- df_all[which(df_all$h_group == "Con"), ]
  pC <- ggplot(data = df_con, aes(x = h_var, y = fit, colour = sex)) +
    geom_line() +  
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 12),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(), 
      legend.position = "none"
    ) +
    coord_flip(xlim = c(cov_min, cov_max)) 
  print(pC)
  
  # ggtitle(paste("Node", id_node, "FA-Behavior Smooth, Exp"))
  df_exp <- df_all[which(df_all$h_group == "Exp"), ]
  pE <- ggplot(data = df_exp, aes(x = h_var, y = fit, colour = sex)) +
    geom_line() +  
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 12),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(), 
      legend.position = "none"
    ) +
    coord_flip(xlim = c(cov_min, cov_max)) 
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
  if(length(sig_nodes) > 0){
    
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
    # ggtitle(paste("Node ", id_node,"FA-Behavior Smooth, Diff"))
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
      theme(
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()
      ) +
      coord_flip()
    # print(pD)
  } else {
    pD <- ggplot(data = p_data, aes(x = h_var, y = est)) +
      geom_hline(yintercept = 0) +
      geom_line() +
      geom_ribbon(aes(ymin = .data$lb, ymax = .data$ub), alpha = 0.2) +
      theme(
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()
      ) +
      coord_flip()
  }
  
  return(list("con" = pC, "exp" = pE, "diff" = pD))
}

pred_group_sex_intx <- function(df_tract, gam_obj, beh){
  # Predict Node-FA-Behavior interactions for each sex.
  #
  # Make predicted interaction plots for Exp and Con groups.
  #
  # Arguments:
  #   df_tract (dataframe) = input data for GAM
  #   gam_obj (gam) = GAM object returned by mgcv, specifically the
  #     function gam_GSintx
  #   beh (str) = column name for covariate of interest
  #
  # Returns:
  #   named list of control (con) and experimental (exp) plots
  
  # setup for predicting con, exp smooths
  node_list <- unique(df_tract$nodeID)
  num_node <- length(node_list)
  sex_list <- unique(as.character(df_tract$sex))
  
  # get range across both sexes
  df_conNode <- df_tract[which(
    df_tract$dx_group == "Con" & 
      df_tract$nodeID == node_list[1]
  ), ]
  
  df_expNode <- df_tract[which(
    df_tract$dx_group == "Exp" & 
      df_tract$nodeID == node_list[1]
  ), ]
  
  cov_min = min(c(df_expNode[, beh], df_conNode[, beh]), na.rm = T)
  cov_max = max(c(df_expNode[, beh], df_conNode[, beh]), na.rm = T)
  
  h_pred_intx <- function(df_tract, sex, beh, node_list, num_node){
    # Title.
    #
    # Desc.
    
    df_con <- 
      df_tract[which(df_tract$dx_group == "Con" & df_tract$sex == sex), ]
    df_con_cov <- df_con[which(df_con$nodeID == node_list[1]), ]
    subj_con <- as.character(df_con_cov$subjectID)
    num_con <- length(subj_con)
    seq_con_cov <- seq(
      min(df_con_cov[, beh]), 
      max(df_con_cov[, beh]), 
      length = num_con
    )
    
    df_exp <- 
      df_tract[which(df_tract$dx_group == "Exp" & df_tract$sex == sex), ]
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
    
    # # limit extrapolation
    # ind_excl <- exclude.too.far(
    #   df_pred_con_intx$nodeID, df_pred_con_intx$h_var,
    #   df_con$nodeID, df_con[, beh],
    #   dist = 0.1
    # )
    # pred_con_fit[ind_excl] <- NA
    df_pred_con_intx <- cbind(df_pred_con_intx, fit = pred_con_fit)
    
    # predict beh smooths, Exp
    df_pred_exp_intx <- data.frame(
      subjectID = rep(subj_exp[1], each = num_node, num_exp),
      sex = df_exp$sex,
      h_group = df_exp$dx_group,
      nodeID = df_exp$nodeID,
      h_var = rep(seq_exp_cov, each = num_node)
    )
    pred_exp_fit <- predict(gam_obj, df_pred_exp_intx)
    
    # # limit extrapolation
    # ind_excl <- exclude.too.far(
    #   df_pred_exp_intx$nodeID, df_pred_exp_intx$h_var,
    #   df_exp$nodeID, df_exp[, beh],
    #   dist = 0.1
    # )
    # pred_exp_fit[ind_excl] <- NA
    df_pred_exp_intx <- cbind(df_pred_exp_intx, fit = pred_exp_fit)
    
    return(list("df_con" = df_pred_con_intx, "df_exp" = df_pred_exp_intx))
  }
  
  df_all <- as.data.frame(matrix(NA, nrow = 1, ncol = 6))
  colnames(df_all) <- c(
    "subjectID", "sex", "h_group", "nodeID", "h_var", "fit"
  )
  
  for(sex in sex_list){
    df_list <- h_pred_intx(df_tract, sex, beh, node_list, num_node)
    df_all <- rbind(df_all, df_list$df_con, df_list$df_exp)
  }
  df_all <- df_all[-1, ]
  
  # df_all <- rbind(df_pred_con_intx, df_pred_exp_intx)
  
  # plot con male, female
  df_cM <- df_all[which(df_all$h_group == "Con" & df_all$sex == "M"), ]
  df_cM$sex <- "Male"
  df_cF <- df_all[which(df_all$h_group == "Con" & df_all$sex == "F"), ]
  df_cF$sex <- "Female"
  
  pCM <- ggplot(df_cM, aes(x = nodeID, y = h_var, z = fit)) +
    geom_tile(aes(fill = fit)) +
    geom_contour(colour = "black") +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    scale_fill_viridis(
      option = "G", 
      name = "Est. FA Fit"
    ) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 12),
      legend.text=element_text(size = 10),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    )
  
  pCF <- ggplot(df_cF, aes(x = nodeID, y = h_var, z = fit)) +
    geom_tile(aes(fill = fit)) +
    geom_contour(colour = "black") +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    scale_fill_viridis(
      option = "F", 
      name = "Est. FA Fit"
    ) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 12),
      legend.text=element_text(size = 10),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    )
  # print(pC)
  
  # plot exp male, female
  df_eM <- df_all[which(df_all$h_group == "Exp" & df_all$sex == "M"), ]
  df_eM$sex <- "Male"
  df_eF <- df_all[which(df_all$h_group == "Exp" & df_all$sex == "F"), ]
  df_eF$sex <- "Female"
  
  pEM <- ggplot(df_eM, aes(x = nodeID, y = h_var, z = fit)) +
    geom_tile(aes(fill = fit)) +
    geom_contour(colour = "black") +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    scale_fill_viridis(
      option = "G", 
      name = "Est. FA Fit"
    ) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 12),
      legend.text=element_text(size = 10),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    )
  
  pEF <- ggplot(df_eF, aes(x = nodeID, y = h_var, z = fit)) +
    geom_tile(aes(fill = fit)) +
    geom_contour(colour = "black") +
    scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
    scale_fill_viridis(
      option = "F", 
      name = "Est. FA Fit"
    ) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 12),
      legend.text=element_text(size = 10),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    )
  # print(pE)
  return(list("conM" = pCM, "conF" = pCF, "expM" = pEM, "expF" = pEF))
}
