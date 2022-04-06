library("tools")
library("dplyr")
library("ggplot2")
library("ez")


# functions ----
adjust_outliers <- function(df, col_name){
  # Replace outliers with max/min values.
  #
  # Arguments:
  #   
  # Returns:
  #
  
  # find min/max
  h_iqr <- IQR(df[, col_name], na.rm = TRUE)
  h_quant <- quantile(df[, col_name], na.rm = TRUE, names = FALSE)
  h_min <- h_quant[2] - (1.5 * h_iqr)
  h_max <- h_quant[4] + (1.5 * h_iqr)
  
  # detect, replace outliers
  ind_out <- which(df[, col_name] < h_min | df[, col_name] > h_max)
  for(ind in ind_out){
    if(df[ind, col_name] > h_max){
      df[ind, col_name] <- h_max
    }else if(df[ind, col_name] < h_min){
      df[ind, col_name] <- h_min
    }
  }
  return(df)
}

linear_model <- function(df, roi, col_name, data_dir){
  # Title.
  #
  # Arguments:
  #   
  # Returns:
  #
  h_fit <- lm(df[, col_name] ~ pscared * dx_group, data = df)
  h_aov <- anova(h_fit)
  out_fit <- paste0(
    data_dir, "/ppi_stats/LM_", roi, "_", col_name, "-pScared.txt"
  )
  writeLines(capture.output(summary(h_fit)), out_fit)
}

switch_roi <- function(h_str) {
  # Title.
  #
  # Arguments:
  #
  # Returns:
  #
  h_roi <- switch(
    h_str,
    "NSlacc" = "ACC",
    "NSldmpfc" = "dmPFC",
    "NSlsfs" = "SFS",
  )
  return(h_roi)
}

switch_beh <- function(h_str) {
  # Title.
  #
  # Arguments:
  #
  # Returns:
  #
  h_beh <- switch(
    h_str,
    "SnegLF" = "Negative",
    "SneuLF" = "Neutral",
    "SPnegLF" = "Prec. Negative",
    "SPneuLF" = "Prec. Neutral",
  )
  return(h_beh)
}

linear_plot <- function(df, roi, col_name, data_dir){
  # Title.
  #
  # Arguments:
  #
  # Returns:
  #
  h_roi <- switch_roi(roi)
  y_lab <- paste(switch_beh(col_name), "Lure FA PPI Term")
  ggplot(data = df, aes(x = pscared, y = df[, col_name], color = dx_group)) +
    geom_point() +
    geom_ribbon(
      stat = "smooth", 
      method = "lm", 
      se = T, 
      alpha=0.1,
      aes(color = NULL, group = factor(dx_group))
    ) +
    geom_line(stat = "smooth", method = "lm", alpha = 1) +
    ggtitle(paste0("L. Amg -- ", h_roi)) +
    ylab(y_lab) +
    xlab("Parent's SCARED") +
    labs(color="Dx Group")
  
  out_plot <- paste0(
    data_dir, "/ppi_plots/LM_", roi, "_", col_name, "-pScared.png"
  )
  ggsave(
    out_plot,
    units = "in",
    width = 6,
    height = 6,
    device = "png"
  )
}

long_format <- function(df, beh_list){
  # Title.
  #
  # Arguments:
  #
  # Returns:
  #
  subj_list <- unique(df$subj)
  num_subj <- length(subj_list)
  num_beh <- length(beh_list)
  
  df_long <- as.data.frame(
    matrix(NA, nrow=num_subj*num_beh, ncol=4)
  )
  colnames(df_long) <- c("subj", "group", "beh", "value")
  df_long$subj <- rep(subj_list, each = num_beh)
  df_long$beh <- rep(beh_list, num_subj)
  
  for(h_row in 1:dim(df)[1]){
    ind_subj <- which(df_long$subj == df[h_row,]$subj)
    df_long[ind_subj, ]$group <- as.character(df[h_row, ]$dx_group)
    for(h_beh in beh_list){
      ind_beh <- which(df_long$subj == df[h_row,]$subj & df_long$beh == h_beh)
      df_long[ind_beh, ]$value <- df[h_row, h_beh]
    }
  }
  df_long$group <- factor(df_long$group)
  return(df_long)
}

anova_model <- function(df_long, roi, beh_list, data_dir){
  # Title.
  #
  # Arguments:
  #
  # Returns:
  #
  h_fit <- ezANOVA(
    df_long, 
    dv = value, 
    wid = subj, 
    within = beh, 
    between = group
  )
  out_fit <- paste0(
    data_dir, "/ppi_stats/AN_", roi, ".txt"
  )
  writeLines(capture.output(h_fit), out_fit)
}

anova_plot <- function(df_long, roi, beh_list, data_dir){
  # Title.
  #
  # Arguments:
  #
  # Returns:
  #
  h_roi <- switch_roi(roi)
  y_lab <- "Lure FA PPI Term"
  ind_a <- which(df_long$beh == beh_list[1])
  ind_b <- which(df_long$beh == beh_list[2])
  df_long[ind_a,]$beh <- switch_beh(beh_list[1])
  df_long[ind_b,]$beh <- switch_beh(beh_list[2])
  
  ggplot(df_long, aes(x = beh, y = value, fill = group)) +
    geom_boxplot() +
    ggtitle(paste0("L. Amg -- ", h_roi)) +
    labs(fill="Dx Group") +
    ylab(y_lab) +
    scale_x_discrete(name = "Stimulus Valence")
  
  out_plot <- paste0(
    data_dir, "/ppi_plots/AN_", roi, ".png"
  )
  ggsave(
    out_plot,
    units = "in",
    width = 6,
    height = 6,
    device = "png"
  )
}


# Set Up ----
data_dir <- file_path_as_absolute(paste0(getwd(), "/../data"))


# amgL NSlacc, NSldmpfc, NSlsfs ----
roi_list <- c("NSlacc", "NSldmpfc", "NSlsfs")
beh_list <- c("SPnegLF", "SPneuLF")
sess <- "ses-S1"
task <- "task-study"

# omnibus
df_a <- read.csv(
  paste0(data_dir, "/df_", sess, "_", task, "_", "amgL-", roi_list[1], ".csv")
)
df_a$roi <- roi_list[1]
df_b <- read.csv(
  paste0(data_dir, "/df_", sess, "_", task, "_", "amgL-", roi_list[2], ".csv")
)
df_b$roi <- roi_list[2]
df_c <- read.csv(
  paste0(data_dir, "/df_", sess, "_", task, "_", "amgL-", roi_list[3], ".csv")
)
df_c$roi <- roi_list[3]

for(beh in beh_list){
  df_a <- adjust_outliers(df_a, beh)
  df_b <- adjust_outliers(df_b, beh)
  df_c <- adjust_outliers(df_c, beh)
}
df_master <- rbind(df_a, df_b, df_c)
df_master <- na.omit(df_master)
rm(df_a)
rm(df_b)
rm(df_c)

num_roi <- length(roi_list)
num_beh <- length(beh_list)
subj_list <- unique(df_master$subj)
num_subj <- length(subj_list)
col_names <- c("subj", "group", "roi", "beh", "value")
df_long <- as.data.frame(
  matrix(NA, nrow=num_subj*num_roi*num_beh, ncol=length(col_names))
)
colnames(df_long) <- col_names

df_long$subj <- rep(subj_list, each = num_roi * num_beh)
df_long$roi <- rep(roi_list, each = num_beh, times = num_subj)
df_long$beh <- rep(beh_list, times = num_roi * num_subj)
for(subj in subj_list){
  for(roi in roi_list){
    ind_wide <- which(df_master$subj == subj & df_master$roi == roi)
    for(beh in beh_list){
      ind_long <- which(
        df_long$subj == subj & df_long$roi == roi & df_long$beh == beh
      )
      df_long[ind_long,]$value <- df_master[ind_wide, beh]
      df_long[ind_long,]$group <- df_master[ind_wide, ]$dx_group
    }
  }
}
df_long$roi <- factor(df_long$roi)
df_long$group <- factor(df_long$group)

fit_omni <- ezANOVA(
  df_long, value, wid=subj, within = c(beh, roi), between = group
)
fit_omni
ggplot(df_long, aes(x = beh, y = value, fill = group)) +
  facet_wrap(~roi) +
  geom_boxplot()


# post-hoc
for(roi in roi_list){
  
  # get data
  df <- read.csv(
    paste0(data_dir, "/df_", sess, "_", task, "_", "amgL-", roi, ".csv")
  )
  df$dx_group <- factor(df$dx_group)
  df <- na.omit(df)
  
  # deal w/outliers
  for(beh in beh_list){
    df <- adjust_outliers(df, beh)
  }
  
  # beh x group
  df_long <- long_format(df, beh_list)
  anova_model(df_long, roi, beh_list, data_dir)
  anova_plot(df_long, roi, beh_list, data_dir)
  
  # linear models
  for( beh in beh_list){
    linear_model(df, roi, beh, data_dir)
    linear_plot(df, roi, beh, data_dir)
  }
}

