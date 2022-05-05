library("tools")
library("dplyr")
library("tidyr")
library("ggplot2")
library("ez")


# Set Up ----
data_dir <- file_path_as_absolute(paste0(getwd(), "/../data"))


# ROI Coefs Analysis ----
#
# Test for an interaction between ROI (L, R Amg), scene rating (neg, neu),
# and Group (con, exp) in beta-coefficient.

# use only subjects who have AFQ data
df_afq <- read.csv(paste0(data_dir, "/AFQ_dataframe.csv"))
subj_list <- unique(df_afq$subjectID)
subj_list <- paste0("sub-", subj_list)
num_subj <- length(subj_list)

# get L/R amg data
df_amgL <- read.csv(
  paste0(data_dir, "/df_ses-S1_task-study_decon-rVal_amgL.csv")
)
df_amgR <- read.csv(
  paste0(data_dir, "/df_ses-S1_task-study_decon-rVal_amgR.csv")
)

# construct long df
beh_list <- c("neg", "neu")
roi_list <- c("amgL", "amgR")
num_beh <- length(beh_list)
num_roi <- length(roi_list)

df_long <- as.data.frame(
  matrix(NA, nrow = num_subj * num_beh * num_roi, ncol = 5)
)
colnames(df_long) <- c("subj", "group", "roi", "beh", "coef")
df_long$subj <- rep(subj_list, each = num_beh * num_roi)
df_long$roi <- rep(rep(roi_list, each = num_beh), num_subj)
df_long$beh <- rep(beh_list, num_subj * num_roi)

#  mine data for e/subj
for(subj in subj_list){

  # get group
  ind_long_subj <- which(df_long$subj == subj)
  ind_amgL_subj <- which(df_amgL$subj == subj)
  if(length(ind_amgL_subj) == 0){
    next
  }
  df_long[ind_long_subj, ]$group <- df_amgL[ind_amgL_subj, ]$dx_group

  # get roi beh coefs
  for(roi in roi_list){
    h_df <- get(paste0("df_", roi))
    for(beh in beh_list){
      ind_long <- which(df_long$roi == roi & df_long$beh == beh & df_long$subj == subj)
      ind_beh <- which(h_df$subj == subj)
      df_long[ind_long, ]$coef <- h_df[ind_beh, beh]
    }
    rm(h_df)
  }
}

# clean up df
df_long <- df_long[complete.cases(df_long), ]
df_long$group <- factor(df_long$group)
df_long$roi <- factor(df_long$roi)
df_long$beh <- factor(df_long$beh)

# get final subj num
final_subj_num <- length(unique(df_long$subj))

# omnibus test
fit_omni <- ezANOVA(
  df_long, coef, wid=subj, within = c(beh, roi), between = group
)
fit_omni # ME beh, roi only

# rename vars for pretty plots
ind_neg <- which(df_long$beh == "neg")
ind_neu <- which(df_long$beh == "neu")
df_long$beh <- as.character(df_long$beh)
df_long[ind_neg, ]$beh <- "Neg"
df_long[ind_neu, ]$beh <- "Neu"

# plot data, save
new_labs <- c("** L. Amg", "** R. Amg")
names(new_labs) <- c("amgL", "amgR")

ggplot(df_long, aes(x = beh, y = coef, fill = group)) +
  facet_wrap(~roi, labeller = labeller(roi = new_labs)) +
  geom_boxplot() +
  annotate(
    "segment", x = 1, xend = 2, y = 0.25, yend = 0.25, color = "black"
  ) +
  annotate(
    "segment", x = 1, xend = 1, y = 0.25, yend = 0.23, color = "black"
  ) +
  annotate(
    "segment", x = 2, xend = 2, y = 0.25, yend = 0.23, color = "black"
  ) +
  annotate("text", x = 1.5, y = 0.27, label = "***") +
  labs(x = "Scene Rating", y = "Coefficient") +
  scale_fill_discrete(name = "Group") +
  ggtitle("Scene Valence Ratings, Amygdaloid Signal") +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(size=12)
    )
# ggsave(
#   "/Users/nmuncy/Desktop/roi_amg.png",
#   plot = last_plot(),
#   units = "in",
#   width = 4,
#   height = 3,
#   dpi = 600,
#   device = "png"
# )


# PPI Coefs Analysis ----
#
# Test for an interaction between ROI (acc, dmpfc, sfs), scene
# rating (neg, neu), and Group (con, exp) in PPI correlation term.
#
# Essentially the same as ROI analysis.

# use only subjects who have AFQ data
df_afq <- read.csv(paste0(data_dir, "/AFQ_dataframe.csv"))
subj_list <- unique(df_afq$subjectID)
subj_list <- paste0("sub-", subj_list)
num_subj <- length(subj_list)

# get lacc, ldmpfc, lsfs data
df_lacc <- read.csv(
  paste0(data_dir, "/df_ses-S1_task-study_decon-rVal_amgL-NSlacc.csv")
)
df_ldmpfc <- read.csv(
  paste0(data_dir, "/df_ses-S1_task-study_decon-rVal_amgL-NSldmpfc.csv")
)
df_lsfs <- read.csv(
  paste0(data_dir, "/df_ses-S1_task-study_decon-rVal_amgL-NSlsfs.csv")
)

# construct long df
beh_list <- c("neg", "neu")
roi_list <- c("lacc", "ldmpfc", "lsfs")
num_beh <- length(beh_list)
num_roi <- length(roi_list)

df_long <- as.data.frame(
  matrix(NA, nrow = num_subj * num_beh * num_roi, ncol = 5)
)
colnames(df_long) <- c("subj", "group", "roi", "beh", "coef")
df_long$subj <- rep(subj_list, each = num_beh * num_roi)
df_long$roi <- rep(rep(roi_list, each = num_beh), num_subj)
df_long$beh <- rep(beh_list, num_subj * num_roi)

for(subj in subj_list){

  # get group
  ind_long_subj <- which(df_long$subj == subj)
  ind_lacc_subj <- which(df_lacc$subj == subj)
  if(length(ind_lacc_subj) == 0){
    next
  }
  df_long[ind_long_subj, ]$group <- df_lacc[ind_lacc_subj, ]$dx_group

  # get roi beh coefs
  for(roi in roi_list){
    h_df <- get(paste0("df_", roi))
    for(beh in beh_list){
      ind_long <- which(
        df_long$roi == roi & df_long$beh == beh & df_long$subj == subj
      )
      ind_beh <- which(h_df$subj == subj)
      df_long[ind_long, ]$coef <- h_df[ind_beh, beh]
    }
    rm(h_df)
  }
}

# clean up df, get final count
df_long <- df_long[complete.cases(df_long), ]
df_long$group <- factor(df_long$group)
df_long$roi <- factor(df_long$roi)
df_long$beh <- factor(df_long$beh)
final_subj_num <- length(unique(df_long$subj))

# omnibus test
fit_omni <- ezANOVA(
  df_long, coef, wid=subj, within = c(beh, roi), between = group
)
fit_omni # ME of roi only

# rename vars for pretty plots
ind_lacc <- which(df_long$roi == "lacc")
ind_ldmpfc <- which(df_long$roi == "ldmpfc")
ind_lsfs <- which(df_long$roi == "lsfs")
ind_neg <- which(df_long$beh == "neg")
ind_neu <- which(df_long$beh == "neu")

df_long$roi <- as.character(df_long$roi)
df_long$beh <- as.character(df_long$beh)

df_long[ind_lacc, ]$roi <- "* L. ACC"
df_long[ind_ldmpfc, ]$roi <- "* L. dmPFC"
df_long[ind_lsfs, ]$roi <- "* L. SFS"
df_long[ind_neg, ]$beh <- "Neg"
df_long[ind_neu, ]$beh <- "Neu"

# plot data, save
ggplot(df_long, aes(x = beh, y = coef, fill = group)) +
  facet_wrap(~roi) +
  geom_boxplot() +
  labs(x = "Scene Rating", y = "Coefficient") +
  scale_fill_discrete(name = "Group") +
  ggtitle("Scene Valence Rating, L. Amg PPI") +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(size=12)
  )
ggsave(
  "/Users/nmuncy/Desktop/ppi_amg.png",
  plot = last_plot(),
  units = "in",
  width = 4,
  height = 3,
  dpi = 600,
  device = "png"
)
