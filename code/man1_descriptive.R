library("tidyr")
library("ggplot2")
library("mgcv")
library("parallel")
library("itsadug")
library("mgcViz")
library("fitdistrplus")
library("viridis")
library("devtools")
install_local(path = "./DiffGamm", force = T)
library("DiffGamm")


data_dir <- "/Users/nmuncy/Projects/emu_unc/data"
df_afq <- read.csv(paste0(data_dir, "/AFQ_dataframe.csv"))
df_subset <- df_afq[which(df_afq$tractID == "UNC_L" & df_afq$nodeID == 10), ]
df_subset <- df_node10 %>% drop_na(dx)

# Get Demographics ------
num_subj <- dim(df_subset)[1]
num_female <- length(which(df_subset$sex == "F"))
age_avg <- round(mean(df_subset$age), 2)
age_sd <- round(sd(df_subset$age), 2)

# stat_age <- aov(Age ~ as.factor(Group), data = df_subset)
# summary(stat_age)
# etaSquared(stat_age)
#
# stat_pds <- aov(PDS ~ as.factor(Group), data = df_subset)
# summary(stat_pds)
# etaSquared(stat_pds)


# Generate interaction smooth ----
#
# Simulate data to demonstrate (a) groups smooths differing,
# and (b) groups interacting differently with a continuous
# covariate.
gen_norm <- function(y_scale, n = 50, mean = 1, sd = 0.1) {
  # Generate a normal distribution.
  #
  # Simulate tract FA values by generating a Gaussian distribution
  # of values.
  #
  # Arguments:
  #   y_scale (numeric) = scaling factor to adjust y-range. Used so different
  #     groups have different "heights" if generated tracts/Gaussian values
  #   n (numeric) = number of values to generate
  #   mean (numeric) = mean of norm dist
  #   sd (numeric) = sd of norm dist
  #
  # Returns:
  #   list of values length=n

  # set up range of sequence
  h_min <- mean - (2.5 * sd)
  h_max <- mean + (2.5 * sd)
  h_seq <- seq(h_min, h_max, length = n)

  # generate norm points in sequence with jitter
  set.seed(y_scale)
  dist_out <- y_scale * dnorm(h_seq, mean, sd) + rnorm(n, 0, (0.15 * y_scale))
  return(dist_out)
}

# set number of subjects, name subjects, and generate differing 
# scaling values for each group
set.seed(12)
num_groupA <- num_groupB <- 30
subj_groupA <- 100 + seq(1:num_groupA)
subj_groupB <- 200 + seq(1:num_groupB)
scale_groupA <- runif(num_groupA, 7, 10)
scale_groupB <- runif(num_groupB, 7, 12)

# set up long dataframe
df_long <- as.data.frame(matrix(NA, nrow = 2 * num_groupA * 50, ncol = 6))
colnames(df_long) <- c("subj", "group", "node", "y_scale", "fa", "cov")
df_long$subj <- c(rep(subj_groupA, each = 50), rep(subj_groupB, each = 50))
df_long$group <- c(rep("A", num_groupA * 50), rep("B", num_groupB * 50))
df_long$group <- factor(df_long$group)
df_long$node <- rep(seq(1, 50), 2 * num_groupA)
df_long$y_scale <- c(rep(scale_groupA, each = 50), rep(scale_groupB, each = 50))

# fill dataframe with generated values
c_start <- 1
for (h_seed in c(scale_groupA, scale_groupB)) {
  c_end <- c_start + 49
  df_long[c_start:c_end, ]$fa <- gen_norm(h_seed)
  c_start <- c_end + 1
}

# plot node data
ggplot(data = df_long, aes(x = node, y = fa, color = group)) +
  geom_point()

# generate linear, exponential data for group covariates
seq_values <- 1:num_groupA
lin_cov <- (0.1 * seq_values + 2) + rnorm(length(seq_values), 0, 0.1)
exp_cov <- (((0.1 * seq_values)^2) + 2) + rnorm(length(seq_values), 0, 0.1)
plot(seq_values, lin_cov)
plot(seq_values, exp_cov)

# get group interaction according to ordered scale_y position
sorted_scaleA <- sort(scale_groupA)
sorted_scaleB <- sort(scale_groupB)
c <- 1
while (c <= length(subj_groupA)) {

  # determine subject ids
  subjA <- subj_groupA[c]
  subjB <- subj_groupB[c]

  # find index of y_scale value in sorted scale (position of subj in scale)
  ind_scaleA <- which(grepl(scale_groupA[c], sorted_scaleA) == T)
  ind_scaleB <- which(grepl(scale_groupB[c], sorted_scaleB) == T)

  # determine df_long rows
  ind_dfA <- which(df_long$subj == subjA)
  ind_dfB <- which(df_long$subj == subjB)

  # fill df_long rows with indexed group covariate - so subjs with
  # larger y_scale values also have larger covariates
  df_long[ind_dfA, ]$cov <- lin_cov[ind_scaleA]
  df_long[ind_dfB, ]$cov <- exp_cov[ind_scaleB]
  c <- c + 1
}

# plot intx of scaling and cov
df_ind1 <- df_long[which(df_long$node == 1), ]
ggplot(data = df_ind1, aes(x = y_scale, y = cov, color = group)) +
  geom_point()

# gam via GS method
descdist(df_long$fa, discrete = F)
gam_GS <- bam(fa ~
    s(subj, bs = "re") +
    s(node, bs = "cr", k = 10, m = 2) +
    s(node, group, bs = "fs", k = 10, m = 2),
  data = df_long,
  family = gaussian(),
  method = "fREML"
)
plot(gam_GS)

# use ordered factor to make group difference smooth
df_long$groupOF <- ordered(df_long$group)
gam_OF <- bam(fa ~
    s(subj, bs = "re") +
    s(node, bs = "cr", k = 10, m = 2) +
    s(node, by = groupOF, bs = "cr", k = 10, m = 2),
  data = df_long,
  family = gaussian(),
  method = "fREML"
)
plot(gam_OF)
# h_plot <- getViz(gam_OF)
# plot(sm(h_plot, 2))
# plot(sm(h_plot, 3))

# interaction of groupA - visualize linear fa-cov intx
df_groupA <- df_long[which(df_long$group == "A"), ]
gam_groupA <- bam(fa ~
    s(subj, bs = "re") +
    te(node, cov, bs = c("cr", "tp"), k = c(10, 10)),
  data = df_groupA,
  family = gaussian(),
  method = "fREML"
)
summary(gam_groupA)

# make contour plot
plot(gam_groupA)
h_plot <- getViz(gam_groupA)
plot(sm(h_plot, 2))
p <- plot(sm(h_plot, 2))
# df_A <- p$data$fit

# # make 3d plot
# vis.gam(gam_groupA,
#   view = c("node", "cov"),
#   color = "topo",
#   theta = 330
# )

# interaction of groupB - visualize exponential fa-cov intx
df_groupB <- df_long[which(df_long$group == "B"), ]
gam_groupB <- bam(fa ~
    s(subj, bs = "re") +
    te(node, cov, bs = c("cr", "tp"), k = c(10, 10)),
  data = df_groupB,
  family = gaussian(),
  method = "fREML"
)

# make contour
plot(gam_groupB)
h_plot <- getViz(gam_groupB)
plot(sm(h_plot, 2))
p <- plot(sm(h_plot, 2))
# df_B <- p$data$fit
# 
# df_diff <- df_B - df_A
# df_diff$x <- df_B$x
# df_diff$y <- df_B$y
# ggplot(data = df_diff, aes(x = x, y = y, z = z)) +
#   geom_tile(aes(fill = z)) +
#   geom_contour(colour = "black") +
#   scale_fill_viridis(option = "D", name = "Fit FA") +
#   labs(x = "Covariate Term", y = "NodeID") +
#   ggtitle("Exp group difference interaction smooth")


# full model - model how fa is a function of subject, node, cov, and group
# gam_cov <- bam(fa ~
#     s(subj, bs = "re") +
#     te(node, cov, bs = c("cr", "tp"), k = c(10, 10), m = 2) +
#     t2(
#       node, cov, group,
#       bs = c("cr", "tp", "re"),
#       k = c(10, 10, 2),
#       m = 2,
#       full = T
#     ),
#   data = df_long,
#   family = gaussian(),
#   method = "fREML"
# )
# gam.check(gam_cov, rep = 1000)
# summary(gam_cov)

gam_cov <- bam(fa ~ 
    s(subj, bs = "re") +
    s(node, bs = "cr", k = 10, m = 1) +
    s(cov, by = group, bs = "tp", k = 10, m = 2) +
    ti(
      node, cov, by = group, bs = c("cr", "tp"), k = c(10, 10), m = 2
    ),
  data = df_long,
  family = gaussian(),
  method = "fREML",
  discrete = T
)
plot(gam_cov)
summary(gam_cov)

# visualize
h_plot <- getViz(gam_cov)
plot(sm(h_plot, 2))
plot(sm(h_plot, 3))
plot(sm(h_plot, 4))
plot(sm(h_plot, 5))
plot(sm(h_plot, 6))

# extract plot data to compare
p_con <- plot(sm(h_plot, 5))
p_exp <- plot(sm(h_plot, 6))
df_con <- p_con$data$fit
df_exp <- p_exp$data$fit
df_diff <- df_exp - df_con
df_diff$x <- df_con$x
df_diff$y <- df_con$y
ggplot(data = df_diff, aes(x = x, y = y, z = z)) +
  geom_tile(aes(fill = z)) +
  geom_contour(colour = "black") +
  scale_fill_viridis(option = "D", name = "Fit FA") +
  labs(x = "Covariate Term", y = "NodeID") +
  ggtitle("Exp group difference interaction smooth")


# Use ordered factors to show how group B differs in its interaction term
# from group A
df_long$groupOF <- NA
df_long$groupOF <- factor(df_long$group, ordered = T)
# gam_covOF <- bam(fa ~
#     s(subj, bs = "re") +
#     te(node, cov, bs = c("cr", "tp"), k = c(10, 10), m = 2) +
#     t2(
#       node, cov,
#       by = groupOF,
#       bs = c("cr", "tp"),
#       k = c(10, 10),
#       m = 2,
#       full = T
#     ),
#   data = df_long,
#   family = gaussian(),
#   method = "fREML"
# )
# gam.check(gam_covOF, rep = 1000)

gam_covOF <- bam(fa ~ 
    s(subj, bs = "re") +
    s(node, bs = "cr", k = 10, m = 1) +
    s(cov, by = group, bs = "tp", k = 10, m = 2) +
    ti(node, cov, bs = c("cr", "tp"), k = c(10, 10), m = 2) +
    ti(
      node, cov, by = groupOF, bs = c("cr", "tp"), k = c(10, 10), m = 2
    ),
  data = df_long,
  family = gaussian(),
  method = "fREML",
  discrete = T
)
summary(gam_covOF)

# plot reference smooth
h_plot <- getViz(gam_covOF)
plot(sm(h_plot, 3))
plot(sm(h_plot, 6))

plot(sm(h_plot, 2)) +
  labs(colour = "Fit FA", x = "Node", y = "Covariate term") +
  ggtitle("Group A interaction smooth")

# plot diff smooth, seems inverted
plot(sm(h_plot, 3))

# difference is A-B, so invert for more intuitive plots
p <- plot(sm(h_plot, 3))
p_data <- p$data$fit
p_data$zI <- -1 * p_data$z

ggplot(data = p_data, aes(x = x, y = y, z = zI)) +
  geom_tile(aes(fill = zI)) +
  geom_contour(colour = "black") +
  scale_fill_viridis(option = "D", name = "Fit FA") +
  labs(x = "Node", y = "Covariate term") +
  ggtitle("Group B difference interaction smooth")

# use GI method to get access to plot_diff2 (to verify gam_covOF)
gam_covGI <- bam(fa ~
   s(group, bs = "re") +
   s(subj, bs = "re") +
   s(node, bs = "cr", k = 10, m = 2) +
   te(node, cov, by = group, bs = c("cr", "tp"), k = c(10, 10), m = 1),
  data = df_long,
  family = gaussian(),
  method = "fREML"
)
plot(gam_covGI)
h_plot <- getViz(gam_covGI)
plot(sm(h_plot, 3))
plot(sm(h_plot, 4))
plot(sm(h_plot, 5))

par(cex = 1, mar = c(5.1, 4.1, 4.1, 3.5))
plot_diff2(
  gam_covGI,
  view=c("node", "cov"),
  comp=list(group=c("B", "A")),
  rm.ranef=T
)



# Use real FA data, simulate interaction ----

# get tract data, clean up
df_tract <- df_afq[which(df_afq$tractID == "UNC_L"), ]
df_tract <- df_tract %>% drop_na(dx)
ind_pat <- which(df_tract$dx_group == "Pat")
df_tract[ind_pat, ]$dx_group <- "Exp"
df_tract$dx_group <- factor(df_tract$dx_group)
df_tract$subjectID <- factor(df_tract$subjectID)
df_tract$dx_groupOF <- factor(df_tract$dx_group, ordered = T)

ind_keep <- which(
  df_tract$nodeID >= 10 & df_tract$nodeID <= 89
)
df_tract <- df_tract[ind_keep, ]

num_groupC <- length(which(df_tract$nodeID == 50 & df_tract$dx_group == "Con"))
num_groupE <- length(which(df_tract$nodeID == 50 & df_tract$dx_group == "Exp"))


# Make group x cov interaction ----
# generate linear covariate values for control, 
# exponential for experimental group
set.seed(12)
seq_C <- 1:num_groupC
seq_E <- 1:num_groupE
lin_cov <- (0.1 * seq_C + 2) + rnorm(length(seq_C), 0, 0.1)
exp_cov <- (((0.1 * seq_E)^2) + 2) + rnorm(length(seq_E), 0, 0.1)

plot(seq_C, lin_cov)
plot(seq_E, exp_cov)

plot(seq_C, lin_cov)
points(seq_E, exp_cov)

# groups differ at nodes 30-45, use ordered node 37 to assign simulated
# covariate values -- so participants with larger 37 FA value have larger
# cov value
fa_37C <- sort(df_tract[which(
  df_tract$dx_group == "Con" & df_tract$nodeID == 37
), ]$dti_fa)

fa_37E <- sort(df_tract[which(
  df_tract$dx_group == "Exp" & df_tract$nodeID == 37
), ]$dti_fa)

# assign simulated covariates
df_tract$cov <- NA
subj_C <- unique(df_tract[which(df_tract$dx_group == "Con"), ]$subjectID)
subj_E <- unique(df_tract[which(df_tract$dx_group == "Exp"), ]$subjectID)
for(subj in subj_C){
  val_37 <- df_tract[which(
    df_tract$subjectID == subj & df_tract$nodeID == 37
  ), ]$dti_fa
  ord_37 <- which(grepl(val_37, fa_37C) == T)
  ind_tract <- which(df_tract$subjectID == subj)
  df_tract[ind_tract, ]$cov <- lin_cov[ord_37]
}

for(subj in subj_E){
  val_37 <- df_tract[which(
    df_tract$subjectID == subj & df_tract$nodeID == 37
  ), ]$dti_fa
  ord_37 <- which(grepl(val_37, fa_37E) == T)
  ind_tract <- which(df_tract$subjectID == subj)
  df_tract[ind_tract, ]$cov <- exp_cov[ord_37]
}

# plot intx of group, cov, and node37 FA
df_ind1 <- df_tract[which(df_tract$nodeID == 37), ]
ggplot(data = df_ind1, aes(x = dti_fa, y = cov, color = dx_group)) +
  geom_point()


# Model GS, individual group-cov intxs ----
# gam via GS method
descdist(df_tract$dti_fa, discrete = F)
# gam_GS <- gam_GS_model(df_tract, "gamma", "dx_group")
gam_GS <- readRDS(file = "/Users/nmuncy/Projects/emu_unc/stats/Model_UNC_L_dxGS.Rda")
# gam.check(gam_GS, rep = 1000)
plot_GS <- getViz(gam_GS)
plot(sm(plot_GS, 2))
plot(sm(plot_GS, 3))

# gam via GS method, ordered factor
gam_GSOF <- gam_GSOF_model(df_tract, "gamma", "dx_groupOF")
plot_GSOF <- getViz(gam_GSOF)
plot(sm(plot_GSOF, 1))
plot(sm(plot_GSOF, 2))
plot(sm(plot_GSOF, 3))


# interaction of groupC - visualize linear fa-cov intx
df_groupA <- df_tract[which(df_tract$dx_group == "Con"), ]
gam_groupA <- bam(dti_fa ~
    s(subjectID, bs = "re") +
    te(nodeID, cov, bs = c("cr", "tp"), k = c(50, 10)),
  data = df_groupA,
  family = Gamma(link = "logit"),
  method = "fREML"
)
summary(gam_groupA)

# make contour plot
plot_groupA <- getViz(gam_groupA)
plot(sm(plot_groupA, 2))

# interaction of groupE - visualize exponential fa-cov intx
df_groupB <- df_tract[which(df_tract$dx_group == "Exp"), ]
gam_groupB <- bam(dti_fa ~
    s(subjectID, bs = "re") +
    te(nodeID, cov, bs = c("cr", "tp"), k = c(50, 10)),
  data = df_groupB,
  family = Gamma(link = "logit"),
  method = "fREML"
)
summary(gam_groupB)

# make contour plot
plot_groupB <- getViz(gam_groupB)
plot(sm(plot_groupB, 2))


# Model group-cov interaction ----
# full interaction model
# gam_cov <- gam_intx_model(df_tract, "gamma", "dx_group", "cov")
# saveRDS(gam_cov, file = "/Users/nmuncy/Desktop/gam_intx.Rda")
gam_cov <- readRDS(file = "/Users/nmuncy/Desktop/gam_intx.Rda")
summary(gam_cov)
plot_cov <- getViz(gam_cov)
plot(sm(plot_cov, 1))
plot(sm(plot_cov, 2))

# ordered interaction model to get Exp difference
# gam_covOF <- gam_intxOF_model(df_tract, "gamma", "dx_groupOF", "cov")
# saveRDS(gam_covOF, file = "/Users/nmuncy/Desktop/gam_intxOF.Rda")
gam_covOF <- readRDS(file = "/Users/nmuncy/Desktop/gam_intxOF.Rda")
summary(gam_covOF)
plot_covOF <- getViz(gam_covOF)
plot(sm(plot_covOF, 2))
plot(sm(plot_covOF, 3))

# difference is A-B, so invert for more intuitive plots
p <- plot(sm(plot_covOF, 3))
p_data <- p$data$fit
p_data$zI <- -1 * p_data$z

ggplot(data = p_data, aes(x = x, y = y, z = zI)) +
  geom_tile(aes(fill = zI)) +
  geom_contour(colour = "black") +
  scale_fill_viridis(option = "D", name = "Fit FA") +
  labs(x = "Node", y = "Covariate term") +
  ggtitle("Exp group difference interaction smooth")

# use GI method to get access to plot_diff2 (to verify gam_covOF)
# gam_covGI <- bam(dti_fa ~
#    s(dx_group, bs = "re") +
#    s(subjectID, bs = "re") +
#    s(nodeID, bs = "cr", k = 50, m = 2) +
#    te(nodeID, cov, by = dx_group, bs = c("cr", "tp"), k = c(50, 10), m = 1),
#  data = df_tract,
#  family = Gamma(link = "logit"),
#  method = "fREML"
# )
# saveRDS(gam_covGI, file = "/Users/nmuncy/Desktop/gam_GIintx.Rda")
gam_covGI <- readRDS(file = "/Users/nmuncy/Desktop/gam_GIintx.Rda")
plot(gam_covGI)
h_plot <- getViz(gam_covGI)
plot(sm(h_plot, 3))
plot(sm(h_plot, 4))
plot(sm(h_plot, 5))

par(cex = 1, mar = c(5.1, 4.1, 4.1, 3.5))
plot_diff2(
  gam_covGI,
  view=c("nodeID", "cov"),
  comp=list(dx_group=c("Con", "Exp")),
  rm.ranef=T
)


# New intx method ----
gam_new_intx <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 50, m = 1) +
    s(cov, by = dx_group, bs = "tp", k = 10, m = 2) +
    ti(
      nodeID, cov, by = dx_group, bs = c("cr", "tp"), k = c(50, 10), m = 2
    ),
  data = df_tract,
  family = Gamma(link = "logit"),
  method = "fREML",
  discrete = T
)
gam.check(gam_new_intx, rep = 1000)
summary(gam_new_intx)

plot_new <- getViz(gam_new_intx)
plot(sm(plot_new, 1))
plot(sm(plot_new, 3))
plot(sm(plot_new, 4))
plot(sm(plot_new, 5))
plot(sm(plot_new, 6))

# extract plot data to compare
p_con <- plot(sm(plot_new, 5))
p_exp <- plot(sm(plot_new, 6))
df_con <- p_con$data$fit
df_exp <- p_exp$data$fit
df_diff <- df_exp - df_con
df_diff$x <- df_con$x
df_diff$y <- df_con$y
ggplot(data = df_diff, aes(x = x, y = y, z = z)) +
  geom_tile(aes(fill = z)) +
  geom_contour(colour = "black") +
  scale_fill_viridis(option = "D", name = "Fit FA") +
  labs(x = "Covariate Term", y = "NodeID") +
  ggtitle("Exp group difference interaction smooth")


# ordered factor
gam_new_intxOF <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    #s(dx_group, bs = "re") +
    s(nodeID, bs = "cr", k = 50, m = 1) +
    s(cov, by = dx_group, bs = "tp", k = 10, m = 2) +
    ti(nodeID, cov, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
    ti(
      nodeID, cov, by = dx_groupOF, bs = c("cr", "tp"), k = c(50, 10), m = 2
    ),
  data = df_tract,
  family = Gamma(link = "logit"),
  method = "fREML",
  discrete = T
)
summary(gam_new_intxOF)

plot_newOF <- getViz(gam_new_intxOF)
plot(sm(plot_newOF, 1))
plot(sm(plot_newOF, 2))
plot(sm(plot_newOF, 3))
plot(sm(plot_newOF, 4))
plot(sm(plot_newOF, 5))
plot(sm(plot_newOF, 6))


# invert difference smooth to help interpretation
p <- plot(sm(plot_newOF, 6))
p_data <- p$data$fit
p_data$zI <- -1 * p_data$z

ggplot(data = p_data, aes(x = x, y = y, z = zI)) +
  geom_tile(aes(fill = zI)) +
  geom_contour(colour = "black") +
  scale_fill_viridis(option = "D", name = "Fit FA") +
  labs(x = "Covariate Term", y = "NodeID") +
  ggtitle("Exp group difference interaction smooth")




