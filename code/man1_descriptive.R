library("tidyr")
library("ggplot2")
library("mgcv")
library("itsadug")
library("mgcViz")
library("fitdistrplus")
library("viridis")


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
scale_groupB <- runif(num_groupB, 8, 12)

# set up long dataframe
df_long <- as.data.frame(matrix(NA, nrow = 2 * num_groupA * 50, ncol = 6))
colnames(df_long) <- c("subj", "group", "node", "y_scale", "fa", "cov")
df_long$subj <- c(rep(subj_groupA, each = 50), rep(subj_groupB, each = 50))
df_long$group <- c(rep("A", num_groupA * 50), rep("B", num_groupB * 50))
df_long$group <- factor(df_long$group)
df_long$node <- rep(seq(1, 50), 2 * num_groupA)
df_long$y_scale <- c(rep(scale_groupA, each = 50), rep(scale_groupB, each = 50))

# fill dataframe with generate values
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

# get group interaction according to scale_y position
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
gam_GS <- bam(fa ~
    s(subj, bs = "re") +
    s(node, bs = "cr", k = 10, m = 2) +
    s(node, group, bs = "fs", k = 10, m = 2),
  data = df_long,
  family = gaussian(),
  method = "fREML"
)
plot(gam_GS)

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
h_plot <- getViz(gam_groupA)
plot(sm(h_plot, 2))

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
h_plot <- getViz(gam_groupB)
plot(sm(h_plot, 2))

# full model - model how fa is a function of subject, node, cov, and group
gam_cov <- bam(fa ~
    s(subj, bs = "re") +
    te(node, cov, bs = c("cr", "tp"), k = c(10, 10), m = 2) +
    t2(
      node, cov, group,
      bs = c("cr", "tp", "re"),
      k = c(10, 10, 2),
      m = 2,
      full = T
    ),
  data = df_long,
  family = gaussian(),
  method = "fREML"
)
gam.check(gam_cov, rep = 1000)
summary(gam_cov)

# visualize, but this does not really show the group differences
h_plot <- getViz(gam_cov)
plot(sm(h_plot, 2))

# Use ordered factors to show how group B differs in its interaction term
# from group A
df_long$groupOF <- NA
df_long$groupOF <- factor(df_long$group, ordered = T)
gam_covOF <- bam(fa ~
    s(subj, bs = "re") +
    te(node, cov, bs = c("cr", "tp"), k = c(10, 10), m = 2) +
    t2(
      node, cov,
      by = groupOF,
      bs = c("cr", "tp"),
      k = c(10, 10),
      m = 2,
      full = T
    ),
  data = df_long,
  family = gaussian(),
  method = "fREML"
)
gam.check(gam_covOF, rep = 1000)

# plot reference smooth
h_plot <- getViz(gam_covOF)
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
