library("tidyr")
library("itsadug")
library("fitdistrplus")
source("./diff4_calc_gams.R")
source("./diff4_plot_gams.R")


# General Notes ----
#
# Simulate data for hypothesis visuals, demonstrate analyses. Also
# conduct quick stats needed for manuscript.


# Read in AFQ data ----
data_dir <- "/Users/nmuncy/Projects/emu_unc/data"
df_afq <- read.csv(paste0(data_dir, "/AFQ_dataframe.csv"))


# Get Demographics ------
df_subset <- df_afq[which(df_afq$tractID == "UNC_L" & df_afq$nodeID == 10), ]
df_subset <- df_node10 %>% drop_na(dx)
num_subj <- dim(df_subset)[1]
num_female <- length(which(df_subset$sex == "F"))
age_avg <- round(mean(df_subset$age), 2)
age_sd <- round(sd(df_subset$age), 2)


# Simulate data to demonstrate an interaction smooth ----
#
# Simulate data to demonstrate (a) groups smooths differing from the global,
# and (b) groups interacting differently with a continuous covariate.
#
# May be outdated, I spent more time working with the real FA values and
# simulated covariates.
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
  geom_point() +
  ggtitle("Two simulated tracts")

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
  geom_point() +
  ggtitle("Two simulated interactions with a covariate")

# gam via GS method to get group smooths, plot relevant attributes
descdist(df_long$fa, discrete = F) # just pretend the uniform dist is gaus
gam_GS <- bam(fa ~
s(subj, bs = "re") +
  s(node, bs = "cr", k = 10, m = 2) +
  s(node, group, bs = "fs", k = 10, m = 2),
data = df_long,
family = gaussian(),
method = "fREML"
)
gam.check(gam_GS, rep = 1000) # just pretend this is fine
plot(gam_GS)
h_plot <- getViz(gam_GS)
plot(sm(h_plot, 2)) + ggtitle("Global smooth of simulated tract")
plot(sm(h_plot, 3)) + ggtitle("Group smooths of simulated tract")

# use ordered factor to make group difference smooth, plot
df_long$groupOF <- ordered(df_long$group)
gam_OF <- bam(fa ~
s(subj, bs = "re") +
  s(node, bs = "cr", k = 10, m = 2) +
  s(node, by = groupOF, bs = "cr", k = 10, m = 2),
data = df_long,
family = gaussian(),
method = "fREML"
)
h_plot <- getViz(gam_OF)
plot(sm(h_plot, 3)) + ggtitle("B - A difference smooth")

# model cov-fa interaction of groupA - visualize linear fa-cov intx,
# make a contour plot, 3d plot
df_groupA <- df_long[which(df_long$group == "A"), ]
gam_groupA <- bam(fa ~
s(subj, bs = "re") +
  te(node, cov, bs = c("cr", "tp"), k = c(10, 10)),
data = df_groupA,
family = gaussian(),
method = "fREML"
)
summary(gam_groupA)
h_plot <- getViz(gam_groupA)
plot(sm(h_plot, 2)) + ggtitle("Group A node FA - cov interaction")
# vis.gam(gam_groupA,
#   view = c("node", "cov"),
#   color = "topo",
#   theta = 330
# )

# model cov-fa interaction of groupB - visualize exponential fa-cov intx,
# make a contour plot
df_groupB <- df_long[which(df_long$group == "B"), ]
gam_groupB <- bam(fa ~
s(subj, bs = "re") +
  te(node, cov, bs = c("cr", "tp"), k = c(10, 10)),
data = df_groupB,
family = gaussian(),
method = "fREML"
)
summary(gam_groupB)
h_plot <- getViz(gam_groupB)
plot(sm(h_plot, 2)) + ggtitle("Group B node FA - cov interaction")

# model how fa is a function of subject, node, cov, and group
gam_cov <- bam(fa ~
s(subj, bs = "re") +
  s(node, bs = "cr", k = 10, m = 2) +
  s(cov, by = group, bs = "tp", k = 10, m = 2) +
  ti(
    node, cov,
    by = group, bs = c("cr", "tp"), k = c(10, 10), m = 2
  ),
data = df_long,
family = gaussian(),
method = "fREML",
discrete = T
)
summary(gam_cov)

# visualize the various gam attributes
h_plot <- getViz(gam_cov)
plot(sm(h_plot, 2))
plot(sm(h_plot, 3))
plot(sm(h_plot, 4))
plot(sm(h_plot, 5))
plot(sm(h_plot, 6))

# extract plotting data to calculate difference contour of groupB fa-cov
# interaction from groupA
p_con <- plot(sm(h_plot, 5))
p_exp <- plot(sm(h_plot, 6))
df_con <- p_con$data$fit
df_exp <- p_exp$data$fit
df_diff <- df_con - df_exp
df_diff$x <- df_con$x
df_diff$y <- df_con$y
ggplot(data = df_diff, aes(x = x, y = y, z = z)) +
  geom_tile(aes(fill = z)) +
  geom_contour(colour = "black") +
  scale_fill_viridis(option = "D", name = "Fit FA") +
  labs(x = "node", y = "cov") +
  ggtitle("Group B difference FA - cov interaction smooth (manual)")

# use ordered factors to show how group B differs in its interaction term
# from group A, yields an F-stat
df_long$groupOF <- NA
df_long$groupOF <- factor(df_long$group, ordered = T)
gam_covOF <- bam(fa ~
s(subj, bs = "re") +
  s(node, bs = "cr", k = 10, m = 1) +
  s(cov, by = group, bs = "tp", k = 10, m = 2) +
  ti(node, cov, bs = c("cr", "tp"), k = c(10, 10), m = 2) +
  ti(
    node, cov,
    by = groupOF, bs = c("cr", "tp"), k = c(10, 10), m = 2
  ),
data = df_long,
family = gaussian(),
method = "fREML",
discrete = T
)
summary(gam_covOF)

# visualize
h_plot <- getViz(gam_covOF)
plot(sm(h_plot, 3))
plot(sm(h_plot, 6))

# difference is A-B, so invert for more intuitive plots
p <- plot(sm(h_plot, 6))
p_data <- p$data$fit
p_data$zI <- -1 * p_data$z

ggplot(data = p_data, aes(x = x, y = y, z = zI)) +
  geom_tile(aes(fill = zI)) +
  geom_contour(colour = "black") +
  scale_fill_viridis(option = "D", name = "Fit FA") +
  labs(x = "node", y = "cov") +
  ggtitle("Group B difference FA - cov interaction smooth (ordered GAM)")



# Use real FA data, simulate interaction ----
#
# Using actual L. Unc data, simulate a group x covariate interaction
# at a region that differs between Experimental/Control groups.

# get tract data, clean up
df_tract <- df_afq[which(df_afq$tractID == "UNC_L"), ]
df_tract <- df_tract %>% drop_na(dx)
ind_pat <- which(df_tract$dx_group == "Pat")
df_tract[ind_pat, ]$dx_group <- "Exp"
df_tract$dx_group <- factor(df_tract$dx_group)
df_tract$subjectID <- factor(df_tract$subjectID)
df_tract$dx_groupOF <- factor(df_tract$dx_group, ordered = T)

# clip tails
ind_keep <- which(
  df_tract$nodeID >= 10 & df_tract$nodeID <= 89
)
df_tract <- df_tract[ind_keep, ]

# determine num of subjects in each group
num_groupC <- length(which(df_tract$nodeID == 50 & df_tract$dx_group == "Con"))
num_groupE <- length(which(df_tract$nodeID == 50 & df_tract$dx_group == "Exp"))

# generate linear covariate values for control, exponential
# for experimental group
set.seed(12)
seq_C <- 1:num_groupC
seq_E <- 1:num_groupE
lin_cov <- (0.1 * seq_C + 2) + rnorm(length(seq_C), 0, 0.1)
exp_cov <- (((0.1 * seq_E)^2) + 2) + rnorm(length(seq_E), 0, 0.1)

# groups differ at nodes 30-45, use ordered node 37 to assign simulated
# covariate values -- participants with larger 37 FA value have larger
# cov value
fa_37C <- sort(df_tract[which(
  df_tract$dx_group == "Con" & df_tract$nodeID == 37
), ]$dti_fa)

fa_37E <- sort(df_tract[which(
  df_tract$dx_group == "Exp" & df_tract$nodeID == 37
), ]$dti_fa)

# assign simulated covariates for con/exp subjs
df_tract$cov <- NA
subj_C <- unique(df_tract[which(df_tract$dx_group == "Con"), ]$subjectID)
subj_E <- unique(df_tract[which(df_tract$dx_group == "Exp"), ]$subjectID)

for (subj in subj_C) {

  # get fa value at node 3d
  val_37 <- df_tract[which(
    df_tract$subjectID == subj & df_tract$nodeID == 37
  ), ]$dti_fa

  # determine position of subj of node 37 fa value in sorted list
  ord_37 <- which(grepl(val_37, fa_37C) == T)
  ind_tract <- which(df_tract$subjectID == subj)

  # get covariate from same position
  df_tract[ind_tract, ]$cov <- lin_cov[ord_37]
}

# repeat for exp group
for (subj in subj_E) {
  val_37 <- df_tract[which(
    df_tract$subjectID == subj & df_tract$nodeID == 37
  ), ]$dti_fa
  ord_37 <- which(grepl(val_37, fa_37E) == T)
  ind_tract <- which(df_tract$subjectID == subj)
  df_tract[ind_tract, ]$cov <- exp_cov[ord_37]
}

# plot intx of group, cov, and node37 FA
df_node37 <- df_tract[which(df_tract$nodeID == 37), ]
ggplot(data = df_node37, aes(x = dti_fa, y = cov, color = dx_group)) +
  geom_point() +
  labs(y = "Simulated Covariate", x = "Node 37 FA", colour = "Group") +
  ggtitle("L. Uncinate Group-Covariate Interactions") +
  theme(text = element_text(family = "Times New Roman"))
ggsave(
  filename = "/Users/nmuncy/Desktop/sim_intx.png",
  plot = last_plot(),
  units = "in",
  width = 4,
  height = 3,
  device = "png"
)
rm(df_node37)


# Model Global and Group smooths ----
#
# Conduct GS style GAM to show that groups differ in their smooths.

# gam via GS method, make getViz plot object
descdist(df_tract$dti_fa, discrete = F)
gam_GS <- gam_GS(df_tract, "gamma", "dx_group")
plot_GS <- getViz(gam_GS)

# make, save global tract smooth plot
p <- plot(sm(plot_GS, 2))
p_data <- as.data.frame(p$data$fit)
colnames(p_data) <- c("nodeID", "est", "ty", "se")
p_data$lb <- as.numeric(p_data$est - (2 * p_data$se))
p_data$ub <- as.numeric(p_data$est + (2 * p_data$se))

# sequence range (10-89) results from clipping of tails
ggplot(data = p_data, aes(x = .data$nodeID, y = .data$est)) +
  geom_line() +
  geom_ribbon(aes(ymin = .data$lb, ymax = .data$ub), alpha = 0.2) +
  scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
  ggtitle("L. Uncinate Tract Smooth") +
  ylab("Est. FA Fit") +
  xlab("Tract Node") +
  theme(text = element_text(family = "Times New Roman"))
ggsave(
  filename = "/Users/nmuncy/Desktop/lunc_global.png",
  plot = last_plot(),
  units = "in",
  width = 4,
  height = 3,
  device = "png"
)

# make, save group smooths plot
p <- plot(sm(plot_GS, 3))
p_data <- as.data.frame(p$data$fit)
colnames(p_data) <- c("nodeID", "est", "ty", "Group")

# make, save ggplot
ggplot(
  data = p_data,
  aes(x = .data$nodeID, y = .data$est, group = .data$Group)
) +
  geom_line(aes(color = .data$Group)) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
  ggtitle("L. Uncinate Group Smooths") +
  labs(y = "Est. FA Fit", x = "Tract Node") +
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none"
  )
ggsave(
  filename = "/Users/nmuncy/Desktop/lunc_group.png",
  plot = last_plot(),
  units = "in",
  width = 4,
  height = 3,
  device = "png"
)

# gam via ordered factor method to produce difference smooth
gam_GSOF <- gam_GSOF(df_tract, "gamma", "dx_groupOF")
plot_GSOF <- getViz(gam_GSOF)

# unpack difference smooth data for pretty plotting
p <- plot(sm(plot_GSOF, 3)) +
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

# draw difference smooth with differing nodes highlighted
ggplot(data = p_data, aes(x = .data$nodeID, y = .data$est)) +
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
  ggtitle("L. Uncinate Exp-Con Difference Smooth") +
  ylab("Est. Difference") +
  xlab("Tract Node") +
  theme(text = element_text(family = "Times New Roman"))
ggsave(
  filename = "/Users/nmuncy/Desktop/lunc_group-diff.png",
  plot = last_plot(),
  units = "in",
  width = 4,
  height = 3,
  device = "png"
)


# Model individual group-covariate interactions ----
#
# Model each group separately to check the interaction of node-fa-cov.

# interaction of groupC - visualize linear fa-cov intx
df_groupA <- df_tract[which(df_tract$dx_group == "Con"), ]
gam_groupA <- bam(dti_fa ~
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(cov, bs = "tp", k = 5, m = 2) +
  ti(nodeID, cov, bs = c("cr", "tp"), k = c(50, 5), m = 2),
data = df_groupA,
family = Gamma(link = "logit"),
method = "fREML",
discrete = T
)
summary(gam_groupA)

# make contour plot
plot_groupA <- getViz(gam_groupA)
p <- plot(sm(plot_groupA, 4))
p_data <- p$data$fit
colnames(p_data) <- c("z", "tz", "cov", "node", "se")

ggplot(
  data = p_data, aes(x = .data$node, y = .data$cov, z = .data$z)
) +
  geom_tile(aes(fill = .data$z)) +
  geom_contour(colour = "black") +
  scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
  scale_fill_viridis(option = "D", name = "Est. FA Fit") +
  labs(y = "Simulated Covariate", x = "Tract Node") +
  ggtitle("L. Uncinate Node-FA-Covariate Interaction, Con") +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(size = 12)
  )
ggsave(
  filename = "/Users/nmuncy/Desktop/lunc_gA-intx.png",
  plot = last_plot(),
  units = "in",
  width = 4,
  height = 3,
  device = "png"
)

# interaction of groupE - visualize exponential fa-cov intx
df_groupB <- df_tract[which(df_tract$dx_group == "Exp"), ]
gam_groupB <- bam(dti_fa ~
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(cov, bs = "tp", k = 5, m = 2) +
  ti(nodeID, cov, bs = c("cr", "tp"), k = c(50, 5), m = 2),
data = df_groupB,
family = Gamma(link = "logit"),
method = "fREML",
discrete = T
)
summary(gam_groupB)

# make contour plot
plot_groupB <- getViz(gam_groupB)
p <- plot(sm(plot_groupB, 4))
p_data <- p$data$fit
colnames(p_data) <- c("z", "tz", "cov", "node", "se")

ggplot(
  data = p_data, aes(x = .data$node, y = .data$cov, z = .data$z)
) +
  geom_tile(aes(fill = .data$z)) +
  geom_contour(colour = "black") +
  scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
  scale_fill_viridis(option = "D", name = "Est. FA Fit") +
  labs(y = "Simulated Covariate", x = "Tract Node") +
  ggtitle("L. Uncinate Node-FA-Covariate Interaction, Exp") +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(size = 12)
  )
ggsave(
  filename = "/Users/nmuncy/Desktop/lunc_gB-intx.png",
  plot = last_plot(),
  units = "in",
  width = 4,
  height = 3,
  device = "png"
)


# Model group-cov interaction ----
#
# Model the full dataset, determine how experimental group differs from
# control group.

# full interaction model, visualize
gam_cov <- gam_GSintx(df_tract, "gamma", "dx_group", "cov")
plot_gam_cov <- getViz(gam_cov)

# make node-cov-fa interaction plots for each group (con/exp)
p <- plot(sm(plot_gam_cov, 5))
p_data <- p$data$fit
colnames(p_data) <- c("z", "tz", "cov", "node", "se")

ggplot(
  data = p_data, aes(x = .data$node, y = .data$cov, z = .data$z)
) +
  geom_tile(aes(fill = .data$z)) +
  geom_contour(colour = "black") +
  scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
  scale_fill_viridis(option = "D", name = "Est. FA Fit") +
  labs(y = "Simulated Covariate", x = "Tract Node") +
  ggtitle("L. Uncinate Node-FA-Covariate Interaction, Con") +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(size = 12)
  )
ggsave(
  filename = "/Users/nmuncy/Desktop/lunc_gA-intx.png",
  plot = last_plot(),
  units = "in",
  width = 4,
  height = 3,
  device = "png"
)

p <- plot(sm(plot_gam_cov, 6))
p_data <- p$data$fit
colnames(p_data) <- c("z", "tz", "cov", "node", "se")

ggplot(
  data = p_data, aes(x = .data$node, y = .data$cov, z = .data$z)
) +
  geom_tile(aes(fill = .data$z)) +
  geom_contour(colour = "black") +
  scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
  scale_fill_viridis(option = "D", name = "Est. FA Fit") +
  labs(y = "Simulated Covariate", x = "Tract Node") +
  ggtitle("L. Uncinate Node-FA-Covariate Interaction, Exp") +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(size = 12)
  )
ggsave(
  filename = "/Users/nmuncy/Desktop/lunc_gB-intx.png",
  plot = last_plot(),
  units = "in",
  width = 4,
  height = 3,
  device = "png"
)

# ordered interaction model to get Exp difference smooth
gam_covOF <- gam_GSintxOF(df_tract, "gamma", "dx_group", "dx_groupOF", "cov")
summary(gam_covOF)
plot_gam_covOF <- getViz(gam_covOF)

# plot node-fa-cov exp diff (from con) smooth, invert for ease of
# interpretation
p <- plot(sm(plot_gam_covOF, 6))
p_data <- p$data$fit
p_data$zI <- -1 * p_data$z
colnames(p_data) <- c("z", "tz", "cov", "node", "se", "zI")

ggplot(
  data = p_data, aes(x = .data$node, y = .data$cov, z = .data$zI)
) +
  geom_tile(aes(fill = .data$zI)) +
  geom_contour(colour = "black") +
  scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
  scale_fill_viridis(option = "D", name = "Est. FA Fit") +
  labs(y = "Simulated Covariate", x = "Tract Node") +
  ggtitle("L. Uncinate Node-FA-Covariate Interaction, Diff") +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(size = 12)
  )
ggsave(
  filename = "/Users/nmuncy/Desktop/lunc_gAB-diff.png",
  plot = last_plot(),
  units = "in",
  width = 4,
  height = 3,
  device = "png"
)

# calc cov group diff smooth
p_data <- plot_diff(
  gam_covOF,
  view = c("h_var"),
  comp = list(h_group = c("Con", "Exp")),
  rm.ranef = T
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

# draw smooth, shade diff regions
ggplot(data = p_data, aes(x = h_var, y = est)) +
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
  labs(x = "Simulated Covariate", y = "Est. FA Fit") +
  ggtitle("L. Uncinate Covariate Smooth, Diff") +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(size = 12)
  )
ggsave(
  filename = "/Users/nmuncy/Desktop/lunc_cov-diff.png",
  plot = last_plot(),
  units = "in",
  width = 4,
  height = 3,
  dpi = 600,
  device = "png"
)
