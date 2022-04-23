library("tidyr")
library("itsadug")
library("mgcv")
library("mgcViz")
library("fitdistrplus")
library("ggplot2")
library("viridis")
library("gridExtra")


# General Notes ----
#
# Simulate data for hypothesis visuals, demonstrate analyses. Also
# conduct quick stats needed for manuscript.

data_dir <- "/Users/nmuncy/Projects/emu_unc/data"
out_dir <- "/Users/nmuncy/Desktop"


# # Get Demographics ------
# df_afq <- read.csv(paste0(data_dir, "/AFQ_dataframe.csv"))
# df_subset <- df_afq[which(df_afq$tractID == "UNC_L" & df_afq$nodeID == 10), ]
# df_subset <- df_node10 %>% drop_na(dx)
# num_subj <- dim(df_subset)[1]
# num_female <- length(which(df_subset$sex == "F"))
# age_avg <- round(mean(df_subset$age), 2)
# age_sd <- round(sd(df_subset$age), 2)


# Simulate Data: Make tract and covariates ----
#
# Simulate data to demonstrate (a) groups smooths differing from the global,
# and (b) groups interacting differently with a continuous covariate.

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
  #   list of values length(n)

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
num_con <- num_exp <- 60
subj_con <- 100 + seq(1:num_con)
subj_exp <- 200 + seq(1:num_exp)
scale_con <- runif(num_con, 7, 10)
scale_exp <- runif(num_exp, 7, 12)

# set up long dataframe
df_long <- as.data.frame(matrix(NA, nrow = 2 * num_con * 50, ncol = 6))
colnames(df_long) <- c("subj", "group", "node", "y_scale", "fa", "cov")
df_long$subj <- c(rep(subj_con, each = 50), rep(subj_exp, each = 50))
df_long$group <- c(rep("Con", num_con * 50), rep("Exp", num_exp * 50))
df_long$group <- factor(df_long$group)
df_long$node <- rep(seq(1, 50), 2 * num_con)
df_long$y_scale <- c(rep(scale_con, each = 50), rep(scale_exp, each = 50))

# fill dataframe with generated values
c_start <- 1
for (h_seed in c(scale_con, scale_exp)) {
  c_end <- c_start + 49
  df_long[c_start:c_end, ]$fa <- gen_norm(h_seed)
  c_start <- c_end + 1
}

# put simulated values in FA range
df_long$fa <- (df_long$fa / 100) + 0.3

# plot node data
# guides(colour = guide_legend(override.aes = list(alpha = 1))) +
pA1 <- ggplot(data = df_long, aes(x = node, y = fa, colour = group)) +
  facet_wrap(~group) +
  geom_point(size = 1, alpha = 0.3) +
  theme(legend.position = "none") +
  labs(x = "Tract Node", y = "Simulated FA") +
  ggtitle("Simulated Group Tracts") +
  theme(
    text = element_text(family = "Times New Roman")
  )
print(pA1)

# generate two covariate distributions for the groups - the con
# dist is normal, exp dist is left skewed with larger mean
set.seed(12)
cov_con <- sort(rbeta(num_con, 50, 30))
cov_exp <- sort(rbeta(num_exp, 8, 2))
# hist(cov_con)
# hist(cov_exp)
# plot(1:num_exp, cov_exp)
# points(1:num_con, cov_con)

# get group interaction according to ordered scale_y position
sorted_scaleA <- sort(scale_con)
sorted_scaleB <- sort(scale_exp)
c <- 1
while (c <= length(subj_con)) {

  # determine subject ids
  subjA <- subj_con[c]
  subjB <- subj_exp[c]

  # find index of y_scale value in sorted scale (position of subj in scale)
  ind_scaleA <- which(grepl(scale_con[c], sorted_scaleA) == T)
  ind_scaleB <- which(grepl(scale_exp[c], sorted_scaleB) == T)

  # determine df_long rows
  ind_dfA <- which(df_long$subj == subjA)
  ind_dfB <- which(df_long$subj == subjB)

  # fill df_long rows with indexed group covariate - so subjs with
  # larger y_scale values also have larger covariates
  df_long[ind_dfA, ]$cov <- cov_con[ind_scaleA]
  df_long[ind_dfB, ]$cov <- cov_exp[ind_scaleB]
  c <- c + 1
}

# plot intx of scaling and cov
df_ind1 <- df_long[which(df_long$node == 1), ]
pB1 <- ggplot(data = df_ind1, aes(x = y_scale, y = cov, color = group)) +
  geom_violin(fill = "gray", colour = "gray", alpha = 0.2) +
  geom_point() +
  facet_wrap(~group) +
  labs(
    x = "Tract Scaling Factor",
    y = "Simulated Covariate",
    colour = "Group"
  ) +
  ggtitle("Tract-Covariate Interaction by Group") +
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none"
  )
print(pB1)
rm(df_ind1)


# Simulate Data: Model Global and Group Smooths ----
#
# Fit data with GS model to demonstrate how Exp/Con groups differ
# in tract.

# gam via GS method to get group smooths
descdist(df_long$fa, discrete = F) # just pretend the uniform dist is gaus
gam_GS <- bam(fa ~
  s(subj, bs = "re") +
  s(node, bs = "cr", k = 30, m = 2) +
  s(node, group, bs = "fs", k = 30, m = 2),
data = df_long,
family = gaussian(),
method = "fREML"
)
gam.check(gam_GS, rep = 1000) # just pretend this is fine
plot_GS <- getViz(gam_GS)

# make, save global tract smooth plot
p <- plot(sm(plot_GS, 2))
p_data <- as.data.frame(p$data$fit)
colnames(p_data) <- c("nodeID", "est", "ty", "se")
p_data$lb <- as.numeric(p_data$est - (2 * p_data$se))
p_data$ub <- as.numeric(p_data$est + (2 * p_data$se))

pA2 <- ggplot(data = p_data, aes(x = .data$nodeID, y = .data$est)) +
  geom_line() +
  geom_ribbon(aes(ymin = .data$lb, ymax = .data$ub), alpha = 0.2) +
  scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  ggtitle("Simulated Tract Smooth") +
  ylab("Est. FA Fit") +
  xlab("Tract Node") +
  theme(text = element_text(family = "Times New Roman"))
print(pA2)

# make, save group smooths plot
p <- plot(sm(plot_GS, 3))
p_data <- as.data.frame(p$data$fit)
colnames(p_data) <- c("nodeID", "est", "ty", "Group")

# make, save ggplot
pA3 <- ggplot(
  data = p_data,
  aes(x = .data$nodeID, y = .data$est, group = .data$Group)
) +
  facet_wrap(~Group) +
  geom_line(aes(color = .data$Group)) +
  scale_y_continuous(limits = c(-0.05, 0.05)) +
  scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  ggtitle("Simulated Group Smooths") +
  labs(y = "Est. FA Fit", x = "Tract Node") +
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none"
  )
print(pA3)

# clean
rm(gam_GS)
rm(p)
rm(p_data)
rm(plot_GS)



# Simulate Data: Model Node-FA-Cov Interactions ----
#
# Model the interaction of the simulated tracts and covariate by group.
gam_cov <- bam(fa ~
  s(subj, bs = "re") +
  s(node, bs = "cr", k = 30, m = 2) +
  s(cov, by = group, bs = "tp", k = 10, m = 2) +
  ti(
    node, cov,
    by = group, bs = c("cr", "tp"), k = c(30, 10), m = 2
  ),
data = df_long,
family = gaussian(),
method = "fREML",
discrete = T
)
gam.check(gam_cov, rep = 1000)
summary(gam_cov)

# make covariate sequences for each group
cov_con_seq <- seq(min(cov_con), max(cov_con), length = 60)
cov_exp_seq <- seq(min(cov_exp), max(cov_exp), length = 60)

# predict control cov smooth
df_con_cov <- df_long[which(df_long$group == "Con" & df_long$node == 25), ]
df_pred_con_cov <- data.frame(
  subj = df_con_cov$subj,
  group = df_con_cov$group,
  node = df_con_cov$node,
  cov = cov_con_seq
)
pred_con_cov <- predict(gam_cov, df_pred_con_cov, se.fit = T)
pred_con_cov <- transform(
  pred_con_cov, ub = fit + (2*se.fit), lb = fit - (2*se.fit)
)
df_pred_con_cov <- cbind(df_pred_con_cov, pred_con_cov)

# predict experimental cov smooth
df_exp_cov <- df_long[which(df_long$group == "Exp" & df_long$node == 25), ]
df_pred_exp_cov <- data.frame(
  subj = df_exp_cov$subj,
  group = df_exp_cov$group,
  node = df_exp_cov$node,
  cov = cov_exp_seq
)
pred_exp_cov <- predict(gam_cov, df_pred_exp_cov, se.fit = T)
pred_exp_cov <- transform(
  pred_exp_cov, ub = fit + (2*se.fit), lb = fit - (2*se.fit)
)
df_pred_exp_cov <- cbind(df_pred_exp_cov, pred_exp_cov)

# stitch prediction dfs together, draw plot
df_pred_cov <- rbind(df_pred_con_cov, df_pred_exp_cov)
pB2 <- ggplot(data = df_pred_cov, aes(x = cov, y = fit)) +
  geom_line(aes(color = group)) +
  facet_wrap(~group) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2) +
  labs(x = "Covariate Value", y = "Est. FA Fit") +
  ggtitle("Simulated Covariate Smooths") +
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none"
  )
print(pB2)

# clean
rm(df_con_cov)
rm(df_exp_cov)
rm(pred_con_cov)
rm(pred_exp_cov)
rm(df_pred_con_cov)
rm(df_pred_exp_cov)
rm(df_pred_cov)

# predict node-fa-cov interaction, Con group
df_con <- df_long[which(df_long$group == "Con"), ]
df_pred_con_intx <- data.frame(
  subj = df_con$subj,
  group = df_con$group,
  node = df_con$node,
  cov = rep(cov_con_seq, each = 50)
)
df_pred_con_intx$fit <- predict(gam_cov, df_pred_con_intx)

pC1 <- ggplot(df_pred_con_intx, aes(x = node, y = cov, z = fit)) +
  geom_tile(aes(fill = fit)) +
  geom_contour(colour = "black") +
  scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  scale_fill_viridis(
    option = "D", 
    name = "Est. FA Fit", 
    limits = c(0.3, 0.8)
  ) +
  labs(y = "Simulated Covariate", x = "Tract Node") +
  ggtitle("Tract Node-FA-Covariate Interaction, Con") +
  theme(
    text = element_text(family = "Times New Roman")
  )
print(pC1)

# predict node-fa-cov interaction, Exp group
df_exp <- df_long[which(df_long$group == "Exp"), ]
df_pred_exp_intx <- data.frame(
  subj = df_exp$subj,
  group = df_exp$group,
  node = df_exp$node,
  cov = rep(cov_exp_seq, each = 50)
)
df_pred_exp_intx$fit <- predict(gam_cov, df_pred_exp_intx)

pC2 <- ggplot(df_pred_exp_intx, aes(x = node, y = cov, z = fit)) +
  geom_tile(aes(fill = fit)) +
  geom_contour(colour = "black") +
  scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  scale_fill_viridis(
    option = "D", 
    name = "Est. FA Fit", 
    limits = c(0.3, 0.8)
  ) +
  labs(y = "Simulated Covariate", x = "Tract Node") +
  ggtitle("Tract Node-FA-Covariate Interaction, Exp") +
  theme(
    text = element_text(family = "Times New Roman")
  )
print(pC2)

# clean
rm(df_con)
rm(df_exp)
rm(df_pred_con_intx)
rm(df_pred_exp_intx)
rm(gam_cov)

# use ordered factors to show how group B differs in its interaction term
# from group A, yields an F-stat
df_long$groupOF <- NA
df_long$groupOF <- factor(df_long$group, ordered = T)
gam_covOF <- bam(fa ~
  s(subj, bs = "re") +
  s(node, bs = "cr", k = 30, m = 2) +
  s(cov, by = group, bs = "tp", k = 10, m = 2) +
  ti(node, cov, bs = c("cr", "tp"), k = c(30, 10), m = 2) +
  ti(
    node, cov,
    by = groupOF, bs = c("cr", "tp"), k = c(30, 10), m = 2
  ),
data = df_long,
family = gaussian(),
method = "fREML",
discrete = T
)
gam.check(gam_covOF, rep = 1000)
summary(gam_covOF)

# predict node-fa-cov difference interaction
df_exp <- df_long[which(df_long$group == "Exp"), ]
df_pred_diff <- data.frame(
  subj = df_exp$subj,
  group = df_exp$group,
  groupOF = df_exp$groupOF,
  node = df_exp$node,
  cov = rep(cov_exp_seq, each = 50)
)
pred_intx_diff <- as.data.frame(predict.gam(
  gam_covOF, df_pred_diff, type = "terms"
))
df_pred_diff <- cbind(df_pred_diff, pred_intx_diff[6])
colnames(df_pred_diff) <- c(
  "subj", "group", "groupOF", "node", "cov", "fit"
)
pC3 <- ggplot(df_pred_diff, aes(x = node, y = cov, z = fit)) +
  geom_tile(aes(fill = fit)) +
  geom_contour(colour = "black") +
  scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  scale_fill_viridis(option = "D", name = "Est. FA Fit") +
  labs(y = "Simulated Covariate", x = "Tract Node") +
  ggtitle("Tract Node-FA-Covariate Interaction, Diff") +
  theme(
    text = element_text(family = "Times New Roman")
  )
print(pC3)

# calc cov group diff smooth
p_data <- plot_diff(
  gam_covOF,
  view = c("cov"),
  comp = list(group = c("Exp", "Con")),
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
sig_nodes <- p_data[sig_rows, ]$cov

# find start, end points of sig regions
vec_start <- sig_nodes[1]
vec_end <- vector()
y_min <- min(p_data$lb)
num_rows <- length(sig_rows)
c <- 2
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
pB3 <- ggplot(data = p_data, aes(x = cov, y = est)) +
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
  ggtitle("Tract Covariate Smooth, Diff") +
  theme(
    text = element_text(family = "Times New Roman")
  )
print(pB3)


# make multipane plot for figure
pAll <- grid.arrange(
  pA1, pB1, pC1, 
  pA2, pB2, pC2, 
  pA3, pB3, pC3,
  nrow = 3,
  ncol = 3
  )

ggsave(
  paste0(out_dir, "/sim_data.png"),
  plot = pAll,
  units = "in",
  width = 12,
  height = 9,
  dpi = 600,
  device = "png"
)
