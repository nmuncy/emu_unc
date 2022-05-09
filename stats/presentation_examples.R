library("ggplot2")
library("mgcv")
library("dplyr")
library("tidyr")
library("gridExtra")

# Set Up ----
proj_dir <- "/Users/nmuncy/Projects/emu_unc"
data_dir <- paste0(proj_dir, "/data")
out_dir <- "/Users/nmuncy/Desktop"

# import data, setup factors
df_afq <- read.csv(paste0(data_dir, "/AFQ_dataframe.csv"))
df_afq$sex <- factor(df_afq$sex)
ind_exp <- which(df_afq$dx_group == "Pat")
df_afq[ind_exp, ]$dx_group <- "Exp"
df_afq$dx_group <- factor(df_afq$dx_group)
df_afq$subjectID <- factor(df_afq$subjectID)
df_afq$dx_groupOF <- factor(df_afq$dx_group, ordered = T)

# clip tails
ind_keep <- which(
  df_afq$nodeID >= 10 & df_afq$nodeID <= 89
)
df_afq <- df_afq[ind_keep, ]
rm(ind_keep)
rm(ind_exp)


# Plot tract points ----
#
# Show raw data, linear fit, GAM smooths
df_tract <- df_afq[which(df_afq$tractID == "UNC_L"), ]
df_tract <- df_tract %>% drop_na(dx)

p1 <- ggplot(df_tract, aes(x = nodeID, y = dti_fa)) +
  geom_point(size = 0.1) +
  labs(
    x = "Node", 
    y = "Fractional Anisotropy", 
    title = "L. Uncinate - Data Points"
    ) +
  scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
  theme(text = element_text(family = "Times New Roman"))
ggsave(
  paste0(out_dir, "/Plot_raw_tract.png"),
  plot = p1,
  units = "in",
  height = 3,
  width = 4,
  dpi = 300,
  device = "png"
)

p2 <- ggplot(df_tract, aes(x = nodeID, y = dti_fa)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", formula = y ~ x, color = "blue") +
  labs(
    x = "Node", 
    y = "Fractional Anisotropy", 
    title = "L. Uncinate - Linear Model"
    ) +
  scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
  theme(text = element_text(family = "Times New Roman"))
ggsave(
  paste0(out_dir, "/Plot_linear_tract.png"),
  plot = p2,
  units = "in",
  height = 3,
  width = 4,
  dpi = 300,
  device = "png"
)

p3 <- ggplot(df_tract, aes(x = nodeID, y = dti_fa)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "blue") +
  labs(x = "Node", y = "Fractional Anisotropy", title = "L. Uncinate - GAM") +
  scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
  theme(text = element_text(family = "Times New Roman"))
ggsave(
  paste0(out_dir, "/Plot_gam_tract.png"),
  plot = p3,
  units = "in",
  height = 3,
  width = 4,
  dpi = 300,
  device = "png"
)

p4 <- ggplot(df_tract, aes(x = nodeID, y = dti_fa, color = dx_group)) +
  geom_point(size = 0.1, color = "black") +
  geom_smooth(method = "gam", formula = y ~ s(x)) +
  labs(
    x = "Node", 
    y = "Fractional Anisotropy", 
    title = "L. Uncinate - GAM by Group", 
    color = "Group"
    ) +
  scale_x_continuous(breaks = c(seq(10, 89, by = 10), 89)) +
  scale_color_manual(values = c("red", "blue")) +
  theme(
    text = element_text(family = "Times New Roman"), 
    legend.position = c(0.9, 0.8)
    )
ggsave(
  paste0(out_dir, "/Plot_gamm_tract.png"),
  plot = p4,
  units = "in",
  height = 3,
  width = 4,
  dpi = 300,
  device = "png"
)

pOut <- grid.arrange(
  p1, p2,
  p3, p4,
  nrow = 2,
  ncol = 2
)
ggsave(
  paste0(out_dir, "/Plot_raw_data.png"),
  plot = pOut,
  units = "in",
  height = 6,
  width = 9,
  dpi = 300,
  device = "png"
)


# Simulate data for basis function illustration ----
set.seed(1)
x <- seq(0, pi * 2, 0.1)
sin_x <- sin(x)
y <- sin_x + rnorm(n = length(x), mean = 0, sd = sd(sin_x / 2))
df_sim <- data.frame(y, x)

# gam fit for thin plates, cubic splines
gam_tp <- gam(y ~ s(x, bs = "tp", k = 5), method = "REML")
gam_pred <- predict(gam_tp, newdata = data.frame(x = df_sim[, 2]))
basis_tp <- predict(gam_tp, type = "lpmatrix")

gam_cr <- gam(y ~ s(x, bs = "cr", k = 5), method = "REML")
gam_pred <- predict(gam_cr, newdata = data.frame(x = df_sim[, 2]))
basis_cr <- predict(gam_cr, type = "lpmatrix")

# thin plates
# png(
#   file = paste0(out_dir, "/Plot_sim_data.png"),
#   width = 12, height = 4, units = "in", res = 300
# )
par(family = "Times New Roman")
# par(mfrow = c(1, 4))
plot(y ~ x, main = "Distribution of Points")

plot(y ~ x, main = "GAM Fit")
lines(x, gam_pred, col = "blue", lwd = 2)

plot(y ~ x, main = "GAM Fit & Basis Functions, Thin Plate")
lines(x, gam_pred, col = "blue", lwd = 2)
matplot(x, basis_tp[, -1], type = "l", lty = 1, add = T)

# cubic splines
plot(y ~ x, main = "GAM Fit & Basis Functions, Cubic")
lines(x, gam_pred, col = "blue", lwd = 2)
matplot(x, basis_cr[, -1], type = "l", lty = 1, add = T)
# dev.off()
