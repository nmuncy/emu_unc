library("ggplot2")
library("fitdistrplus")
library("mgcv")
library("itsadug")
library("tidymv")
library("dplyr")
library("mgcViz")
library("tools")
library("viridis")
library("cowplot")


# Set Up ----
data_dir <- file_path_as_absolute(paste0(getwd(), "/../data"))
# capture.output(sessionInfo(), file = paste0(data_dir, "R_session_info.txt"))
df_afq <- read.csv(paste0(data_dir, "/AFQ_dataframe.csv"))
df_afq$sex <- factor(df_afq$sex)
df_afq$dx_group <- factor(df_afq$dx_group)
df_afq$subjectID <- factor(df_afq$subjectID)


# L. Unc ----
#
#

# subset df_afq, take complete cases, make factors
df_tract <- df_afq[which(df_afq$tractID == "UNC_L"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

# determine distribution
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete = F)

# build gam
lunc_gaus <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_gaus, rep = 1000)
summary(lunc_gaus)
plot(lunc_gaus)

# pds effect
lunc_pds <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50) +
  s(pds, by = sex),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_pds, rep = 1000)
summary(lunc_pds)
plot(lunc_pds)
compareML(lunc_gaus, lunc_pds) # gaus preferred

# L. Unc: GS
lunc_dxGS <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_dxGS, rep = 1000)
compareML(lunc_gaus, lunc_dxGS) # lunc_dxGS preferred
summary(lunc_dxGS)
plot(lunc_dxGS)

# L. Unc: GI
lunc_dxGI <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(dx_group, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_group, bs = "cr", k = 50, m = 1),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_dxGI, rep = 1000)
summary(lunc_dxGI)
plot(lunc_dxGI)
compareML(lunc_dxGI, lunc_dxGS) # lunc_dxGS preferred

# L. Unc: GS difference via ordered factors
df_tract$dx_groupOF <- factor(df_tract$dx_group, ordered = T)
lunc_dxGS_OF <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_groupOF, bs = "cr", k = 50, m = 2),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_dxGS_OF, rep = 1000)
plot(lunc_dxGS_OF)
summary(lunc_dxGS_OF)
plot_lunc_dxGS_OF <- getViz(lunc_dxGS_OF)
plot(sm(plot_lunc_dxGS_OF, 2))
plot(sm(plot_lunc_dxGS_OF, 3))

# unpack difference smooth data
p <- plot(sm(plot_lunc_dxGS_OF, 3)) +
  geom_hline(yintercept = 0)
p_data <- as.data.frame(p$data$fit)

# find sig nodes
p_data$lb <- as.numeric(p_data$y - (2 * p_data$se))
p_data$ub <- as.numeric(p_data$y + (2 * p_data$se))
sig_nodes <- which(
  (p_data$y < 0 & p_data$ub < 0) |
    (p_data$y > 0 & p_data$lb > 0)
)

# find start, end points
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

# make df for drawing rectangles, adjust for 0-index nodeID
d_rect <- data.frame(
  x_start = vec_start,
  x_end = vec_end,
  y_start = rep(y_min, length(vec_start)),
  y_end = rep(0, length(vec_start))
)
d_rect$x_start <- d_rect$x_start - 1
d_rect$x_end <- d_rect$x_end - 1

# draw
p + annotate(
  "rect",
  xmin = c(d_rect$x_start),
  xmax = c(d_rect$x_end),
  ymin = c(d_rect$y_start),
  ymax = c(d_rect$y_end),
  alpha = 0.2,
  fill = "red"
) +
  ggtitle("Difference Smooth, Patient-Control") +
  ylab("Est. Difference") +
  xlab("Tract Node")

# L. Unc: GS intx neg LGI
# s(nodeID, bs = "cr", k = 50, m = 2) +
# s(nodeID, dx_group, bs = "fs", k = 50, m = 2) +
lunc_dxGS_neg <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  te(nodeID, lgi_neg, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
  t2(
    nodeID, lgi_neg, dx_group,
    bs = c("cr", "tp", "re"),
    k = c(50, 10, 2),
    m = 2,
    full = TRUE
  ),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_dxGS_neg, rep = 1000)
compareML(lunc_dxGS, lunc_dxGS_neg)
summary(lunc_dxGS_neg)
plot(lunc_dxGS_neg)
plot_lunc_dxGS_neg <- getViz(lunc_dxGS_neg)
plot(sm(plot_lunc_dxGS_neg, 2))


df_pred <- transform(
  df_tract,
  dxGS_neg = predict(lunc_dxGS_neg, type = "response")
)
ggplot(
  data = df_pred,
  aes(
    x = nodeID,
    y = lgi_neg,
    fill = dxGS_neg,
    color = dxGS_neg,
    height = lgi_neg
  )
) +
  geom_tile() +
  facet_wrap(~dx_group, ncol = 2) +
  scale_fill_viridis("dti_fa") +
  scale_color_viridis("dti_fa") +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 50, 99)) +
  labs(x = "nodeID", y = "Negative LGI") +
  theme(legend.position = "right")

# L. Unc: GS intx neg LGI, decompose tract curve
lunc_dxGS_neg_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50) +
  s(lgi_neg) +
  te(nodeID, lgi_neg, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
  t2(
    nodeID, lgi_neg, dx_group,
    bs = c("cr", "tp", "re"),
    k = c(50, 10, 2),
    m = 2,
    full = TRUE
  ),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_dxGS_neg_decomp, rep = 1000)
summary(lunc_dxGS_neg_decomp)
plot(lunc_dxGS_neg_decomp)
plot_lunc_dxGS_neg_decomp <- getViz(lunc_dxGS_neg_decomp)
plot(sm(plot_lunc_dxGS_neg_decomp, 4))

# L. Unc: GS intx neu LGI
lunc_dxGS_neu <- bam(dti_fa ~ sex +
   s(subjectID, bs = "re") +
   te(nodeID, lgi_neu, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
   t2(
     nodeID, lgi_neu, dx_group,
     bs = c("cr", "tp", "re"),
     k = c(50, 10, 2),
     m = 2,
     full = TRUE
   ),
 data = df_tract,
 family = gaussian(),
 method = "fREML"
)
gam.check(lunc_dxGS_neu, rep = 1000)
summary(lunc_dxGS_neu)
plot(lunc_dxGS_neu)
plot_lunc_dxGS_neu <- getViz(lunc_dxGS_neu)
plot(sm(plot_lunc_dxGS_neu, 2))

df_pred <- transform(
  df_tract,
  dxGS_neu = predict(lunc_dxGS_neu, type = "response")
)
ggplot(
  data = df_pred,
  aes(
    x = nodeID,
    y = lgi_neg,
    fill = dxGS_neu,
    color = dxGS_neu,
    height = lgi_neu
  )
) +
  geom_tile() +
  facet_wrap(~dx_group, ncol = 2) +
  scale_fill_viridis("dti_fa") +
  scale_color_viridis("dti_fa") +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 50, 99)) +
  labs(x = "nodeID", y = "Neutral LGI") +
  theme(legend.position = "right")

# L. Unc: GS intx neu LGI, decompose tract curve
lunc_dxGS_neu_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50) +
  s(lgi_neu) +
  te(nodeID, lgi_neu, bs = c("cr", "tp"), k = c(50, 10), m = 2) +
  t2(
    nodeID, lgi_neu, dx_group,
    bs = c("cr", "tp", "re"),
    k = c(50, 10, 2),
    m = 2,
    full = TRUE
  ),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_dxGS_neu_decomp, rep = 1000)
summary(lunc_dxGS_neu_decomp)


# R. Unc ----
#
#

# subset df_afq, take complete cases, make factors
df_tract <- df_afq[which(df_afq$tractID == "UNC_R"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

# determine distribution
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete = F)

# build gam
runc_gaus <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(runc_gaus, rep = 1000)
summary(runc_gaus)
plot(runc_gaus)

# R. Unc: GS
runc_dxGS <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(runc_dxGS, rep = 1000)
compareML(runc_gaus, runc_dxGS) # runc_dxGS preferred
summary(runc_dxGS)
plot(runc_dxGS)

# R. Unc: GI
runc_dxGI <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(dx_group, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_group, bs = "cr", k = 50, m = 1),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(runc_dxGI, rep = 1000)
summary(runc_dxGI)
plot(runc_dxGI)
compareML(runc_dxGI, runc_dxGS) # runc_dxGS preferred

# R. Unc: GS difference via ordered factors
df_tract$dx_groupOF <- factor(df_tract$dx_group, ordered = T)
runc_dxGS_OF <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_groupOF, bs = "cr", k = 50, m = 2),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(runc_dxGS_OF, rep = 1000)
plot(runc_dxGS_OF)
summary(runc_dxGS_OF)
plot_runc_dxGS_OF <- getViz(runc_dxGS_OF)
plot(sm(plot_runc_dxGS_OF, 2))
plot(sm(plot_runc_dxGS_OF, 3)) +
  geom_hline(yintercept = 0)


# L. Cing ----
#
#

# subset df_afq, take complete cases, make factors
df_tract <- df_afq[which(df_afq$tractID == "CGC_L"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

# determine distribution
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete = F)
ggplot(df_tract, aes(y = dti_fa, x = nodeID)) +
  geom_point()

# build gam
lcgc_beta <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = betar(link = "logit"),
method = "fREML"
)
gam.check(lcgc_beta, rep = 1000)

lcgc_gamma <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
gam.check(lcgc_gamma, rep = 1000)
compareML(lcgc_beta, lcgc_gamma)

# L. Cing: GS
lcgc_dxGS <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
gam.check(lcgc_dxGS, rep = 1000)
compareML(lcgc_gamma, lcgc_dxGS) # almost equal
summary(lcgc_dxGS)
plot(lcgc_dxGS)

# L. Cing: GI
lcgc_dxGI <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(dx_group, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_group, bs = "cr", k = 50, m = 1),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
gam.check(lcgc_dxGI, rep = 1000)
summary(lcgc_dxGI)
plot(lcgc_dxGI)
compareML(lcgc_dxGI, lcgc_dxGS) # lcgc_dxGS preferred


# R. Cing ----
#
#

# subset df_afq, take complete cases, make factors
df_tract <- df_afq[which(df_afq$tractID == "CGC_R"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

# determine distribution
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete = F)
ggplot(df_tract, aes(y = dti_fa, x = nodeID)) +
  geom_point()

# build gam
rcgc_beta <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = betar(link = "logit"),
method = "fREML"
)
gam.check(rcgc_beta, rep = 1000)

rcgc_gamma <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
gam.check(rcgc_gamma, rep = 1000)
compareML(rcgc_beta, rcgc_gamma)

# R. Cing: GS
rcgc_dxGS <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, dx_group, bs = "fs", k = 50, m = 2),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
gam.check(rcgc_dxGS, rep = 1000)
compareML(rcgc_gamma, rcgc_dxGS) # almost equal
summary(rcgc_dxGS)
plot(rcgc_dxGS)

# R. Cing: GI
rcgc_dxGI <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(dx_group, bs = "re") +
  s(nodeID, bs = "cr", k = 50, m = 2) +
  s(nodeID, by = dx_group, bs = "cr", k = 50, m = 1),
data = df_tract,
family = Gamma(link = "logit"),
method = "fREML"
)
gam.check(rcgc_dxGI, rep = 1000)
summary(rcgc_dxGI)
plot(rcgc_dxGI)
compareML(rcgc_dxGI, rcgc_dxGS) # rcgc_dxGS preferred
