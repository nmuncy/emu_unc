library("ggplot2")
library("fitdistrplus")
library("mgcv")
library("itsadug")
library("tidymv")
library("dplyr")
library("mgcViz")
library("tools")


# Set Up ----
data_dir <- file_path_as_absolute(paste0(getwd(), "/../data"))
# capture.output(sessionInfo(), file = paste0(data_dir, "R_session_info.txt"))
df_afq <- read.csv(paste0(data_dir, "/AFQ_dataframe.csv"))
df_afq$sex <- factor(df_afq$sex)
df_afq$subjectID <- factor(df_afq$subjectID)


# L. Unc: Determine data distribution, base GAM function ----
#
#

# subset df_afq, take complete cases, make factors
df_tract <- df_afq[which(df_afq$tractID == "UNC_L"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

# determine distribution
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete = F)

# build gam
fit_normal <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", by = sex, k = 40),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(fit_normal, rep = 1000)
summary(fit_normal)
plot_normal <- getViz(fit_normal)
plot(sm(plot_normal, 1))
plot(sm(plot_normal, 2))
plot(sm(plot_normal, 3))

# get stat of male smooth diff from female
df_tract$sexOF <- factor(df_tract$sex, ordered = T)
fit_normalOF <- bam(dti_fa ~ sexOF +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 40) +
  s(nodeID, by = sexOF, bs = "cr", k = 40),
data = df_tract,
family = gaussian(),
method = "fREML"
)
summary(fit_normalOF)
plot(fit_normalOF)
plot_normalOF <- getViz(fit_normalOF)
plot(sm(plot_normalOF, 3)) + geom_hline(yintercept = 0)
gamtabs(fit_normalOF)
report_stats(fit_normalOF)

# add pds covariate
fit_pds <- bam(dti_fa ~ sex +
  s(pds, by = sex) +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 40),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(fit_pds, rep = 1000)
compareML(fit_normal, fit_pds) # fit normal preferred
summary(fit_pds)

plot_pds <- getViz(fit_pds)
plot(sm(plot_pds, 4))
gamtabs(fit_pds)
report_stats(fit_pds)


# L. Unc: GAM with continuous Parent's SCARED ----
#
# Rather than looking for group differences, see if adding anxiety measure
# changes model fit. Then plot those two splines to find differences.

# intx of scared with fa
fit_pscared_intx <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  te(nodeID, pscared, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(fit_pscared_intx, rep = 1000)
k.check(fit_pscared_intx)

summary(fit_pscared_intx)
compareML(fit_normal, fit_pscared_intx)

# topo plot
plot(fit_pscared_intx)
plot_pscared_intx <- getViz(fit_pscared_intx)
plot(sm(plot_pscared_intx, 2))

# 3d plot
vis.gam(
  fit_pscared_intx,
  view = c("nodeID", "pscared"),
  color = "topo",
  theta = 330,
  phi = 40
)
vis.gam(fit_pscared_intx, view = c("nodeID", "pscared"), theta = 60)


# decompose intx term
fit_pscared_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 40) +
  s(pscared, k = 10) +
  ti(nodeID, pscared, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(fit_pscared_decomp, rep = 1000)
k.check(fit_pscared_decomp)

summary(fit_pscared_decomp)
compareML(fit_normal, fit_pscared_decomp)

# topo plot
plot(fit_pscared_decomp)
plot_pscared_decomp <- getViz(fit_pscared_decomp)
plot(sm(plot_pscared_decomp, 2))
plot(sm(plot_pscared_decomp, 3))
plot(sm(plot_pscared_decomp, 4))

vis.gam(
  fit_pscared_decomp,
  view = c("nodeID", "pscared"),
  color = "topo",
  phi = 40,
  theta = 330
)
# vis.gam(fit_pscared_decomp, view = c("nodeID", "pscared"), plot.type = "contour", color = "topo", too.far = 0.1)
# vis.gam(fit_pscared_decomp, view = c("nodeID", "pscared"), plot.type = "contour", color = "topo", too.far = 0.05)


# L. Unc: GAM with continuous PARS-6 ----
#
#

# intx
fit_pars6_intx <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  te(nodeID, pars6, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(fit_pars6_intx, rep = 1000)
k.check(fit_pars6_intx)

summary(fit_pars6_intx)
compareML(fit_normal, fit_pars6_intx)

plot_pars6_intx <- getViz(fit_pars6_intx)
plot(sm(plot_pars6_intx, 2))

# decompose
fit_pars6_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 40) +
  s(pars6) +
  ti(nodeID, pars6, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(fit_pars6_decomp, rep = 1000)
k.check(fit_pars6_decomp)

summary(fit_pars6_decomp)
compareML(fit_normal, fit_pars6_decomp)

plot_pars6_decomp <- getViz(fit_pars6_decomp)
plot(sm(plot_pars6_decomp, 2))
plot(sm(plot_pars6_decomp, 3))
plot(sm(plot_pars6_decomp, 4))

vis.gam(
  fit_pars6_decomp,
  view = c("nodeID", "pars6"),
  color = "topo",
  phi = 40,
  theta = 330
)


# L. Unc: GAM with continuous Child's SCARED ----
#

# intx of scared with fa
fit_cscared_intx <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  te(nodeID, cscared, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(fit_cscared_intx, rep = 1000)
k.check(fit_cscared_intx)

summary(fit_cscared_intx)
compareML(fit_normal, fit_cscared_intx)

# topo plot
plot(fit_cscared_intx)
plot_cscared_intx <- getViz(fit_cscared_intx)
plot(sm(plot_cscared_intx, 2))

# 3d plot
vis.gam(fit_cscared_intx, view = c("nodeID", "cscared"), color = "topo", theta = 330, phi = 40)
vis.gam(fit_cscared_intx, view = c("nodeID", "cscared"), theta = 60)

# parent's scared - decompose intx term
fit_cscared_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 40) +
  s(cscared, k = 10) +
  ti(nodeID, cscared, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(fit_cscared_decomp, rep = 1000)
k.check(fit_cscared_decomp)

summary(fit_cscared_decomp)
compareML(fit_normal, fit_cscared_decomp)

# topo plot
plot(fit_cscared_decomp)
plot_cscared_decomp <- getViz(fit_cscared_decomp)
plot(sm(plot_cscared_decomp, 2))
plot(sm(plot_cscared_decomp, 3))
plot(sm(plot_cscared_decomp, 4))

vis.gam(
  fit_cscared_decomp,
  view = c("nodeID", "cscared"),
  color = "topo",
  phi = 40,
  theta = 330
)
# vis.gam(fit_cscared_decomp, view = c("nodeID", "cscared"), plot.type = "contour", color = "topo", too.far = 0.1)
# vis.gam(fit_cscared_decomp, view = c("nodeID", "cscared"), plot.type = "contour", color = "topo", too.far = 0.05)


# L. Unc: GAM with summed C/P-SCARED ----
#

df_tract$sum_scared <- df_tract$cscared + df_tract$pscared

# intx of scared with fa
fit_sumscared_intx <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, sum_scared, bs = c("cr", "tp"), k = c(40, 10)),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
)
gam.check(fit_sumscared_intx, rep = 1000)
summary(fit_sumscared_intx)
plot_sumscared_intx <- getViz(fit_sumscared_intx)
plot(sm(plot_sumscared_intx, 2))

# parent's scared - decompose intx term
fit_sumscared_decomp <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 40) +
    s(sum_scared, k = 10) +
    ti(nodeID, sum_scared, bs = c("cr", "tp"), k = c(40, 10)),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
)
summary(fit_sumscared_decomp)

# topo plot
plot(fit_sumscared_decomp)
plot_sumscared_decomp <- getViz(fit_sumscared_decomp)
plot(sm(plot_sumscared_decomp, 2))
plot(sm(plot_sumscared_decomp, 3))
plot(sm(plot_sumscared_decomp, 4))



# A. Forceps: Determine data distribution, base GAM function ----
# subset df_afq, take complete cases, make factors
df_tract <- df_afq[which(df_afq$tractID == "FA"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

# determine distribution
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete = F)

# build gam
fit_normal <- bam(dti_fa ~
s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 40),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(fit_normal, rep = 1000)

# fit_gamma <- bam(dti_fa ~
#                     s(subjectID, bs = "re") +
#                     s(nodeID, bs = "cr", k = 40),
#                   data = df_tract,
#                   family = Gamma(link = "logit"),
#                   method = "fREML"
# )
# gam.check(fit_gamma, rep = 1000)
# compareML(fit_gamma, fit_normal)
#
# fit_beta <- bam(dti_fa ~
#                    s(subjectID, bs = "re") +
#                    s(nodeID, bs = "cr", k = 40),
#                  data = df_tract,
#                  family = betar(link = "logit"),
#                  method = "fREML"
# )
# gam.check(fit_beta, rep = 1000)
# compareML(fit_beta, fit_normal)

# add sex
fit_normal <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 40),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(fit_normal, rep = 1000)
plot(fit_normal)

summary(fit_normal)
plot_normal <- getViz(fit_normal)
plot(sm(plot_normal, 1))
plot(sm(plot_normal, 2))


# A. Forceps: GAM with continuous Parent's SCARED ----

# intx of scared with fa
fit_pscared_intx <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  te(nodeID, pscared, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(fit_pscared_intx, rep = 1000)

summary(fit_pscared_intx)
compareML(fit_normal, fit_pscared_intx)

# topo plot
plot(fit_pscared_intx)
plot_pscared_intx <- getViz(fit_pscared_intx)
plot(sm(plot_pscared_intx, 2))

# decompose intx term
fit_pscared_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 40) +
  s(pscared, k = 10) +
  ti(nodeID, pscared, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(fit_pscared_decomp, rep = 1000)
summary(fit_pscared_decomp)
compareML(fit_normal, fit_pscared_decomp)

# topo plot
plot(fit_pscared_decomp)
plot_pscared_decomp <- getViz(fit_pscared_decomp)
plot(sm(plot_pscared_decomp, 2))
plot(sm(plot_pscared_decomp, 3))
plot(sm(plot_pscared_decomp, 4))

vis.gam(
  fit_pscared_decomp,
  view = c("nodeID", "pscared"),
  color = "topo",
  phi = 40,
  theta = 330
)


# L. Cing: Determine data distribution, base GAM function ----
#
#

# subset df_afq, take complete cases, make factors
df_tract <- df_afq[which(df_afq$tractID == "CGC_L"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

# determine distribution
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete = F)

hist(log10(max(df_tract$dti_fa + 1) - df_tract$dti_fa))
descdist(log10(max(df_tract$dti_fa + 1) - df_tract$dti_fa), discrete = F)


# # build gam
# fit_normal <- bam(dti_fa ~
#                     s(subjectID, bs = "re") +
#                     s(nodeID, bs = "cr", k = 40),
#                   data = df_tract,
#                   family = gaussian(),
#                   method = "fREML"
# )
# gam.check(fit_normal, rep = 1000)
#
# fit_beta <- bam(dti_fa ~
#                     s(subjectID, bs = "re") +
#                     s(nodeID, bs = "cr", k = 40),
#                   data = df_tract,
#                   family = betar(link = "logit"),
#                   method = "fREML"
# )
# gam.check(fit_beta, rep = 1000)
# compareML(fit_normal, fit_beta) # normal preferred
#
# fit_gamma <- bam(dti_fa ~
#                   s(subjectID, bs = "re") +
#                   s(nodeID, bs = "cr", k = 40),
#                 data = df_tract,
#                 family = Gamma(link = "logit"),
#                 method = "fREML"
# )
# gam.check(fit_gamma, rep = 1000)
# compareML(fit_normal, fit_gamma) # normal preferred
#
# fit_normal <- bam(dti_fa ~ sex +
#                     s(subjectID, bs = "re") +
#                     s(nodeID, by= sex, bs = "cr", k = 40),
#                   data = df_tract,
#                   family = gaussian(),
#                   method = "fREML"
# )
# gam.check(fit_normal, rep = 1000)
# plot(fit_normal)
#
#
# summary(fit_normal)
# plot_normal <- getViz(fit_normal)
# plot(sm(plot_normal, 1))
# plot(sm(plot_normal, 2))
# plot(sm(plot_normal, 3))
#
# # get stat of male smooth diff from female
# df_tract$sexOF <- factor(df_tract$sex, ordered = T)
# fit_normalOF <- bam(dti_fa ~ sexOF +
#                       s(subjectID, bs = "re") +
#                       s(nodeID, bs = "cr", k = 40) +
#                       s(nodeID, by = sexOF, bs = "cr", k = 40),
#                     data = df_tract,
#                     family = gaussian(),
#                     method = "fREML"
# )
# summary(fit_normalOF)
# plot(fit_normalOF)
# plot_normalOF <- getViz(fit_normalOF)
# plot(sm(plot_normalOF, 3)) + geom_hline(yintercept = 0)
# gamtabs(fit_normalOF)
# report_stats(fit_normalOF)
#
# # add pds covariate
# fit_pds <- bam(dti_fa ~ sex +
#                  s(pds, by = sex) +
#                  s(subjectID, bs = "re") +
#                  s(nodeID, bs = "cr", k = 40),
#                data = df_tract,
#                family = gaussian(),
#                method = "fREML"
# )
# gam.check(fit_pds, rep = 1000)
# compareML(fit_normal, fit_pds) # fit normal preferred
# summary(fit_pds)
#
# plot_pds <- getViz(fit_pds)
# plot(sm(plot_pds, 4))
# gamtabs(fit_pds)
# report_stats(fit_pds)
