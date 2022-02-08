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
df_afq$dx <- factor(df_afq$dx)
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
lunc_normal <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", by = sex, k = 40),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_normal, rep = 1000)
summary(lunc_normal)
plot_lunc_normal <- getViz(lunc_normal)
plot(sm(plot_lunc_normal, 1))
plot(sm(plot_lunc_normal, 2))
plot(sm(plot_lunc_normal, 3))

# get stat of male smooth diff from female
df_tract$sexOF <- factor(df_tract$sex, ordered = T)
lunc_normalOF <- bam(dti_fa ~ sexOF +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 40) +
  s(nodeID, by = sexOF, bs = "cr", k = 40),
data = df_tract,
family = gaussian(),
method = "fREML"
)
summary(lunc_normalOF)
plot(lunc_normalOF)
plot_lunc_normalOF <- getViz(lunc_normalOF)
plot(sm(plot_lunc_normalOF, 3)) + geom_hline(yintercept = 0)
gamtabs(lunc_normalOF)
report_stats(lunc_normalOF)

# add pds covariate
lunc_pds <- bam(dti_fa ~ sex +
  s(pds, by = sex) +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 40),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_pds, rep = 1000)
compareML(lunc_normal, lunc_pds) # fit normal preferred
summary(lunc_pds)

plot_lunc_pds <- getViz(lunc_pds)
plot(sm(plot_lunc_pds, 4))
gamtabs(lunc_pds)
report_stats(lunc_pds)


# L. Unc: GAM with continuous P-SCARED ----
#
# Rather than looking for group differences, see if adding anxiety measure
# changes model fit. Then plot those two splines to find differences.

# intx of scared with fa
lunc_pscared_intx <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  te(nodeID, pscared, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_pscared_intx, rep = 1000)
k.check(lunc_pscared_intx)

summary(lunc_pscared_intx)
compareML(lunc_normal, lunc_pscared_intx)

# topo plot
plot(lunc_pscared_intx)
plot_lunc_pscared_intx <- getViz(lunc_pscared_intx)
plot(sm(plot_lunc_pscared_intx, 2))

# 3d plot
vis.gam(
  lunc_pscared_intx,
  view = c("nodeID", "pscared"),
  color = "topo",
  theta = 330,
  phi = 40
)
vis.gam(lunc_pscared_intx, view = c("nodeID", "pscared"), theta = 60)


# decompose intx term
lunc_pscared_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 40) +
  s(pscared, k = 10) +
  ti(nodeID, pscared, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_pscared_decomp, rep = 1000)
k.check(lunc_pscared_decomp)

summary(lunc_pscared_decomp)
compareML(lunc_normal, lunc_pscared_decomp)

# topo plot
plot(lunc_pscared_decomp)
plot_lunc_pscared_decomp <- getViz(lunc_pscared_decomp)
plot(sm(plot_lunc_pscared_decomp, 2))
plot(sm(plot_lunc_pscared_decomp, 3))
plot(sm(plot_lunc_pscared_decomp, 4))

vis.gam(
  lunc_pscared_decomp,
  view = c("nodeID", "pscared"),
  color = "topo",
  phi = 40,
  theta = 330
)
# vis.gam(lunc_pscared_decomp, view = c("nodeID", "pscared"), plot.type = "contour", color = "topo", too.far = 0.1)
# vis.gam(lunc_pscared_decomp, view = c("nodeID", "pscared"), plot.type = "contour", color = "topo", too.far = 0.05)


# L. Unc: GAM with continuous PARS-6 ----
#
#

# intx
lunc_pars6_intx <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  te(nodeID, pars6, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_pars6_intx, rep = 1000)
k.check(lunc_pars6_intx)

summary(lunc_pars6_intx)
compareML(lunc_normal, lunc_pars6_intx)

plot_lunc_pars6_intx <- getViz(lunc_pars6_intx)
plot(sm(plot_lunc_pars6_intx, 2))

# decompose
lunc_pars6_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 40) +
  s(pars6) +
  ti(nodeID, pars6, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_pars6_decomp, rep = 1000)
k.check(lunc_pars6_decomp)

summary(lunc_pars6_decomp)
compareML(lunc_normal, lunc_pars6_decomp)

plot_lunc_pars6_decomp <- getViz(lunc_pars6_decomp)
plot(sm(plot_lunc_pars6_decomp, 2))
plot(sm(plot_lunc_pars6_decomp, 3))
plot(sm(plot_lunc_pars6_decomp, 4))

vis.gam(
  lunc_pars6_decomp,
  view = c("nodeID", "pars6"),
  color = "topo",
  phi = 40,
  theta = 330
)


# L. Unc: GAM with continuous C-SCARED ----
#

# intx of scared with fa
lunc_cscared_intx <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  te(nodeID, cscared, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_cscared_intx, rep = 1000)
k.check(lunc_cscared_intx)

summary(lunc_cscared_intx)
compareML(lunc_normal, lunc_cscared_intx)

# topo plot
plot(lunc_cscared_intx)
plot_lunc_cscared_intx <- getViz(lunc_cscared_intx)
plot(sm(plot_lunc_cscared_intx, 2))

# 3d plot
vis.gam(lunc_cscared_intx, view = c("nodeID", "cscared"), color = "topo", theta = 330, phi = 40)
vis.gam(lunc_cscared_intx, view = c("nodeID", "cscared"), theta = 60)

# parent's scared - decompose intx term
lunc_cscared_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 40) +
  s(cscared, k = 10) +
  ti(nodeID, cscared, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(lunc_cscared_decomp, rep = 1000)
k.check(lunc_cscared_decomp)

summary(lunc_cscared_decomp)
compareML(lunc_normal, lunc_cscared_decomp)

# topo plot
plot(lunc_cscared_decomp)
plot_lunc_cscared_decomp <- getViz(lunc_cscared_decomp)
plot(sm(plot_lunc_cscared_decomp, 2))
plot(sm(plot_lunc_cscared_decomp, 3))
plot(sm(plot_lunc_cscared_decomp, 4))

vis.gam(
  lunc_cscared_decomp,
  view = c("nodeID", "cscared"),
  color = "topo",
  phi = 40,
  theta = 330
)
# vis.gam(lunc_cscared_decomp, view = c("nodeID", "cscared"), plot.type = "contour", color = "topo", too.far = 0.1)
# vis.gam(lunc_cscared_decomp, view = c("nodeID", "cscared"), plot.type = "contour", color = "topo", too.far = 0.05)


# L. Unc: GAM with summed C/P-SCARED ----
#

df_tract$sum_scared <- df_tract$cscared + df_tract$pscared

# intx of scared with fa
lunc_sumscared_intx <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, sum_scared, bs = c("cr", "tp"), k = c(40, 10)),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
)
gam.check(lunc_sumscared_intx, rep = 1000)
summary(lunc_sumscared_intx)
plot_lunc_sumscared_intx <- getViz(lunc_sumscared_intx)
plot(sm(plot_lunc_sumscared_intx, 2))

# parent's scared - decompose intx term
lunc_sumscared_decomp <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 40) +
    s(sum_scared, k = 10) +
    ti(nodeID, sum_scared, bs = c("cr", "tp"), k = c(40, 10)),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
)
summary(lunc_sumscared_decomp)

# topo plot
plot(lunc_sumscared_decomp)
plot_lunc_sumscared_decomp <- getViz(lunc_sumscared_decomp)
plot(sm(plot_lunc_sumscared_decomp, 2))
plot(sm(plot_lunc_sumscared_decomp, 3))
plot(sm(plot_lunc_sumscared_decomp, 4))



# L. Unc: GAM with continuous P-SCARED by dx ----
lunc_dx_hgam <- bam(dti_fa ~ sex +
    s(subjectID, bs = "re") +
    te(nodeID, pscared, bs = c("cr", "tp"), k = c(40, 10), m = 2) +
    t2(
      nodeID, pscared, dx,
      bs = c("cr", "tp", "re"),
      k = c(40, 10, 3),
      m = 2,
      full = T
    ),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
)
summary(lunc_dx_hgam)
plot(lunc_dx_hgam)


lunc_dx_intx <- bam(dti_fa ~ sex + dx +
   s(subjectID, bs = "re") +
   te(nodeID, pscared, by = dx, bs = c("cr", "tp"), k = c(40, 10)),
 data = df_tract,
 family = gaussian(),
 method = "fREML"
)
gam.check(lunc_dx_intx)
summary(lunc_dx_intx)

plot(lunc_dx_intx)
plot_lunc_dx_intx <- getViz(lunc_dx_intx)
plot(sm(plot_lunc_dx_intx, 2))
plot(sm(plot_lunc_dx_intx, 3))
plot(sm(plot_lunc_dx_intx, 4))


lunc_dx_decomp <- bam(dti_fa ~ sex + dx +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 40) +
    s(pscared, k = 10) +
    ti(nodeID, pscared, by = dx, bs = c("cr", "tp"), k = c(40, 10)),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
)
summary(lunc_dx_decomp)
plot(lunc_dx_decomp)
plot_lunc_dx_decomp <- getViz(lunc_dx_decomp)
plot(sm(plot_lunc_dx_decomp, 2))
plot(sm(plot_lunc_dx_decomp, 3))
plot(sm(plot_lunc_dx_decomp, 4))
plot(sm(plot_lunc_dx_decomp, 5))
plot(sm(plot_lunc_dx_decomp, 6))




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
af_pscared_intx <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  te(nodeID, pscared, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(af_pscared_intx, rep = 1000)

summary(af_pscared_intx)
compareML(af_normal, af_pscared_intx)

# topo plot
plot(af_pscared_intx)
plot_af_pscared_intx <- getViz(af_pscared_intx)
plot(sm(plot_af_pscared_intx, 2))

# decompose intx term
af_pscared_decomp <- bam(dti_fa ~ sex +
  s(subjectID, bs = "re") +
  s(nodeID, bs = "cr", k = 40) +
  s(pscared, k = 10) +
  ti(nodeID, pscared, bs = c("cr", "tp"), k = c(40, 10)),
data = df_tract,
family = gaussian(),
method = "fREML"
)
gam.check(af_pscared_decomp, rep = 1000)
summary(af_pscared_decomp)
compareML(af_normal, af_pscared_decomp)

# topo plot
plot(af_pscared_decomp)
plot_af_pscared_decomp <- getViz(af_pscared_decomp)
plot(sm(plot_af_pscared_decomp, 2))
plot(sm(plot_af_pscared_decomp, 3))
plot(sm(plot_af_pscared_decomp, 4))

vis.gam(
  af_pscared_decomp,
  view = c("nodeID", "pscared"),
  color = "topo",
  phi = 40,
  theta = 330
)

# intx by dx
af_dx_intx <- bam(dti_fa ~ sex + dx +
     s(subjectID, bs = "re") +
     te(nodeID, pscared, by = dx, bs = c("cr", "tp"), k = c(40, 10)),
   data = df_tract,
   family = gaussian(),
   method = "fREML"
)
plot(af_dx_intx)
plot_af_dx_intx <- getViz(af_dx_intx)
plot(sm(plot_af_dx_intx, 2))
plot(sm(plot_af_dx_intx, 3))
plot(sm(plot_af_dx_intx, 4))

af_dx_decomp <- bam(dti_fa ~ sex + dx +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 40) +
    s(pscared, k = 10) +
    ti(nodeID, pscared, by = dx, bs = c("cr", "tp"), k = c(40, 10)),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
)
summary(af_dx_decomp)
plot(af_dx_decomp)
plot_af_dx_decomp <- getViz(af_dx_decomp)
plot(sm(plot_af_dx_decomp, 2))
plot(sm(plot_af_dx_decomp, 3))
plot(sm(plot_af_dx_decomp, 4))
plot(sm(plot_af_dx_decomp, 5))
plot(sm(plot_af_dx_decomp, 6))



# P. Forceps: Determine data distribution, base GAM function ----
# subset df_afq, take complete cases, make factors
df_tract <- df_afq[which(df_afq$tractID == "CST_L"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

# determine distribution
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete = F)

# build gam
fit_gamma <- bam(dti_fa ~
                    s(subjectID, bs = "re") +
                    s(nodeID, bs = "cr", k = 40),
                  data = df_tract,
                  family = Gamma(link = "logit"),
                  method = "fREML"
)
gam.check(fit_gamma, rep = 1000)

fit_beta <- bam(dti_fa ~
                   s(subjectID, bs = "re") +
                   s(nodeID, bs = "cr", k = 40),
                 data = df_tract,
                 family = betar(link = "logit"),
                 method = "fREML"
)
gam.check(fit_beta, rep = 1000)
compareML(fit_beta, fit_gamma)

# add sex
pf_beta <- bam(dti_fa ~ sex +
                    s(subjectID, bs = "re") +
                    s(nodeID, bs = "cr", k = 40),
                  data = df_tract,
                  family = betar(link = "logit"),
                  method = "fREML"
)
plot(pf_beta)
summary(pf_beta)


# P. Forceps: GAM with continuous Parent's SCARED ----

# intx of scared with fa
pf_pscared_intx <- bam(dti_fa ~ sex +
                         s(subjectID, bs = "re") +
                         te(nodeID, pscared, bs = c("cr", "tp"), k = c(40, 10)),
                       data = df_tract,
                       family = betar(link = "logit"),
                       method = "fREML"
)
gam.check(pf_pscared_intx, rep = 1000)
summary(pf_pscared_intx)
compareML(pf_beta, pf_pscared_intx)

plot(pf_pscared_intx)
plot_pf_pscared_intx <- getViz(pf_pscared_intx)
plot(sm(plot_pf_pscared_intx, 2))

# decompose intx term
pf_pscared_decomp <- bam(dti_fa ~ sex +
                           s(subjectID, bs = "re") +
                           s(nodeID, bs = "cr", k = 40) +
                           s(pscared, k = 10) +
                           ti(nodeID, pscared, bs = c("cr", "tp"), k = c(40, 10)),
                         data = df_tract,
                         family = gaussian(),
                         method = "fREML"
)
gam.check(pf_pscared_decomp, rep = 1000)
summary(pf_pscared_decomp)

# topo plot
plot(pf_pscared_decomp)
plot_pf_pscared_decomp <- getViz(pf_pscared_decomp)
plot(sm(plot_pf_pscared_decomp, 2))
plot(sm(plot_pf_pscared_decomp, 3))
plot(sm(plot_pf_pscared_decomp, 4))


# intx by dx
pf_dx_intx <- bam(dti_fa ~ sex + dx +
                    s(subjectID, bs = "re") +
                    te(nodeID, pscared, by = dx, bs = c("cr", "tp"), k = c(40, 10)),
                  data = df_tract,
                  family = gaussian(),
                  method = "fREML"
)
plot(pf_dx_intx)
plot_pf_dx_intx <- getViz(pf_dx_intx)
plot(sm(plot_pf_dx_intx, 2))
plot(sm(plot_pf_dx_intx, 3))
plot(sm(plot_pf_dx_intx, 4))

pf_dx_decomp <- bam(dti_fa ~ sex + dx +
                      s(subjectID, bs = "re") +
                      s(nodeID, bs = "cr", k = 40) +
                      s(pscared, k = 10) +
                      ti(nodeID, pscared, by = dx, bs = c("cr", "tp"), k = c(40, 10)),
                    data = df_tract,
                    family = gaussian(),
                    method = "fREML"
)
summary(pf_dx_decomp)
plot(pf_dx_decomp)
plot_pf_dx_decomp <- getViz(pf_dx_decomp)
plot(sm(plot_pf_dx_decomp, 2))
plot(sm(plot_pf_dx_decomp, 3))
plot(sm(plot_pf_dx_decomp, 4))
plot(sm(plot_pf_dx_decomp, 5))
plot(sm(plot_pf_dx_decomp, 6))



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
