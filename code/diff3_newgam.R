library("mgcViz")

df_afq <- read.csv("/Users/nmuncy/Desktop/AFQ_dataframe.csv")

tract <- "UNC_L"
df_tract <- df_afq[which(df_afq$tractID == tract), ]
df_tract$sex <- factor(df_tract$sex)


fit_cov_normal <- bam(dti_fa ~ sex + sex * pds +
    s(subjectID, bs = "re") +
    s(nodeID, k = 40),
  data = df_tract,
  family = gaussian(link = "logit"),
  method = "REML"
)

fit_pscared <- bam(dti_fa ~ sex + sex * pds +
   s(subjectID, bs = "re") +
   te(nodeID, pscared, bs = c("cr", "tp"), k = c(40, 10)),
 data = df_tract,
 family = gaussian(),
 method = "REML"
)
gam.check(fit_pscared, rep = 1000)
fit_pscared <- getViz(fit_pscared)
plot(sm(fit_pscared, 2)) + l_fitRaster() + l_fitContour()

fit_pars <- bam(dti_fa ~ sex + sex * pds +
   s(subjectID, bs = "re") +
   te(nodeID, pars6, bs = c("cr", "tp"), k = c(50, 10)),
 data = df_tract,
 family = gaussian(),
 method = "REML"
)
gam.check(fit_pars, rep = 1000)
fit_pars <- getViz(fit_pars)
plot(sm(fit_pars, 2)) + l_fitRaster() + l_fitContour()
