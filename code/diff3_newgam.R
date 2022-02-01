library("mgcViz")
library("gratia")
library("itsadug")

df_afq <- read.csv("/Users/nmuncy/Desktop/AFQ_dataframe.csv")

tract <- "UNC_L"
df_tract <- df_afq[which(df_afq$tractID == tract), ]
df_tract$sex <- factor(df_tract$sex)
df_tract$subjectID <- factor(df_tract$subjectID, ordered = F)
df_tract <- na.omit(df_tract)



fit_covG <- bam(dti_fa ~ sex + sex * pds +
   s(nodeID, bs = "cr", k = 40) +
   s(subjectID, bs = "re"),
  data = df_tract,
  family = gaussian(link = "logit"),
  method = "REML"
)
gam.check(fit_covG, rep = 1000)
# draw(fit_cov)
fit_norm <- getViz(fit_covG)
plot(sm(fit_norm, 1))

pred_covG <- predict(
  fit_covG, 
  exclude_terms = c("pds", "sex"),
  values = list(pds = NULL, sex = NULL),
  se.fit=T, 
  type="response")

df_plot <- data.frame(nodeID = df_tract$nodeID, subjectID = df_tract$subjectID)
df_plot$fit <- pred_covG$fit
df_plot$se.fit <- pred_covG$se.fit
ggplot(data = df_tract, aes(x=nodeID, y=dti_fa, group=subjectID)) +
  facet_wrap(~subjectID) +
  geom_ribbon(aes(ymin=(fit-2*se.fit), ymax=(fit+2*se.fit), x=nodeID),
              data=df_plot, alpha=0.3, inherit.aes=F) +
  geom_line(aes(y=fit), data=df_plot) +
  geom_point()



fit_covGS <- bam(dti_fa ~ sex + sex * pds +
   s(nodeID, k = 40, m = 2) +
   s(nodeID, subjectID, k = 40, bs = "fs", m = 2),
 data = df_tract,
 family = gaussian(link = "logit"),
 method = "REML"
)
gam.check(fit_covGS, rep = 1000)
fit_GS <- getViz(fit_covGS)
plot(sm(fit_GS, 1))
compareML(fit_covG, fit_covGS)

pred_covGS <- predict(
  fit_covGS, 
  exclude_terms = c("pds", "sex"),
  values = list(pds = NULL, sex = NULL),
  se.fit=T, 
  type="response")

df_plot <- data.frame(nodeID = df_tract$nodeID, subjectID = df_tract$subjectID)
df_plot$fit <- pred_covGS$fit
df_plot$se.fit <- pred_covGS$se.fit
df_plot$dti_fa <- df_tract$dti_fa

ggplot(data = df_plot, aes(x=nodeID, y=dti_fa, group=subjectID)) +
  facet_wrap(~subjectID) +
  geom_ribbon(aes(ymin=(fit-2*se.fit), ymax=(fit+2*se.fit)),
              alpha=0.3) +
  geom_line(aes(y=fit)) +
  geom_point()

CO2_modGS_pred <- predict(CO2_modGS, se.fit=TRUE)
CO2 <- transform(CO2, 
                 modGS = CO2_modGS_pred$fit, 
                 modGS_se = CO2_modGS_pred$se.fit)

ggplot(data=CO2, aes(x=conc, y=uptake, group=Plant_uo)) +
  facet_wrap(~Plant_uo) +
  geom_ribbon(aes(ymin=exp(modGS-2*modGS_se),
                  ymax=exp(modGS+2*modGS_se)), alpha=0.25) +
  geom_line(aes(y=exp(modGS))) +
  geom_point()



fit_pscared <- bam(dti_fa ~ sex + sex * pds +
   s(subjectID, bs = "re") +
   t2(nodeID, pscared, bs = c("cr", "tp"), k = c(15, 10), m = 2, full = T),
 data = df_tract,
 family = gaussian(),
 method = "REML"
)

draw(fit_pscared)

gam.check(fit_pscared, rep = 1000)
k.check(fit_pscared)
summary(fit_pscared)
fit_pscared <- getViz(fit_pscared)
plot(sm(fit_pscared, 2)) + l_fitRaster() + l_fitContour()

compareML(fit_cov_normal, fit_pscared)




fit_pars <- bam(dti_fa ~ sex + sex * pds +
   s(subjectID, bs = "re") +
   t2(nodeID, pars6, bs = c("cr", "tp"), k = c(50, 10)),
 data = df_tract,
 family = gaussian(),
 method = "REML"
)
gam.check(fit_pars, rep = 1000)
fit_pars <- getViz(fit_pars)
plot(sm(fit_pars, 2)) + l_fitRaster() + l_fitContour()
