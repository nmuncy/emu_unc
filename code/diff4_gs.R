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

# subset df_afq, take complete cases, make factors
df_tract <- df_afq[which(df_afq$tractID == "UNC_L"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

# build gam
lunc_G <- bam(dti_fa ~ sex +
                     s(subjectID, bs = "re") +
                     s(nodeID, bs = "cr", k = 40),
                   data = df_tract,
                   family = gaussian(),
                   method = "fREML"
)
# gam.check(lunc_normal, rep = 1000)
# summary(lunc_normal)
plot_lunc_G <- getViz(lunc_G)
plot(sm(plot_lunc_G, 2))

# assess fit for e/subj
df_pred <- transform(
  df_tract, pred = predict(lunc_G, type = "response", se.fit = T)
)
ggplot(data=df_pred, aes(x=nodeID, y=dti_fa, group=subjectID)) +
  facet_wrap(~subjectID) +
  geom_ribbon(
    aes(ymin=pred.fit - 2*pred.se.fit, ymax=pred.fit + 2*pred.se.fit, x=nodeID),
    data=df_pred, 
    alpha=0.3, 
    inherit.aes=FALSE
  ) +
  geom_line(aes(y=pred.fit), data=df_pred) +
  geom_point()

# # subset data
# sub_list <- seq(9, 89, by = 10)
# df_sub <- df_tract[seq(10, dim(df_tract)[1], by = 10),]
# df_sub <- df_sub[-which(df_sub$nodeID == 99),]
# 
# lunc_GS <- bam(dti_fa ~ sex +
#                 s(nodeID, bs = "cr", k = 9, m = 2) +
#                 s(nodeID, subjectID, bs = "fs", k= 9, m = 2),
#               data = df_sub,
#               family = gaussian(),
#               method = "fREML",
#               discrete = T
# )
# gam.check(lunc_GS, rep = 1000)
# plot(lunc_GS)
# 
# lunc_GS_pred <- predict(lunc_GS, se.fit=T)
# df_sub <- transform(df_sub, 
#                     modGS = lunc_GS_pred$fit,
#                     modGS_se= lunc_GS_pred$se.fit)
# 
# ggplot(data=df_sub, aes(x=nodeID, y=dti_fa, group=subjectID)) +
#   facet_wrap(~subjectID) +
#   geom_ribbon(aes(ymin=modGS-2*modGS_se, 
#               ymax=modGS+2*modGS_se), alpha = 0.3) +
#   geom_line(aes(y=modGS)) +
#   geom_point()
# 
# 
# # full data
# lunc_GS <- bam(dti_fa ~ sex +
#                  s(nodeID, bs = "cr", k = 90, m = 2) +
#                  s(nodeID, subjectID, bs = "fs", k = 90, m = 2),
#                data = df_tract,
#                family = gaussian(),
#                method = "fREML",
#                discrete = T
# )
# gam.check(lunc_GS, rep = 1000)
# saveRDS(lunc_GS, file="/Users/nmuncy/Desktop/lunc_GS.rds")
# plot(lunc_GS)
# 
# lunc_GS_pred <- predict(lunc_GS, se.fit=T)
# df_tract <- transform(df_tract, 
#                     modGS = lunc_GS_pred$fit,
#                     modGS_se= lunc_GS_pred$se.fit)
# 
# ggplot(data=df_tract, aes(x=nodeID, y=dti_fa, group=subjectID)) +
#   facet_wrap(~subjectID) +
#   geom_ribbon(aes(ymin=modGS-2*modGS_se, 
#                   ymax=modGS+2*modGS_se), alpha = 0.3) +
#   geom_line(aes(y=modGS)) +
#   geom_point()


# Smooth for group
lunc_groupGS <- bam(dti_fa ~ sex +
                 s(dx, bs = "re") +
                 s(subjectID, bs = "re") +
                 s(nodeID, bs = "cr", k = 50, m = 2) +
                 s(nodeID, dx, bs = "fs", k = 50, m = 2),
               data = df_tract,
               family = gaussian(),
               method = "fREML"
)
gam.check(lunc_groupGS, rep = 1000)
summary(lunc_groupGS)
plot(lunc_groupGS)

plot_lunc_groupGS <- getViz(lunc_groupGS)
plot(sm(plot_lunc_groupGS, 1))
plot(sm(plot_lunc_groupGS, 2))
plot(sm(plot_lunc_groupGS, 3))
plot(sm(plot_lunc_groupGS, 4))

# TODO split dataframes by dx



# smooth for group * pscared
# ti(nodeID, pscared, by = dx, bs = c("cr", "tp"), k = c(50, 10), m = 1)
lunc_groupGS_pscared <- bam(dti_fa ~ sex +
    s(dx, bs = "re") +
    s(subjectID, bs = "re") +
    s(pscared, k = 10) +
    s(nodeID, bs = "cr", k = 50, m = 2) +
    s(nodeID, dx, bs = "fs", k = 50, m = 2) +
    ti(nodeID, pscared, bs = c("cr", "tp"), k = c(50, 10), m = 2),
  data = df_tract,
  family = gaussian(),
  method = "fREML"
)
gam.check(lunc_groupGS_pscared, rep = 1000)
summary(lunc_groupGS_pscared)
plot(lunc_groupGS_pscared)

plot_lunc_groupGS_pscared <- getViz(lunc_groupGS_pscared)
plot(sm(plot_lunc_groupGS_pscared, 1))
plot(sm(plot_lunc_groupGS_pscared, 2))
plot(sm(plot_lunc_groupGS_pscared, 3))
plot(sm(plot_lunc_groupGS_pscared, 4))
plot(sm(plot_lunc_groupGS_pscared, 5))
plot(sm(plot_lunc_groupGS_pscared, 6))
plot(sm(plot_lunc_groupGS_pscared, 7))
plot(sm(plot_lunc_groupGS_pscared, 8))
plot(sm(plot_lunc_groupGS_pscared, 9))


# ordered dx factors
# df_tract$dxOF <- factor(df_tract$dx, ordered = T)
lunc_groupGS_pscared_OF <- bam(dti_fa ~ sex + dx +
    s(subjectID, bs = "re") +
    s(nodeID, bs = "cr", k = 50, m = 2) +
    s(nodeID, dx, bs = "fs", k = 50, m = 2) +
    ti(nodeID, pscared, by = dxOF, bs = c("cr", "tp"), k = c(50, 3), m = 2) ,
  data = df_tract,
  family = gaussian(),
  method = "fREML"
)
gam.check(lunc_groupGS_pscared_OF, rep = 1000)
summary(lunc_groupGS_pscared_OF)
plot(lunc_groupGS_pscared_OF)

plot_lunc_groupGS_pscared_OF <- getViz(lunc_groupGS_pscared_OF)
plot(sm(plot_lunc_groupGS_pscared_OF, 1))
plot(sm(plot_lunc_groupGS_pscared_OF, 2))
plot(sm(plot_lunc_groupGS_pscared_OF, 3))
plot(sm(plot_lunc_groupGS_pscared_OF, 4))
plot(sm(plot_lunc_groupGS_pscared_OF, 5))
plot(sm(plot_lunc_groupGS_pscared_OF, 6))
plot(sm(plot_lunc_groupGS_pscared_OF, 7))
