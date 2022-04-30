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

ggplot(data = df_tract, aes(x = nodeID, y = dti_fa, color = dx_group)) +
  geom_point()

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
h_plot <- getViz(lunc_gaus)
plot(sm(h_plot, 2)) + ggtitle("L. Unc G")

# bad k
lunc_gaus <- bam(dti_fa ~ sex +
                   s(subjectID, bs = "re") +
                   s(nodeID, bs = "cr", k = 100),
                 data = df_tract,
                 family = gaussian(),
                 method = "fREML"
)
plot(lunc_gaus)
h_plot <- getViz(lunc_gaus)
plot(sm(h_plot, 2)) + ggtitle("L. Unc G")

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
h_plot <- getViz(lunc_dxGS)
plot(sm(h_plot, 2)) + ggtitle("L. Unc GS")
plot(sm(h_plot, 3)) + ggtitle("L. Unc GS")

# L. Unc: GI
# 
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
h_plot <- getViz(lunc_dxGI)
plot(sm(h_plot, 3)) + ggtitle("L. Unc GI")
plot(sm(h_plot, 4)) + ggtitle("L. Unc GI")

# L. Unc: GS difference via ordered factors
df_tract$dx_groupOF <- factor(df_tract$dx_group, ordered = T)
lunc_dxGS_OF <- bam(dti_fa ~ sex + dx_groupOF +
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
plot(sm(plot_lunc_dxGS_OF, 3)) +
  geom_hline(yintercept = 0)
