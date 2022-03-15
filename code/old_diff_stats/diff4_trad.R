library("tools")
library("ggplot2")

data_dir <- file_path_as_absolute(paste0(getwd(), "/../data"))
# capture.output(sessionInfo(), file = paste0(data_dir, "R_session_info.txt"))
df_afq <- read.csv(paste0(data_dir, "/AFQ_dataframe.csv"))
df_afq$sex <- factor(df_afq$sex)
df_afq$dx <- factor(df_afq$dx)
df_afq$subjectID <- factor(df_afq$subjectID)

df_tract <- df_afq[which(df_afq$tractID == "UNC_L"), ]
df_tract <- df_tract[complete.cases(df_tract), ]

subj_list <- unique(df_tract$subjectID)
df_tract$avgFA <- NA
for (subj in subj_list){
  ind_subj <- which(df_tract$subjectID == subj)
  df_tract[ind_subj, ]$avgFA <- mean(df_tract[ind_subj, ]$dti_fa)
}
df_trad <- df_tract[which(df_tract$nodeID == 0),]

fit_lm <- lm(avgFA ~ pscared * dx, data = df_trad)
summary(fit_lm)
ggplot(data = df_trad, aes(x = pscared, y = avgFA, color = dx)) +
  geom_ribbon(
    stat='smooth', method = "lm", se=TRUE, alpha=0.1, 
    aes(color = NULL, group = factor(dx))) +
  geom_line(stat='smooth', method = "lm", alpha=1) +
  geom_point()

