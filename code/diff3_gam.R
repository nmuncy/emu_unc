library("ggplot2")
library("fitdistrplus")
library("mgcv")
library("itsadug")


one_dir <- "/Users/nmuncy/Florida International University/EMU Study - Documents/EMU Data/current_working_data/datasets/_full_working_dataset/"
proj_dir <- "/Volumes/homes/MaDLab/projects/McMakin_EMUR01/"
data_dir <- paste0(proj_dir, "derivatives/emu_unc/reco_afq/")
out_dir <- paste0(proj_dir, "derivatives/emu_unc/analyses/")

# capture.output(sessionInfo(), file = paste0(data_dir, "R_session_info.txt"))

# get dataframes
df_summary <- read.csv(paste0(one_dir, "emuR01_summary_latest.csv"))
df_afq <- read.csv(paste0(data_dir, "tract_profiles.csv"))

# start age, sex, pds, parent's scared columns
df_afq$age <- df_afq$sex <- df_afq$pds <- 
  df_afq$pars6 <- df_afq$pars6_group <-
  df_afq$pscared <- df_afq$pscared_group <- NA

# get list of subjects
subj_list <- unique(df_afq$subjectID)

# fill data for each subject
for(subj in subj_list){
  
  # determine indices of subjects in dfs
  ind_afq <- which(df_afq$subjectID == subj)
  ind_summ <- which(df_summary$emu_study_id == subj)
  if(length(ind_summ) == 0){
    print(paste("Missing summary data for", subj))
    next
  }
  
  # fill age, sex, pds
  df_afq[ind_afq,]$age <- df_summary[ind_summ,]$pinf_age_mo
  df_afq[ind_afq,]$sex <- df_summary[ind_summ,]$pinf_gender
  df_afq[ind_afq,]$pds <- df_summary[ind_summ,]$pds_shirtcliff
  
  # get, determine parent's scared group
  num_pscared <- df_summary[ind_summ,]$scaredp_sum
  if(num_pscared <= 10){
    group_pscared <- "Low"
  }else if(num_pscared > 10 & num_pscared < 25){
    group_pscared <- "Med"
  }else if(num_pscared >= 25){
    group_pscared <- "High"
  }
  
  # fill pscared info
  df_afq[ind_afq,]$pscared <- num_pscared
  df_afq[ind_afq,]$pscared_group <- group_pscared
  
  # get, determine pars-6 scared group
  num_pars6 <- df_summary[ind_summ,]$pars_6
  if (num_pars6 <= 3) {
    group_pars6 <- "Low"
  } else if (num_pars6 > 3 & num_pars6 <= 12) {
    group_pars6 <- "Med"
  } else if (num_pars6 > 12) {
    group_pars6 <- "High"
  }
  
  # fill pars info
  df_afq[ind_afq,]$pars6 <- num_pars6
  df_afq[ind_afq,]$pars6_group <- group_pars6
}

# save df
write_out <- paste0(out_dir, "AFQ_dataframe.csv")
write.table(df_afq, write_out, sep = ",", row.names=F)

# tract list
tract_list <- c("UNC_L", "CGC_L")

# subset df_afq
tract <- "UNC_L"
df_tract <- df_afq[which(df_afq$tractID == tract),]
df_tract <- df_tract[complete.cases(df_tract),]

# set factors
#   sex: 1 = F, 2 = M
#   pscared_group: 1 = High, 2 = Low, 3 = Med
#   pars6_group: 1 = High, 2 = Low, 3 = Med
df_tract$sex <- factor(df_tract$sex)
df_tract$pscared_group <- factor(df_tract$pscared_group)
df_tract$pars6_group <- factor(df_tract$pars6_group)

# plot raw
ggplot(data = df_tract) +
  geom_point(mapping = aes(x=nodeID, y=dti_fa, color=pscared_group), size=0.3) +
  geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=pscared_group))

# determine distribution
hist(df_tract$dti_fa)
descdist(df_tract$dti_fa, discrete=F) # normal or logistic

# gam
fit_normal <- bam(dti_fa ~ pscared_group +
    sex +
    s(nodeID, by = pscared_group, k = 20) +
    s(subjectID, bs = "re"),
  data = df_tract,
  family = gaussian(),
  method = "REML"
)
gam.check(fit_normal, rep = 1000)

# gam w/cov
fit_cov_normal <- bam(dti_fa ~ pscared_group +
    sex +
    s(nodeID, by = pscared_group, k = 30) +
    s(pds, by = sex) +
    s(subjectID, bs = "re"),
  data = df_tract,
  family = gaussian(),
  method = "REML"
)
gam.check(fit_cov_normal, rep = 1000)

# compare normal w/ cov_normal
compareML(fit_normal, fit_cov_normal) # cov better

# save better gam
gam_file = paste0(out_dir, "GAM_", tract, "_normal_cov_PScared.Rda")
saveRDS(fit_cov_normal, file = gam_file)

# generate predictions
df_pred <- predict.bam(
  fit_cov_normal,
  exclude_terms = c("pds", "sex", "subjectID"),
  values = list(pds = NULL, sex = NULL),
  se.fit = T,
  type = "response"
)

# convert predictions to dataframe
df_pred <- data.frame(
  PScared = df_tract$pscared_group,
  sex = df_tract$sex,
  subjectID = df_tract$subjectID,
  pscared = df_tract$pscared,
  nodeID = df_tract$nodeID,
  fit = df_pred$fit,
  se.fit = df_pred$se.fit
)


h_title <- paste0("GAM Fit of L. Uncinate FA Values")

# draw plot
ggplot(data = df_pred) +
  geom_smooth(mapping = aes(x = nodeID, y = fit, color = PScared)) +
  ggtitle(h_title) +
  ylab("Fit FA") +
  xlab("Tract Node") +
  theme(text = element_text(
    family = "Times New Roman", face = "bold", size = 14
  ))



# look at pars


# gam w/cov
fit_pars6 <- bam(dti_fa ~ pars6_group +
    sex +
    s(nodeID, by = pars6_group, k = 40) +
    s(pds, by = sex) +
    s(subjectID, bs = "re"),
  data = df_tract,
  family = gaussian(),
  method = "REML"
)
gam.check(fit_pars6, rep = 1000)

# compare normal w/ cov_normal
compareML(fit_pars6, fit_cov_normal) # cov better

# generate predictions
df_pred <- predict.bam(
  fit_pars6,
  exclude_terms = c("pds", "sex", "subjectID"),
  values = list(pds = NULL, sex = NULL),
  se.fit = T,
  type = "response"
)

# convert predictions to dataframe
df_pred <- data.frame(
  Pars6 = df_tract$pars6_group,
  sex = df_tract$sex,
  subjectID = df_tract$subjectID,
  pars6 = df_tract$pars6,
  nodeID = df_tract$nodeID,
  fit = df_pred$fit,
  se.fit = df_pred$se.fit
)


h_title <- paste0("GAM Fit of L. Uncinate FA Values")

# draw plot
ggplot(data = df_pred) +
  geom_smooth(mapping = aes(x = nodeID, y = fit, color = Pars6)) +
  ggtitle(h_title) +
  ylab("Fit FA") +
  xlab("Tract Node") +
  theme(text = element_text(
    family = "Times New Roman", face = "bold", size = 14
  ))


