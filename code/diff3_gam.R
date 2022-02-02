library("ggplot2")
library("fitdistrplus")
library("mgcv")
library("itsadug")
library("tidymv")
library("dplyr")
library("mgcViz")


# Set paths ----
one_dir <- "/Users/nmuncy/Florida International University/EMU Study - Documents/EMU Data/current_working_data/datasets/_full_working_dataset/"
proj_dir <- "/Volumes/homes/MaDLab/projects/McMakin_EMUR01/"
data_dir <- paste0(proj_dir, "derivatives/emu_unc/reco_afq/")
out_dir <- paste0(proj_dir, "derivatives/emu_unc/analyses/")

# capture.output(sessionInfo(), file = paste0(data_dir, "R_session_info.txt"))

# for testing
out_dir <- data_dir <- proj_dir <- "/Users/nmuncy/Desktop/"



# Functions ----
make_dataframe <- function(one_dir, data_dir, out_dir) {
  # Add demographic, participant info to AFQ dataframe.
  # 
  # Added values are: age in month, sex, PARS-6, and
  # Parent's SCARED values.
  #
  # PARS-6 groups:
  #   Low <= 3
  #   3 < Med <= 12
  #   High > 12
  #
  # Parent's SCARED groups:
  #   Low <= 10
  #   10 < Med < 25
  #   High >= 25
  #
  # Arguments:
  #   one_dir (str) = path to directory containing emuR01_summary_latest.csv
  #   data_dir (str) = path to directory containing AFQ tract_profiles.csv
  #   out_dir (str) = output location for dataframe
  #
  # Returns:
  #   df_afq (dataframe) = long-formatted, AFQ + summary info
  #
  # Writes:
  #   out_dir/AFQ_dataframe.csv

  # get summary and afq dataframes
  df_summary <- read.csv(paste0(one_dir, "emuR01_summary_latest.csv"))
  df_afq <- read.csv(paste0(data_dir, "tract_profiles.csv"))

  # start age, sex, pds, parent's scared columns
  df_afq$age <- df_afq$sex <- df_afq$pds <-
    df_afq$pars6 <- df_afq$pars6_group <-
    df_afq$pscared <- df_afq$pscared_group <- NA

  # get list of afq subjects
  subj_list <- unique(df_afq$subjectID)

  # fill data for each subject
  for (subj in subj_list) {

    # determine indices of subjects in afq, summary dfs
    ind_afq <- which(df_afq$subjectID == subj)
    ind_summ <- which(df_summary$emu_study_id == subj)

    # make sure subj has summary info
    if (length(ind_summ) == 0) {
      print(paste("Missing summary data for", subj))
      next
    }

    # fill age, sex, pds
    df_afq[ind_afq, ]$age <- df_summary[ind_summ, ]$pinf_age_mo
    df_afq[ind_afq, ]$sex <- df_summary[ind_summ, ]$pinf_gender
    df_afq[ind_afq, ]$pds <- df_summary[ind_summ, ]$pds_shirtcliff

    # get, determine parent's scared number, group
    num_pscared <- df_summary[ind_summ, ]$scaredp_sum
    if (num_pscared <= 10) {
      group_pscared <- "Low"
    } else if (num_pscared > 10 & num_pscared < 25) {
      group_pscared <- "Med"
    } else if (num_pscared >= 25) {
      group_pscared <- "High"
    }

    # fill pscared info
    df_afq[ind_afq, ]$pscared <- num_pscared
    df_afq[ind_afq, ]$pscared_group <- group_pscared

    # get, determine pars-6 scared number, group
    num_pars6 <- df_summary[ind_summ, ]$pars_6
    if (num_pars6 <= 3) {
      group_pars6 <- "Low"
    } else if (num_pars6 > 3 & num_pars6 <= 12) {
      group_pars6 <- "Med"
    } else if (num_pars6 > 12) {
      group_pars6 <- "High"
    }

    # fill pars info
    df_afq[ind_afq, ]$pars6 <- num_pars6
    df_afq[ind_afq, ]$pars6_group <- group_pars6
  }
  
  # set factors
  #   sex: 1 = F, 2 = M
  #   pscared_group: 1 = High, 2 = Low, 3 = Med
  #   pars6_group: 1 = High, 2 = Low, 3 = Med
  df_tract$sex <- factor(df_tract$sex)
  df_tract$pscared_group <- factor(df_tract$pscared_group)
  df_tract$pars6_group <- factor(df_tract$pars6_group)

  # save df
  write_out <- paste0(out_dir, "AFQ_dataframe.csv")
  write.table(df_afq, write_out, sep = ",", row.names = F)

  return(df_afq)
}


# Determine data distribution, base GAM function ----
#
#

# get data
make_new_df <- F
if (make_new_df) {
  df_afq <- make_dataframe(one_dir, data_dir, out_dir)
} else {
  df_afq <- read.csv(paste0(out_dir, "AFQ_dataframe.csv"))
}

# tracts of interest
tract_list <- c("UNC_L", "CGC_L")

# subset df_afq, take complete cases, make factors
tract <- "UNC_L"
df_tract <- df_afq[which(df_afq$tractID == tract), ]
df_tract <- df_tract[complete.cases(df_tract), ]
df_tract$sex <- factor(df_tract$sex)
df_tract$subjectID <- factor(df_tract$subjectID)

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



# GAM with continuous Parent's SCARED ----
#
# Rather than looking for group differences, see if adding anxiety measure
# changes model fit. Then plot those two splines to find differences.

# intx of scared with fa
fit_pscared <- bam(dti_fa ~ sex +
                          s(subjectID, bs = "re") +
                          te(nodeID, pscared, bs = c("cr", "tp"), k = c(40, 10)),
                        data = df_tract,
                        family = gaussian(),
                        method = "fREML"
)
gam.check(fit_pscared, rep = 1000)
k.check(fit_pscared)

summary(fit_pscared)
compareML(fit_normal, fit_pscared)

# topo plot
plot(fit_pscared)
plot_pscared <- getViz(fit_pscared)
plot(sm(plot_pscared, 2))

# 3d plot
vis.gam(fit_pscared, view=c("nodeID", "pscared"), color = "topo", theta=330, phi=40)
vis.gam(fit_pscared, view=c("nodeID", "pscared"), theta=60)


# parent's scared - decompose intx term
fit_pscared_intx <- bam(dti_fa ~ sex +
                     s(subjectID, bs = "re") +
                     s(nodeID, bs = "cr", k = 40) +
                     s(pscared, k = 10) +
                     ti(nodeID, pscared, bs = c("cr", "tp"), k = c(40, 10)),
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
plot(sm(plot_pscared_intx, 3))
plot(sm(plot_pscared_intx, 4))

# vis.gam(fit_pscared_intx, view = c("nodeID", "pscared"), plot.type = "contour", color = "topo")
# vis.gam(fit_pscared_intx, view = c("nodeID", "pscared"), plot.type = "contour", color = "topo", too.far = 0.1)
# vis.gam(fit_pscared_intx, view = c("nodeID", "pscared"), plot.type = "contour", color = "topo", too.far = 0.05)


# GAM with continuous PARS-6 ----
#
#

fit_pars6 <- bam(dti_fa ~ sex +
      s(subjectID, bs = "re") +
      te(nodeID, pars6, bs = c("cr", "tp"), k = c(40, 10)),
    data = df_tract,
    family = gaussian(),
    method = "fREML"
)
gam.check(fit_pars6, rep = 1000)
k.check(fit_pars6)

summary(fit_pars6)
compareML(fit_normal, fit_pars6)

plot_pars6 <- getViz(fit_pars6)
plot(sm(plot_pars6, 2))


fit_pars6_intx <- bam(dti_fa ~ sex +
                   s(subjectID, bs = "re") +
                   s(nodeID, bs = "cr", k = 40) +
                   s(pars6) +
                   ti(nodeID, pars6, bs = c("cr", "tp"), k = c(40, 10)),
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
plot(sm(plot_pars6_intx, 3))
plot(sm(plot_pars6_intx, 4))

