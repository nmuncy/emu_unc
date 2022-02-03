library(tools)

# General Notes ----
#
#

make_dataframe <- function(one_dir, data_dir) {
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
  #
  # Returns:
  #   df_afq (dataframe) = long-formatted, AFQ + summary info
  #
  # Writes:
  #   data_dir/AFQ_dataframe.csv

  # get summary and afq dataframes
  df_summary <- read.csv(paste0(one_dir, "emuR01_summary_latest.csv"))
  df_afq <- read.csv(paste0(data_dir, "/tract_profiles.csv"))

  # start age, sex, pds, parent's scared columns
  df_afq$age <- df_afq$sex <- df_afq$pds <-
    df_afq$pars6 <- df_afq$pars6_group <-
    df_afq$pscared <- df_afq$pscared_group <- df_afq$cscared <- NA

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
    
    # get child's scared
    df_afq[ind_afq, ]$cscared <- df_summary[ind_summ, ]$scaredc_sum

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
  df_afq$sex <- factor(df_afq$sex)
  df_afq$pscared_group <- factor(df_afq$pscared_group)
  df_afq$pars6_group <- factor(df_afq$pars6_group)

  # save df
  write_out <- paste0(data_dir, "/AFQ_dataframe.csv")
  write.table(df_afq, write_out, sep = ",", row.names = F)

  return(df_afq)
}

# Make dataframe ----
one_dir <- "/Users/nmuncy/Florida International University/EMU Study - Documents/EMU Data/current_working_data/datasets/_full_working_dataset/"
data_dir <- file_path_as_absolute(paste0(getwd(), "/../data"))
df_afq <- make_dataframe(one_dir, data_dir)
