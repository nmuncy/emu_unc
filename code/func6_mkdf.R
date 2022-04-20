library("tools")
library("stringr")

# General Notes ----
#
# Make dataframes of coefficients referencing output of func3 and func5
# ROI analyses.
#
# Notes:
#   - Requires OneDrive access for EMU summary dataset, and assumes amgL PPI seed.
#   - Omit "S" for PPI beh name, e.g. use "neg" for "Sneg" sub-brick
#
# Writes:
#   <data_dir>/df_<sess>_<task>_<decon>_ROI.csv
#
# Receives three positional arguments:
#   [1] = BIDS session (ses-S1)
#   [2] = BIDS task (task-study)
#   [3] = decon str (decon-rVal)
#   [4] = sub-brick behavior string (neg)
#   [5] = sub-brick behavior string (neu)
#
# Usage:
#   Rscript func6_mkdf.R ses-S1 task-study decon-rVal neg neu


# Receive Args ----
get_args <- commandArgs(trailingOnly = T)
if (length(get_args) != 5) {
  stop("Please see General Notes. Exiting.")
}
sess <- get_args[1]
task <- get_args[2]
decon <- get_args[3]
beh_A <- get_args[4]
beh_B <- get_args[5]

# # for testing
# sess <- "ses-S1"
# task <- "task-study"
# decon <- "decon-rVal"
# beh_A <- "neg"
# beh_B <- "neu"


# Functions ----
make_dataframe <- function(raw_file, sess, task, beh_A, beh_B, roi, one_dir, data_dir) {
  # Organize func3/func5 output, make dataframes.
  #
  # Added values are: age in month, sex, PARS-6, Child's SCARED,
  # Parent's SCARED, PARS/SCARED groups, primary diagnosis, diagnosis groups.
  #
  # Arguments:
  #   raw_file (str) = path to output text file of func4_roiAnalysis
  #   sess (str) = BIDS session
  #   task (str) = BIDS task
  #   beh_A (str) = sub-brick behavior to search
  #   beh_B (str) = sub-brick behavior to search
  #   ppi_seed (str) =  name of seed
  #   roi (str) = name of ROI, Seed-ROI for PPI
  #     e.g. decon ROI analysis = amgL, PPI ROI analysis= amgL-NSlacc
  #   one_dir (str) = path to directory containing emuR01_summary_latest.csv
  #   data_dir (str) = path to directory containing AFQ tract_profiles.csv
  #
  # Returns:
  #   df_out (dataframe) = wide-formatted, coefs & summary info
  #
  # Writes:
  #   data_dir/df_<sess>_<task>_<ppi_seed>-<roi>.csv

  # get ppi coef info
  df_raw <- read.delim(raw_file, header = F, sep = "\t")
  subj_list <- df_raw[str_detect(df_raw$V2, "File"), 1]
  beh_vecA <- df_raw[str_detect(df_raw$V2, beh_A), 3]
  beh_vecB <- df_raw[str_detect(df_raw$V2, beh_B), 3]
  df_out <- data.frame(
    "subj" <- subj_list,
    behA <- as.numeric(beh_vecA),
    behB <- as.numeric(beh_vecB)
  )
  colnames(df_out) <- c("subj", beh_A, beh_B)

  # add demographic, measure info
  df_out$sex <- df_out$age <- df_out$pars <- df_out$pscared <-
    df_out$cscared <- df_out$lgi <- df_out$dx <- df_out$dx_group <- NA
  df_summary <- read.csv(paste0(one_dir, "/emuR01_summary_latest.csv"))
  for (subj in subj_list) {

    # censure subjects who had <15 events modeled for both desired behaviors
    ind_out <- grep(subj, df_out$subj)
    
    subj_time_dir <- paste(data_dir, "timing_files", subj, sess, sep = "/")
    subj_count <- read.delim(
      paste0(subj_time_dir, "/beh_", decon, "_counts.tsv"),
      header = T,
      sep = "\t"
    )

    # check that all desired behs exist
    if (!beh_A %in% colnames(subj_count) || !beh_B %in% colnames(subj_count)) {
      df_out[ind_out, 2:4] <- NA
      next
    } else {

      # count beh events
      count_A <- as.numeric(subj_count[1, beh_A])
      count_B <- as.numeric(subj_count[1, beh_B])
      if (count_A < 15 || count_B < 15) {
        df_out[ind_out, 2:4] <- NA
        next
      }
    }

    # mine df_summary
    subj_int <- as.integer(gsub("sub-", "", subj))
    ind_summ <- grep(subj_int, df_summary$emu_study_id)
    if (length(ind_summ) == 0) {
      next
    }
    df_out[ind_out, ]$sex <- df_summary[ind_summ, ]$pinf_gender
    df_out[ind_out, ]$age <- df_summary[ind_summ, ]$pinf_age_mo
    df_out[ind_out, ]$pars <- df_summary[ind_summ, ]$pars_6
    df_out[ind_out, ]$pscared <- df_summary[ind_summ, ]$scaredp_sum
    df_out[ind_out, ]$cscared <- df_summary[ind_summ, ]$scaredc_sum
    df_out[ind_out, ]$lgi <- df_summary[ind_summ, ]$lgi_all_1WK

    # Get diagnoses:
    #   GAD = GAD is dx 1, or dx GAD but SAD is not dx 1
    #   SAD = SAD, contains separation/social strings
    #   Con = No dx
    df_dx <- cbind(
      df_summary[ind_summ, ]$adis_r_diag1,
      df_summary[ind_summ, ]$adis_r_diag2,
      df_summary[ind_summ, ]$adis_r_diag3,
      df_summary[ind_summ, ]$adis_r_diag4
    )

    gad_prim <- grepl("Gen", df_dx[1, 1])
    gad_occur <- sum(grep("Gen", df_dx[1, ])) != 0

    sad_str <- c("Separation", "Social")
    sad_pos <- grep(paste(sad_str, collapse = "|"), df_dx[1, ])
    sad_prim <- sad_occur <- F
    if (length(sad_pos) != 0) {
      sad_occur <- T
      sad_prim <- sad_occur[1] == 1
    }

    if (gad_prim) {
      h_dx <- "GAD"
    } else if (gad_occur & !sad_prim) {
      h_dx <- "GAD"
    } else if (sad_occur) {
      h_dx <- "SAD"
    } else {
      h_dx <- "Con"
    }
    df_out[ind_out, ]$dx <- h_dx

    # dx group - control or patient/experimental
    df_out[ind_out, ]$dx_group <- ifelse(h_dx == "Con", "Con", "Exp")
  }

  # save df
  write_out <- paste0(
    data_dir, "/df_", sess, "_", task, "_", decon, "_", roi, ".csv"
  )
  write.table(df_out, write_out, sep = ",", row.names = F)

  # return for checking
  return(df_out)
}


# Make dataframes ----
one_dir <- paste(
  Sys.getenv("HOME"),
  "Florida International University",
  "EMU Study - Documents",
  "EMU Data",
  "current_working_data",
  "datasets",
  "_full_working_dataset",
  sep = "/"
)
data_dir <- file_path_as_absolute(paste0(getwd(), "/../data"))

roi_list <- c("amgL", "amgR")
ns_list <- c("NSlacc", "NSldmpfc", "NSlsfs")

# get coefs for decon roi analysis
for (roi in c(roi_list, ns_list)) {
  
  coef_file <- paste0(
    data_dir, "/Coefs_", sess, "_", task, "_", decon, "_", roi, ".txt"
  )
  df_out <- make_dataframe(
    coef_file, sess, task, beh_A, beh_B, roi, one_dir, data_dir
  )
  
  # get coefs for decon amgL-PPI analysis
  if(roi == "amgL"){
    for(h_roi in ns_list){
      seed_roi = paste("amgL", h_roi, sep = "-")
      coef_file <- paste0(
        data_dir, "/Coefs_", sess, "_", task, "_", decon, "_", seed_roi, ".txt"
      )
      df_out <- make_dataframe(
        coef_file, sess, task, beh_A, beh_B, seed_roi, one_dir, data_dir
      )
    }
  }
}


