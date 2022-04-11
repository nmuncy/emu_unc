library("tools")
library("stringr")

# General Notes ----
#
# Make dataframes of PPI coefficients referencing output of func4_roiAnalysis.
#
# Requires OneDrive access for EMU summary dataset.
#
# Writes:
#   <data_dir>/df_<sess>_<task>_<ppi_seed>-ROI.csv
#
# Receives three positional arguments:
#   [1] = ppi seed (amgL)
#   [2] = BIDS session (ses-S1)
#   [3] = BIDS task (task-study)
#   [4] = sub-brick behavior string (SnegLF)
#   [5] = sub-brick behavior string (SneuLF)
#
# Usage:
#   Rscript func5_mkdf.R amgL ses-S1 task-study SPnegLF SPneuLF

# Receive Args ----
get_args <- commandArgs(trailingOnly = T)
if (length(get_args) != 5) {
  stop("Please see General Notes. Exiting.")
}
ppi_seed <- get_args[1]
sess <- get_args[2]
task <- get_args[3]
beh_A <- get_args[4]
beh_B <- get_args[5]

# # for testing
# ppi_seed <- "amgL"
# sess <- "ses-S1"
# task <- "task-study"
# beh_A <- "SPnegLF"
# beh_B <- "SPneuLF"


# Functions ----
beh_switch <- function(h_beh) {
  x_beh <- switch(h_beh,
    "SPnegLF" = "PnegLF",
    "SPneuLF" = "PneuLF"
  )
  return(x_beh)
}

make_dataframe <- function(raw_file, sess, task, beh_A, beh_B, ppi_seed, roi, one_dir, data_dir) {
  # Organize func4_roiAnalysis output, make dataframes.
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
  #   roi (str) = name of ROI
  #   one_dir (str) = path to directory containing emuR01_summary_latest.csv
  #   data_dir (str) = path to directory containing AFQ tract_profiles.csv
  #
  # Returns:
  #   df_out (dataframe) = wide-formatted, PPI coefs + summary info
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
  df_out[, paste0(beh_A, "-", beh_B)] <- NA
  df_out[, paste0(beh_A, "-", beh_B)] <- df_out[, beh_A] - df_out[, beh_B]
  rm(df_raw)

  # add demographic, measure info
  df_out$sex <- df_out$age <- df_out$pars <- df_out$pscared <-
    df_out$cscared <- df_out$lgi <- df_out$dx <- df_out$dx_group <- NA
  df_summary <- read.csv(paste0(one_dir, "/emuR01_summary_latest.csv"))
  for (subj in subj_list) {

    # censure subjects who had <7 behavior events modeled
    ind_out <- grep(subj, df_out$subj)

    subj_count <- read.delim(
      paste(data_dir, "timing_files", subj, sess, "beh_counts.tsv", sep = "/"),
      header = T,
      sep = "\t"
    )
    h_colA <- beh_switch(beh_A)
    h_colB <- beh_switch(beh_B)

    # check for behavior columns
    if (!h_colA %in% colnames(subj_count) || !h_colB %in% colnames(subj_count)) {
      df_out[ind_out, 2:4] <- NA
      next
    } else {

      # count beh events
      count_A <- as.numeric(subj_count[1, h_colA])
      count_B <- as.numeric(subj_count[1, h_colB])
      if (count_A < 7 || count_B < 7) {
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
    data_dir, "/df_", sess, "_", task, "_", ppi_seed, "-", roi, ".csv"
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

roi_list <- c("NSlacc", "NSldmpfc", "NSlsfs")
for (roi in roi_list) {
  coef_file <- paste0(
    data_dir, "/Coefs_", sess, "_", task, "_", ppi_seed, "-", roi, ".txt"
  )
  df_out <- make_dataframe(
    coef_file, sess, task, beh_A, beh_B, ppi_seed, roi, one_dir, data_dir
  )
  
  # write list of participants with sufficient number of trials
  df_complete <- df_out[complete.cases(df_out), ]
  subj_complete <- df_complete$subj
  subj_out <- paste0(
    data_dir, "/df_", sess, "_", task, "_", 
    ppi_seed, "-", roi, "_subjs-sufficient-trials.txt"
  )
  lapply(subj_complete, write, subj_out, append=TRUE, ncolumns=1)
}
