library("tools")
library("stringr")

# General Notes ----
#
# Make dataframes of PPI coefficients referencing output of func3_roiAnalysis.
#
# Requires OneDrive access for EMU summary dataset.
#
# Writes: 
#   <data_dir>/df_<ppi_seed>-ROI.csv
#
# Receives three positional arguments:
#   [1] = ppi seed (amgL)
#   [2] = sub-brick behavior string (SnegLF)
#   [3] = sub-brick behavior string (SneuLF)
# 
# Usage:
#   Rscript func4_mkdf.R amgL SnegLF SneuLF

# Receive Args ----
get_args <- commandArgs(trailingOnly = T)
if(length(get_args) != 3){
  stop("Please see General Notes. Exiting.")
}
ppi_seed <- get_args[1]
beh_A <- get_args[2]
beh_B <- get_args[3]

# # for testing
# ppi_seed <- "amgL"
# beh_A <- "SnegLF"
# beh_B <- "SneuLF"


# Functions ----
make_dataframe <- function(
  raw_file, beh_A, beh_B, ppi_seed, roi, one_dir, data_dir
  ) {
  # Organize func3_roiAnalysis output, make dataframes.
  # 
  # Added values are: age in month, sex, PARS-6, Child's SCARED, and
  # Parent's SCARED values.
  #
  # Arguments:
  #   raw_file (str) = path to output text file of func3_roiAnalysis
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
  #   data_dir/df_<ppi_seed>-<roi>.csv

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
    df_out$cscared <- df_out$lgi <- NA
  df_summary <- read.csv(paste0(one_dir, "/emuR01_summary_latest.csv"))
  for (subj in subj_list) {
    subj_int <- as.integer(gsub("sub-", "", subj))
    ind_summ <- grep(subj_int, df_summary$emu_study_id)
    if (length(ind_summ) == 0) {
      next
    }
    ind_out <- grep(subj, df_out$subj)
    df_out[ind_out, ]$sex <- df_summary[ind_summ, ]$pinf_gender
    df_out[ind_out, ]$age <- df_summary[ind_summ, ]$pinf_age_mo
    df_out[ind_out, ]$pars <- df_summary[ind_summ, ]$pars_6
    df_out[ind_out, ]$pscared <- df_summary[ind_summ, ]$scaredp_sum
    df_out[ind_out, ]$cscared <- df_summary[ind_summ, ]$scaredc_sum
    df_out[ind_out, ]$lgi <- df_summary[ind_summ, ]$lgi_all_1WK
  }

  # save df
  write_out <- paste0(data_dir, "/df_", ppi_seed, "-", roi, ".csv")
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

vmpfc_file <- paste0(data_dir, "/Coefs_", ppi_seed, "-NSvmpfcL.txt")
df_out <- make_dataframe(
  vmpfc_file, beh_A, beh_B, ppi_seed, "lvmPFC", one_dir, data_dir
)

dmpfc_file <- paste0(data_dir, "/Coefs_", ppi_seed, "-NSdmpfcL.txt")
df_out <- make_dataframe(
  dmpfc_file, beh_A, beh_B, ppi_seed, "ldmPFC", one_dir, data_dir
)
