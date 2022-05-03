# Notes ----
#
# This script will make AFNI-styled timing files from
# <proj_dir>/dset/**/func events.tsv files. Output files
# will be written to <write_dir>.
#
# Timing files made for study encoding trials preceding test
# behavioral response. Study non-responses will be "NR", study
# responses preceding a test non-response will be "PNR". Study
# response for stimuli not found in test will be "noTest"
#
# Update - timing files now also made for participant judgments in Study.
# Also updated timing file name to include decon string.
#
# Positional Arguments:
#   [1] = BIDS project directory
#   [2] = BIDS subject string
#   [3] = subject output directory
#   [4-5] = path to events files for study runs 1-2
#   [6-8] = path to events files for test runs 1-3


# Set Up -----
#
# Receive wrapped variables, set lists, and make switches.

study_list <- vector()
test_list <- vector()

args <- commandArgs(trailingOnly = T)
proj_dir <- args[1]
subj <- args[2]
write_dir <- args[3]

study_list[1] <- args[4]
study_list[2] <- args[5]

# test_list[1] <- args[6]
# test_list[2] <- args[7]
# test_list[3] <- args[8]

# # For testing
# proj_dir <- "/Volumes/homes/MaDLab/projects/McMakin_EMUR01/dset"
# subj <- "sub-4030"
# write_dir <- "/Users/nmuncy/Projects/emu_unc/data/timing_files/sub-4030/ses-S1"
# study_list[1] <- "/Volumes/homes/MaDLab/projects/McMakin_EMUR01/dset/sub-4030/ses-S1/func/sub-4030_ses-S1_task-study_run-1_events.tsv"
# study_list[2] <- "/Volumes/homes/MaDLab/projects/McMakin_EMUR01/dset/sub-4030/ses-S1/func/sub-4030_ses-S1_task-study_run-2_events.tsv"
# test_list[1] <- "/Volumes/homes/MaDLab/projects/McMakin_EMUR01/dset/sub-4030/ses-S2/func/sub-4030_ses-S2_task-test_run-1_events.tsv"
# test_list[2] <- "/Volumes/homes/MaDLab/projects/McMakin_EMUR01/dset/sub-4030/ses-S2/func/sub-4030_ses-S2_task-test_run-2_events.tsv"
# test_list[3] <- "/Volumes/homes/MaDLab/projects/McMakin_EMUR01/dset/sub-4030/ses-S2/func/sub-4030_ses-S2_task-test_run-3_events.tsv"

switch_string <- function(h_str) {
  # Rename Test behaviors to an AFNI length.
  #
  # Arguments:
  #   h_str (str) = behavior string from tsv file
  #
  # Returns:
  #   out_str (str) = AFNI-formatted string
  out_str <- switch(h_str,
    "neg_targ_ht" = "PnegTH",
    "neg_targ_ms" = "PnegTM",
    "neg_lure_cr" = "PnegLC",
    "neg_lure_fa" = "PnegLF",
    "neg_foil_cr" = "PnegFC",
    "neg_foil_fa" = "PnegFF",
    "neu_targ_ht" = "PneuTH",
    "neu_targ_ms" = "PneuTM",
    "neu_lure_cr" = "PneuLC",
    "neu_lure_fa" = "PneuLF",
    "neu_foil_cr" = "PneuFC",
    "neu_foil_fa" = "PneuFF",
    "pos_targ_ht" = "PposTH",
    "pos_targ_ms" = "PposTM",
    "pos_lure_cr" = "PposLC",
    "pos_lure_fa" = "PposLF",
    "pos_foil_cr" = "PposFC",
    "pos_foil_fa" = "PposFF",
    "NaN" = "PNR"
  )
  return(out_str)
}


# Make study, test dataframes ----
#
# Concatenate events files for each session into dataframes. Start
# by reaading in the first run, then rbind. Add run column for
# dividing into timing file rows.

# make master study df
df_study <- read.delim(study_list[1], sep = "\t", header = T)
df_study$run <- 1
for (run in 2:length(study_list)) {
  df <- read.delim(study_list[run], sep = "\t", header = T)
  df$run <- run
  df_study <- rbind(df_study, df)
  rm(df)
}

# # make master test df
# df_test <- read.delim(test_list[1], sep = "\t", header = T)
# df_test$run <- 1
# for (run in 2:length(test_list)) {
#   df <- read.delim(test_list[run], sep = "\t", header = T)
#   df$run <- run
#   df_test <- rbind(df_test, df)
#   rm(df)
# }


# Study judgment timing files ----
#
# Make timing files for participants stimulus valence ratings
# in Study session e.g. participant responds neg -> rNeg.

# get lsits, start count df
run_list <- unique(df_study$run)
beh_list <- unique(df_study$value)
df_count <- as.data.frame(matrix(NA, nrow = 2, ncol = length(beh_list)))
colnames(df_count) <- beh_list

# build timing files by run
for (run in run_list) {
  h_append <- ifelse(run == 1, F, T)
  df_run <- df_study[which(df_study$run == run), ]
  
  # get behavior counts
  for (beh in beh_list) {
    ind_beh <- which(df_run$value == beh)
    df_count[run, beh] <- length(ind_beh)
    
    # write asterisk for no behavior in run, otherwise get onset times
    if (length(ind_beh) == 0) {
      row_out <- "*"
    } else {
      row_out <- round(df_run[ind_beh, ]$onset, 1)
    }
    beh_out <- ifelse(beh == "NaN", "NR", beh)
    out_file <- paste0(
      write_dir, "/", "tf_task-study_decon-rVal_desc-", beh_out, "_events.txt"
    )
    cat(row_out, "\n", file = out_file, append = h_append, sep = "\t")
  }
}

# sum cols for df_count, write
df_count[3, ] <- colSums(df_count)
write.table(
  df_count[3, ],
  file = paste(write_dir, "beh_decon-rVal_counts.tsv", sep = "/"),
  sep = "\t",
  row.names = F
)

tf_list <- list.files(write_dir, pattern = "\\.txt$", full.names = T)
for (h_tf in tf_list) {
  h_check <- read.delim(h_tf, header = F)
  if (length(unique(h_check$V1)) == 1) {
    file.remove(h_tf)
  }
}

quit(save = "no", status = 0, runLast = F)


# Find test response to study stimulus ----
#
# Strip of stimulus similarity identifier to match study stimulus in test,
# get test response and convert it to AFNI format (switch string). Write
# test behavior as new column.

df_study$test_beh <- NA
stim_list <- df_study$stim_file
for (stim in stim_list) {

  # get study stim index, deal with study non responses
  ind_study <- grep(stim, df_study$stim_file)
  if (df_study[ind_study, ]$trial_type == "NaN") {
    df_study[ind_study, ]$test_beh <- "NR"
    next
  }

  # strip simularity identifier, find test response, convert
  stim_strip <- substr(stim, 1, 5)
  ind_test <- grep(stim_strip, df_test$stim_file)
  if (length(ind_test) == 0) {
    df_study[ind_study, ]$test_beh <- "noTest"
  } else {
    df_study[ind_study, ]$test_beh <-
      switch_string(df_test[ind_test, ]$trial_type)
  }
}


# Make Study preceding Test timing files ----
#
# Iterate through runs, behaviors. Find study onset times, write
# out to timing file. Account for no behavior type in run. Also count
# number of behavior types per run, total.

run_list <- unique(df_study$run)
beh_list <- unique(df_study$test_beh)

df_count <- as.data.frame(matrix(NA, nrow = 2, ncol = length(beh_list)))
colnames(df_count) <- beh_list

for (run in run_list) {
  h_append <- ifelse(run == 1, F, T)
  df_run <- df_study[which(df_study$run == run), ]
  for (beh in beh_list) {
    ind_beh <- which(df_run$test_beh == beh)
    df_count[run, beh] <- length(ind_beh)
    if (length(ind_beh) == 0) {
      row_out <- "*"
    } else {
      row_out <- round(df_run[ind_beh, ]$onset, 1)
    }
    out_file <- paste0(
      write_dir, "/", "tf_task-study_decon-precTest_desc-", beh, "_events.txt"
    )
    cat(row_out, "\n", file = out_file, append = h_append, sep = "\t")
  }
}
df_count[3, ] <- colSums(df_count)

# write df_count
write.table(
  df_count[3, ],
  file = paste(write_dir, "beh_decon-precTest_counts.tsv", sep = "/"),
  sep = "\t",
  row.names = F
)


# Remove empty timing files -----
#
# Remove files containing only asterisks, since
# AFNI will have a cow with those. Technically,
# I actually check to make sure there is more than
# one unique value in the first column. Surely this
# will always work.

tf_list <- list.files(write_dir, pattern = "\\.txt$", full.names = T)
for (h_tf in tf_list) {
  h_check <- read.delim(h_tf, header = F)
  if (length(unique(h_check$V1)) == 1) {
    file.remove(h_tf)
  }
}
