
# Set up ----
code_dir <- dirname(getwd())
data_dir <- paste0(code_dir, "/data")
source(paste0(code_dir, "/resources/gams/mk_dfs.R"))

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


# Make AFQ dataframe ----
make_afq_df(one_dir, data_dir)


# Make EPI dataframes ----
roi_list <- c("amgL", "amgR")
ns_list <- c("NSlacc", "NSldmpfc", "NSlsfs")
sess <- "ses-S1"
task <- "task-study"
decon <- "decon-rVal"
beh_A <- "neg"
beh_B <- "neu"

for (roi in c(roi_list, ns_list)) {
  
  # get coefs for decon roi analysis
  coef_file <- paste0(
    data_dir, "/Coefs_", sess, "_", task, "_", decon, "_", roi, ".txt"
  )
  make_epi_df(coef_file, sess, task, beh_A, beh_B, roi, one_dir, data_dir)
  
  # get coefs for decon amgL-PPI analysis
  if(roi == "amgL"){
    for(h_roi in ns_list){
      seed_roi = paste("amgL", h_roi, sep = "-")
      coef_file <- paste0(
        data_dir, "/Coefs_", sess, "_", task, "_", decon, "_", seed_roi, ".txt"
      )
      make_epi_df(
        coef_file, sess, task, beh_A, beh_B, seed_roi, one_dir, data_dir
      )
    }
  }
}
