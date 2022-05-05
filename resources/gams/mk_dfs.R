library("tools")
library("stringr")

make_afq_df <- function(one_dir, data_dir) {
  # Add demographic, participant info to AFQ dataframe.
  # 
  # Added values are: age in month, sex, PARS-6,
  # Parent's/Child's SCARED, PARS/SCARED groups, 
  # primary diagnosis, diagnosis group.
  #
  # Arguments:
  #   one_dir (str) = path to directory containing emuR01_summary_latest.csv
  #   data_dir (str) = path to directory containing AFQ tract_profiles.csv
  #
  # Writes:
  #   data_dir/AFQ_dataframe.csv
  
  # get summary and afq dataframes
  df_summary <- read.csv(paste0(one_dir, "/emuR01_summary_latest.csv"))
  df_afq <- read.csv(paste0(data_dir, "/tract_profiles.csv"))
  
  # start age, sex, pds, parent's scared columns
  df_afq$age <- df_afq$sex <- df_afq$pds <-
    df_afq$pars6 <- df_afq$pars6_group <-
    df_afq$pscared <- df_afq$pscared_group <- 
    df_afq$cscared <- df_afq$dx <- df_afq$dx_group <-
    df_afq$lgi_neg <- df_afq$lgi_neu <- NA
  
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
    
    # fill age, sex, pds, lgi
    df_afq[ind_afq, ]$age <- df_summary[ind_summ, ]$pinf_age_mo
    df_afq[ind_afq, ]$sex <- df_summary[ind_summ, ]$pinf_gender
    df_afq[ind_afq, ]$pds <- df_summary[ind_summ, ]$pds_shirtcliff
    df_afq[ind_afq, ]$lgi_neg <- df_summary[ind_summ, ]$lgi_neg_1WK
    df_afq[ind_afq, ]$lgi_neu <- df_summary[ind_summ, ]$lgi_neu_1WK
    
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
    
    # Get diagnoses:
    #   GAD = GAD is dx 1, or dx GAD but SAD is not dx 1
    #   SAD = SAD, contains separation/social strings
    #   Con = No dx
    df_dx <- cbind(
      df_summary[ind_summ,]$adis_r_diag1,
      df_summary[ind_summ,]$adis_r_diag2,
      df_summary[ind_summ,]$adis_r_diag3,
      df_summary[ind_summ,]$adis_r_diag4
    )
    
    gad_prim <- grepl("Gen", df_dx[1, 1])
    gad_occur <- sum(grep("Gen", df_dx[1, ])) != 0
    
    sad_str <- c("Separation", "Social")
    sad_pos <- grep(paste(sad_str, collapse = "|"), df_dx[1,])
    sad_prim <- sad_occur <- F
    if (length(sad_pos) != 0){
      sad_occur <- T
      sad_prim <- sad_occur[1] == 1
    }
    
    if(gad_prim){
      h_dx <- "GAD"
    } else if (gad_occur & !sad_prim){
      h_dx <- "GAD"
    } else if (sad_occur){
      h_dx <- "SAD"
    }else{
      h_dx <- "Con"
    }
    df_afq[ind_afq, ]$dx <- h_dx
    
    # Set diagnosis group:
    #   Con = No dx
    #   Pat = GAD/SAD
    df_afq[ind_afq, ]$dx_group <- ifelse(h_dx == "Con", "Con", "Pat")
  }
  
  # save df, return for review
  write_out <- paste0(data_dir, "/AFQ_dataframe.csv")
  write.table(df_afq, write_out, sep = ",", row.names = F)
  return(df_afq)
}

make_epi_df <- function(raw_file, sess, task, beh_A, beh_B, roi, one_dir, data_dir) {
  # Organize roi/ppi output, make dataframes.
  #
  # Very similar to make_afq_df. Added values are: age in month, sex, PARS-6, 
  # Child's SCARED, Parent's SCARED, PARS/SCARED groups, primary diagnosis, 
  # diagnosis groups.
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
}

