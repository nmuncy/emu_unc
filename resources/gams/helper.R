tract_fam <- function(tract) {
  # Set family for each tract.
  #
  # These were determined through comparing various models
  # via hist(), fitdistrplus::descdist, and itsadug::compareML.
  #
  # Arguments:
  #   tract (str) = AFQ tract name
  #
  # Returns:
  #   family (str) for use in DiffGamm::switch_family
  h_fam <- switch(tract,
                  "UNC_L" = "gamma",
                  "UNC_R" = "gaus",
                  "CGC_L" = "gaus",
                  "CGC_R" = "gaus",
  )
  return(h_fam)
}

switch_names <- function(name) {
  # Switch tract, Y-axis title names.
  #
  # Convert tract names to long form for titles,
  # convert memory, ppi names to long forms for
  # y-axis titles.
  #
  # Arguments:
  #   name (str) = AFQ tract name, behavior name, PPI seed
  #
  # Returns:
  #   name_long (str) = converted name
  x_name <- switch(name,
                   "UNC_L" = "L. Uncinate",
                   "UNC_R" = "R. Uncinate",
                   "CGC_L" = "L. Cingulum",
                   "CGC_R" = "R. Cingulum",
                   "lgi_neg" = "Negative LGI",
                   "lgi_neu" = "Neutral LGI",
                   "amgL" = "L. Amygdala",
                   "amgR" = "R. Amygdala",
                   "neg" = "Negative Scenes",
                   "neu" = "Neutral Scenes",
  )
  return(x_name)
}

write_gam_stats <- function(gam_obj, out_dir, gam_type, tract) {
  # Write summary stats of GAM object.
  #
  # Arguments:
  #   gam_obj (obj) = GAM object returned by mgcv
  #   out_dir (str) = path to output dir
  #   gam_type (str) = diffentiate type of model used
  #   tract (str) = AFQ tract name
  capture.output(
    summary(gam_obj),
    file = paste0(
      out_dir, "/Stats_GAM_", tract, "_", gam_type, ".txt"
    )
  )
}

write_compare_stats <- function(model_a, model_b, tract, out_dir, out_str) {
  # Write model comparison stats of two GAMs.
  #
  # Arguments:
  #   model_a (obj) = GAM object returned by mgcv
  #   model_b (obj) = another GAM object returned by mgcv
  #   tract (str) = AFQ tract name
  #   out_dir (str) = path to output dir
  #   out_str (str) = extra string for specifying name of out file
  capture.output(
    compareML(model_a, model_b),
    file = paste0(
      out_dir, "/Stats_GAM_", tract, "_compare_", out_str, ".txt"
    )
  )
}

tract_node <- function(tract) {
  # Match tract to node, for interaction investigation.
  #
  # Arguments:
  #   tract (str) = AFQ tract name
  #
  # Returns:
  #   id_node (int) = nodeID number
  id_node <- switch(tract,
                    "UNC_L" = 37,
                    "UNC_R" = 39,
                    "CGC_L" = 57,
                    "CGC_R" = 29
  )
  return(id_node)
}