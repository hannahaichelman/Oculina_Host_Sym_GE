#### Custom functions for calculating gene expression plasticity based on Canonical Correspondence Analysis (CCA):
### This script can be sourced in your main script to run these functions. 
### Written by Colleen Bove, please reach out with questions (colleenbove@gmail.com)!



################################################################################
### Custom function to calculate PCA distances for plasticity

## To run the function, enter the following objects:
# PCAplast(pca = XXX, # the CCA dataframe containing the CCA eigenvalues
#          data = XXX, # the condition/treatment data corresponding to samples
#          sample_ID = "XXX", # the name of column that provide unique ID per sample (if blank, will pull rownames for this)
#          num_pca =  "XXX", # the number of CCAs to include in analysis (default is 'all', but you can specify another number with a minimum of 2 CCAs)
#          control_col = "XXX", # what the 'treatment' column is called
#          control_lvl = "XXX") # level of the treatment that is the 'control'


PCAplast <- function(pca, data, sample_ID = NA, num_pca = "all", control_col, control_lvl) {
  
  # rename the user input info
  pca_dist <- pca
  data_df <- data
  control_name <- control_col
  control_lvl <- control_lvl
  
  # check for correct number of CCAs provided
  if(class(num_pca) == "numeric") {
    if(num_pca < 2) { # will throw error if too few CCAs requested
      stop("please select more than 2 PCAs to calculate distance")
    } 
  }
  
  if(class(num_pca) == "numeric") {
    if(num_pca > (pca_dist %>% dplyr::select(starts_with("PC")) %>% ncol())) { # will throw error if too many PCAs requested
      stop(paste(num_pca, "PCAs requested for calculation, but only", (pca_dist %>% dplyr::select(starts_with("PC")) %>% ncol()), "PCAs available. Select appropriate number of PCAs for calculation."))
    } 
  }
  
  
  # oder the dataframe to ensure correct pairing after calculating distance 
  if(sample_ID %in% colnames(data_df)){
    data_df <- data_df[order(data_df[[sample_ID]]),]
  } else {
    data_df <- data_df[order(row.names(data_df)),]
  }
  
  
  # combine the datasets
  dist_df <- cbind(data_df, pca_dist) 
  
  # make dataframe of control colonies only
  mean_control <- dist_df %>%
    filter(dist_df[[control_name]] == list(control_lvl)[[1]]) %>% 
    rename_with(tolower) %>% # renames all pca's with lowercase 'PCA' (just to differentiate from all sample PCAs)
    dplyr::select(colnames(dist_df[control_name]), starts_with("pc")) %>% # select just the treatment and PCAs
    summarise_if(is.numeric, mean)
  
  # add the control CCA values to all samples 
  dist_df2 <- merge(dist_df, mean_control, all=TRUE)
  
  # again, reorder data
  if(sample_ID %in% colnames(data_df)){
    dist_df2 <- dist_df2[order(dist_df2[[sample_ID]]),]
  } else {
    rownames(dist_df2) <- rownames(data_df)
    dist_df2 <- dist_df2[order(row.names(dist_df2)),]
  }
  
  
  ### Calculate sample (PCA) distances from control (pca) using all PCAs
  # make dataframe to populate with pca distances
  full_calc_dist <- data.frame(control_name = dist_df2[control_name])
  
  if(num_pca == "all") {
    ## forloop that will calculate distances between control and sample for all PCAs (n will be total number)
    for(n in 1:(dist_df %>% dplyr::select(starts_with("PC")) %>% ncol())){
      # makes the PCA column name for control (lowercase) and sample (uppercase)
      PCA_col <- paste0("PC", n)
      pca_col <- paste0("pc", n)
      
      # pulls the PCA column for control (lowercase) and sample (uppercase)
      PCAx <- dist_df2[PCA_col]
      pcax <- dist_df2[pca_col]
      
      pca_calc_dist <- data.frame((PCAx - pcax)^2) # calculates the distance between 2 PCAs
      full_calc_dist <- cbind(full_calc_dist, pca_calc_dist) # add that distance to running dataframe
    }
  } else {
    ## forloop that will calculate distances between control and sample for SPECIFIED # of CCAs (n will be total number)
    for(n in 1:as.numeric(num_pca)){
      # makes the PCA column name for control (lowercase) and sample (uppercase)
      PCA_col <- paste0("PC", n)
      pca_col <- paste0("pc", n)
      
      # pulls the CCA column for control (lowercase) and sample (uppercase)
      PCAx <- dist_df2[PCA_col]
      pcax <- dist_df2[pca_col]
      
      pca_calc_dist <- data.frame((PCAx - pcax)^2) # calculates the distance between 2 CCAs
      full_calc_dist <- cbind(full_calc_dist, pca_calc_dist) # add that distance to running dataframe
    }
  }

  
  ## final distance calculation (adds all CCA distances and takes squareroot)
  distance <- full_calc_dist %>% 
    mutate(dis_sum = rowSums(across(where(is.numeric)))) %>% 
    mutate(dist = sqrt(dis_sum)) %>% 
    dplyr::select(matches("dist"))
  
  
  ## combine the calculated distance with the metadata and remove controls for final dataframe
  dist_df <- data_df %>% 
    bind_cols(distance)
}



