####
#### Ogura Rice 2017
#### Compile data for MDR S-map
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.08.03
library(phyloseq); packageVersion("phyloseq") # 1.46.0, 2024.11.12
library(macam); packageVersion("macam") # 0.1.6, 2024.12.02, for MDR S-map

# Set random seeds (for reproduction)
set.seed(1234)

# Create output directory
(outdir <- macam::outdir_create())
dir.create("00_SessionInfo")


# ------------------------------------------------ #
# Load all data
# ------------------------------------------------ #
# All data
d_all <- readRDS("01_LoadAllDataOut/d_all.obj") %>% data.frame
ps_ecol <- readRDS("01_LoadAllDataOut/ps_ecol.obj")
d_ecol <- otu_table(ps_ecol) %>% data.frame

# Compile sample names
d_all[,"Sample_Name2"][is.na(d_all$Sample_Name2)] <- sprintf("NoDNA%02s", 1:sum(is.na(d_all$Sample_Name2)))
rownames(d_all) <- d_all$Sample_Name2

# Compiled UIC results
# Factors => Rice
uic_hgt <- readRDS("05_MDRsmapPrepOut/uic_hgt.obj")
uic_spd <- readRDS("05_MDRsmapPrepOut/uic_spd.obj")
uic_stm <- readRDS("05_MDRsmapPrepOut/uic_stm.obj")
# Rice => Factors
uic_hgt2 <- readRDS("05_MDRsmapPrepOut/uic_hgt2.obj")
uic_spd2 <- readRDS("05_MDRsmapPrepOut/uic_spd2.obj")
uic_stm2 <- readRDS("05_MDRsmapPrepOut/uic_stm2.obj")


# ------------------------------------------------ #
# MDR S-map
# ------------------------------------------------ #
# Construct the data.frame for MDR S-map
dim(d_all); rownames(d_all); dim(d_ecol); rownames(d_ecol)
## Compile data.frame
mdr_df <- d_all
edna_na_df <- matrix(NA, nrow = dim(d_all)[1], ncol = dim(d_ecol)[2]) %>% data.frame
colnames(edna_na_df) <- colnames(d_ecol)
rownames(edna_na_df) <- rownames(d_all)
## Assign eDNA values
edna_na_df[rownames(d_ecol), colnames(d_ecol)] <- d_ecol
mdr_df <- cbind(d_all, edna_na_df)

## MDR S-map loop
E_uic <- 0:5
# Tentative parameters
theta_range <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
lambda_range <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
alpha_range <- 0 #c(0,1)
regularized_i <- TRUE


# MDR S-map Loop
for (tar_var in c("mean_height_diff",
                  "mean_spad_diff",
                  "mean_n_stem_diff")) {
  # Set a target variable
  if (tar_var == "mean_height_diff") {
    uic_data <- uic_hgt
  } else if (tar_var == "mean_spad_diff") {
    uic_data <- uic_spd
  } else if (tar_var == "mean_n_stem_diff") {
    uic_data <- uic_stm
  }
  
  # Exclude several climate factors from UIC data
  clim_exclude <- c("clim_mean_temp","clim_max_temp","clim_min_temp","clim_daily_rainfall",
                    "clim_sunlight_hours","clim_cum_7days_rainfall")
  uic_data <- uic_data[!(unique(uic_data$cause_var) %in% clim_exclude),]
  
  # Make the data.frame for MDR S-map
  ## Standardize the data.frame
  mdr_df_i_tar <- mdr_df %>% select(as.name(tar_var)) %>% 
    apply(2, function(x) as.numeric(scale(x))) %>% data.frame
  mdr_df_i_var <- mdr_df %>% select(unique(uic_data$cause_var)) %>% 
    apply(2, function(x) as.numeric(scale(x))) %>% data.frame
  mdr_df_i <- cbind(mdr_df_i_tar, mdr_df_i_var)
  
  # Determine the optimal E of the target variable
  simp_x <- rUIC::simplex(mdr_df_i_tar, lib_var = tar_var, E = E_uic, tp = 1)
  (Ex <- with(simp_x, max(c(0, E[pval < 0.05]))))
  
  # Step 3: Make block to calculate multiview distance
  block_mvd <- make_block_mvd(mdr_df_i, uic_data, tar_var, E_effect_var = Ex,
                              th_var_colname = "fdr", p_threshold = 0.05,
                              include_var = "strongest_only", tp_adjust = 0)
  
  # Step. 4: Compute multiview distance
  multiview_dist <- compute_mvd(block_mvd, tar_var, E = Ex, tp = 1, distance_only = FALSE)
  
  # Step. 5: Do MDR S-map
  mdr_res_stats <- data.frame(N = NA, rho = NA, mae = NA, rmse = NA,
                              theta = NA, lambda = NA)
  mdr_res_stats <- NULL
  for (theta_i in theta_range) {
    for (lambda_i in lambda_range) {
      # Record time
      start_time <- proc.time()[3]
      
      # MDR S-map
      mdr_res_i <- s_map_mdr(block_mvd,
                             dist_w = multiview_dist$multiview_dist,
                             theta = theta_i,
                             lambda = lambda_i,
                             alpha = alpha_range,
                             tp = 1,
                             regularized = TRUE,
                             save_smap_coefficients = TRUE)
      mdr_res_i <- mdr_res_i$stats
      mdr_res_i$theta <- theta_i
      mdr_res_i$lambda <- lambda_i
      
      # Combine results
      mdr_res_stats <- rbind(mdr_res_stats, mdr_res_i)
      
      # Output message
      time_elapsed <- proc.time()[3] - start_time
      message(sprintf("%s, theta = %s and lambda = %s: Finished in %0.1f sec", tar_var, theta_i, lambda_i, time_elapsed))
    }
  }
  
  # Identify the best model and perform MDR S-map
  theta_best <- mdr_res_stats[which.min(mdr_res_stats$rmse),"theta"]
  lambda_best <- mdr_res_stats[which.min(mdr_res_stats$rmse),"lambda"]
  mdr_res <- s_map_mdr(block_mvd,
                       dist_w = multiview_dist$multiview_dist,
                       theta = theta_best,
                       lambda = lambda_best,
                       alpha = alpha_range,
                       tp = 1,
                       regularized = TRUE,
                       save_smap_coefficients = TRUE)
  
  # Add the original data.frame
  mdr_res$block_mvd <- block_mvd
  mdr_res$multiview_dist <- multiview_dist
  
  # Save results
  saveRDS(mdr_res, sprintf("06_MDRsmapOut/mdr_res_%s.obj", tar_var))
}

# ------------------------------------------------ #
# Save compiled data
# ------------------------------------------------ #
# Save workspace and session information
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))
macam::save_session_info()

