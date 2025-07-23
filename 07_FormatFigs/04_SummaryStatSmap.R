####
#### Ogura Rice 2017
#### Summary stats
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.08.03
library(ggforce); packageVersion("ggforce") # 0.4.2, 2024.12.09
library(mgcv); packageVersion("mgcv") # 1.9.1, 2025.07.03
library(nlme); packageVersion("nlme") # 3.1.166, 2025.07.03

# Create output directory
dir.create("00_SessionInfo")


# ------------------------------------------------ #
# Load all figure data
# ------------------------------------------------ #
# MDR S-map results
s1 <- readRDS("data_robj/Smap_Height.obj")
s2 <- readRDS("data_robj/Smap_SPAD.obj")
s3 <- readRDS("data_robj/Smap_Nstem.obj")


# ------------------------------------------------ #
# Simple statistical test for the summary table
# ------------------------------------------------ #
# Parameter names
parm_names <- c("effect_var", "cause_var","int_est","int_tval","int_pval","nofert_est","nofert_tval","nofert_pval")


# ------------------------------------------------ #
# Effects on Rice Height
# ------------------------------------------------ #
# Extract data
s1_1_df <- s1[[1]]$data %>% data.frame
s1_1_cause_var <- s1_1_df$cause_var %>% unique

# Calculate overall effects
s1_1_df$overall_effect <- s1_1_df$smapc * s1_1_df$raw_value

# Prepare an output data.frame
s1_1_coefs_df <- data.frame(NULL)
s1_2_coefs_df <- data.frame(NULL)
s1_3_coefs_df <- data.frame(NULL)

# Main loop
for (i in s1_1_cause_var) {
  # Extract and compile the data
  df_tmp <- s1_1_df %>% filter(cause_var == i)
  df_tmp$date_num <- df_tmp$date %>% as.numeric
  df_tmp$date_num <- df_tmp$date_num - min(df_tmp$date_num) + 1
  df_tmp$ind_id_plot <- df_tmp$sample_code %>% str_split("_") %>%
    sapply(`[`, 4) %>% as.numeric %>% sprintf("P%02d", .) %>% factor
  df_tmp$treatment <- df_tmp$treatment %>% factor(levels = c("Conventional", "No Fertilizer"))
  
  ## Smap coefficients: Perform LME
  lme_tmp1 <- lme(smapc ~ treatment, random = list(ind_id_plot=~1), data = df_tmp)
  ## Collect summary information
  sumdf_tmp1 <- summary(lme_tmp1)$tTable
  coefs_tmp1 <- c("smapc", i, sumdf_tmp1[1,c(1,4,5)], sumdf_tmp1[2,c(1,4,5)])
  ## Rename data
  names(coefs_tmp1) <- parm_names
  
  ## Raw values: Perform LME
  lme_tmp2 <- lme(raw_value ~ treatment, random = list(ind_id_plot=~1), data = df_tmp)
  ## Collect summary information
  sumdf_tmp2 <- summary(lme_tmp2)$tTable
  coefs_tmp2 <- c("raw_value", i, sumdf_tmp2[1,c(1,4,5)], sumdf_tmp2[2,c(1,4,5)])
  ## Rename data
  names(coefs_tmp2) <- parm_names
  
  ## Overall effect: Perform LME
  lme_tmp3 <- lme(overall_effect ~ treatment, random = list(ind_id_plot=~1), data = df_tmp)
  ## Collect summary information
  sumdf_tmp3 <- summary(lme_tmp3)$tTable
  coefs_tmp3 <- c("overall_effect", i, sumdf_tmp3[1,c(1,4,5)], sumdf_tmp3[2,c(1,4,5)])
  ## Rename data
  names(coefs_tmp3) <- parm_names
  
  # Bind rows
  s1_1_coefs_df <- s1_1_coefs_df %>% bind_rows(coefs_tmp1)
  s1_2_coefs_df <- s1_2_coefs_df %>% bind_rows(coefs_tmp2)
  s1_3_coefs_df <- s1_3_coefs_df %>% bind_rows(coefs_tmp3)
  # Delete temporal objects
  rm(df_tmp); rm(lme_tmp1); rm(lme_tmp2); rm(lme_tmp3)
  rm(sumdf_tmp1); rm(sumdf_tmp2); rm(sumdf_tmp3)
  rm(coefs_tmp1); rm(coefs_tmp2); rm(coefs_tmp3)
}
# Change data class
s1_1_coefs_df[,3:8] <- apply(s1_1_coefs_df[,3:8], 2, as.numeric)
s1_2_coefs_df[,3:8] <- apply(s1_2_coefs_df[,3:8], 2, as.numeric)
s1_3_coefs_df[,3:8] <- apply(s1_3_coefs_df[,3:8], 2, as.numeric)

# Sort the data
tax_order1 <- s1[[1]]$scales$scales[[3]]$labels %>%
  str_sub(end = -2) %>% str_split("\\(") %>% sapply(`[`, 2) %>% rev
s1_1_coefs_df <- s1_1_coefs_df[match(tax_order1, s1_1_coefs_df$cause_var),]
s1_2_coefs_df <- s1_2_coefs_df[match(tax_order1, s1_2_coefs_df$cause_var),]
s1_3_coefs_df <- s1_3_coefs_df[match(tax_order1, s1_3_coefs_df$cause_var),]
s1_coefs_df <- rbind(s1_1_coefs_df, s1_2_coefs_df, s1_3_coefs_df)


# ------------------------------------------------ #
# Effects on Rice SPAD
# ------------------------------------------------ #
# Extract data
s2_1_df <- s2[[1]]$data %>% data.frame
s2_1_cause_var <- s2_1_df$cause_var %>% unique

# Calculate overall effects
s2_1_df$overall_effect <- s2_1_df$smapc * s2_1_df$raw_value

# Prepare an output data.frame
s2_1_coefs_df <- data.frame(NULL)
s2_2_coefs_df <- data.frame(NULL)
s2_3_coefs_df <- data.frame(NULL)

# Main loop
for (i in s2_1_cause_var) {
  # Extract and compile the data
  df_tmp <- s2_1_df %>% filter(cause_var == i)
  df_tmp$date_num <- df_tmp$date %>% as.numeric
  df_tmp$date_num <- df_tmp$date_num - min(df_tmp$date_num) + 1
  df_tmp$ind_id_plot <- df_tmp$sample_code %>% str_split("_") %>%
    sapply(`[`, 4) %>% as.numeric %>% sprintf("P%02d", .) %>% factor
  df_tmp$treatment <- df_tmp$treatment %>% factor(levels = c("Conventional", "No Fertilizer"))
  
  ## Smap coefficients: Perform LME
  lme_tmp1 <- lme(smapc ~ treatment, random = list(ind_id_plot=~1), data = df_tmp)
  ## Collect summary information
  sumdf_tmp1 <- summary(lme_tmp1)$tTable
  coefs_tmp1 <- c("smapc", i, sumdf_tmp1[1,c(1,4,5)], sumdf_tmp1[2,c(1,4,5)])
  ## Rename data
  names(coefs_tmp1) <- parm_names
  
  ## Raw values: Perform LME
  lme_tmp2 <- lme(raw_value ~ treatment, random = list(ind_id_plot=~1), data = df_tmp)
  ## Collect summary information
  sumdf_tmp2 <- summary(lme_tmp2)$tTable
  coefs_tmp2 <- c("raw_value", i, sumdf_tmp2[1,c(1,4,5)], sumdf_tmp2[2,c(1,4,5)])
  ## Rename data
  names(coefs_tmp2) <- parm_names
  
  ## Overall effect: Perform LME
  lme_tmp3 <- lme(overall_effect ~ treatment, random = list(ind_id_plot=~1), data = df_tmp)
  ## Collect summary information
  sumdf_tmp3 <- summary(lme_tmp3)$tTable
  coefs_tmp3 <- c("overall_effect", i, sumdf_tmp3[1,c(1,4,5)], sumdf_tmp3[2,c(1,4,5)])
  ## Rename data
  names(coefs_tmp3) <- parm_names
  
  # Bind rows
  s2_1_coefs_df <- s2_1_coefs_df %>% bind_rows(coefs_tmp1)
  s2_2_coefs_df <- s2_2_coefs_df %>% bind_rows(coefs_tmp2)
  s2_3_coefs_df <- s2_3_coefs_df %>% bind_rows(coefs_tmp3)
  # Delete temporal objects
  rm(df_tmp); rm(lme_tmp1); rm(lme_tmp2); rm(lme_tmp3)
  rm(sumdf_tmp1); rm(sumdf_tmp2); rm(sumdf_tmp3)
  rm(coefs_tmp1); rm(coefs_tmp2); rm(coefs_tmp3)
}
# Change data class
s2_1_coefs_df[,3:8] <- apply(s2_1_coefs_df[,3:8], 2, as.numeric)
s2_2_coefs_df[,3:8] <- apply(s2_2_coefs_df[,3:8], 2, as.numeric)
s2_3_coefs_df[,3:8] <- apply(s2_3_coefs_df[,3:8], 2, as.numeric)

# Sort the data
tax_order2 <- s2[[1]]$scales$scales[[3]]$labels %>%
  str_sub(end = -2) %>% str_split("\\(") %>% sapply(`[`, 2) %>% rev
s2_1_coefs_df <- s2_1_coefs_df[match(tax_order2, s2_1_coefs_df$cause_var),]
s2_2_coefs_df <- s2_2_coefs_df[match(tax_order2, s2_2_coefs_df$cause_var),]
s2_3_coefs_df <- s2_3_coefs_df[match(tax_order2, s2_3_coefs_df$cause_var),]
s2_coefs_df <- rbind(s2_1_coefs_df, s2_2_coefs_df, s2_3_coefs_df)


# ------------------------------------------------ #
# Effects on Rice N Stems
# ------------------------------------------------ #
# Extract data
s3_1_df <- s3[[1]]$data %>% data.frame
s3_1_cause_var <- s3_1_df$cause_var %>% unique

# Calculate overall effects
s3_1_df$overall_effect <- s3_1_df$smapc * s3_1_df$raw_value

# Prepare an output data.frame
s3_1_coefs_df <- data.frame(NULL)
s3_2_coefs_df <- data.frame(NULL)
s3_3_coefs_df <- data.frame(NULL)

# Main loop
for (i in s3_1_cause_var) {
  # Extract and compile the data
  df_tmp <- s3_1_df %>% filter(cause_var == i)
  df_tmp$date_num <- df_tmp$date %>% as.numeric
  df_tmp$date_num <- df_tmp$date_num - min(df_tmp$date_num) + 1
  df_tmp$ind_id_plot <- df_tmp$sample_code %>% str_split("_") %>%
    sapply(`[`, 4) %>% as.numeric %>% sprintf("P%02d", .) %>% factor
  df_tmp$treatment <- df_tmp$treatment %>% factor(levels = c("Conventional", "No Fertilizer"))
  
  ## Smap coefficients: Perform LME
  lme_tmp1 <- lme(smapc ~ treatment, random = list(ind_id_plot=~1), data = df_tmp)
  ## Collect summary information
  sumdf_tmp1 <- summary(lme_tmp1)$tTable
  coefs_tmp1 <- c("smapc", i, sumdf_tmp1[1,c(1,4,5)], sumdf_tmp1[2,c(1,4,5)])
  ## Rename data
  names(coefs_tmp1) <- parm_names
  
  ## Raw values: Perform LME
  lme_tmp2 <- lme(raw_value ~ treatment, random = list(ind_id_plot=~1), data = df_tmp)
  ## Collect summary information
  sumdf_tmp2 <- summary(lme_tmp2)$tTable
  coefs_tmp2 <- c("raw_value", i, sumdf_tmp2[1,c(1,4,5)], sumdf_tmp2[2,c(1,4,5)])
  ## Rename data
  names(coefs_tmp2) <- parm_names
  
  ## Overall effect: Perform LME
  lme_tmp3 <- lme(overall_effect ~ treatment, random = list(ind_id_plot=~1), data = df_tmp)
  ## Collect summary information
  sumdf_tmp3 <- summary(lme_tmp3)$tTable
  coefs_tmp3 <- c("overall_effect", i, sumdf_tmp3[1,c(1,4,5)], sumdf_tmp3[2,c(1,4,5)])
  ## Rename data
  names(coefs_tmp3) <- parm_names
  
  # Bind rows
  s3_1_coefs_df <- s3_1_coefs_df %>% bind_rows(coefs_tmp1)
  s3_2_coefs_df <- s3_2_coefs_df %>% bind_rows(coefs_tmp2)
  s3_3_coefs_df <- s3_3_coefs_df %>% bind_rows(coefs_tmp3)
  
  # Delete temporal objects
  rm(df_tmp); rm(lme_tmp1); rm(lme_tmp2); rm(lme_tmp3)
  rm(sumdf_tmp1); rm(sumdf_tmp2); rm(sumdf_tmp3)
  rm(coefs_tmp1); rm(coefs_tmp2); rm(coefs_tmp3)
}
# Change data class
s3_1_coefs_df[,3:8] <- apply(s3_1_coefs_df[,3:8], 2, as.numeric)
s3_2_coefs_df[,3:8] <- apply(s3_2_coefs_df[,3:8], 2, as.numeric)
s3_3_coefs_df[,3:8] <- apply(s3_3_coefs_df[,3:8], 2, as.numeric)

# Sort the data
tax_order3 <- s3[[1]]$scales$scales[[3]]$labels %>%
  str_sub(end = -2) %>% str_split("\\(") %>% sapply(`[`, 2) %>% rev
s3_1_coefs_df <- s3_1_coefs_df[match(tax_order3, s3_1_coefs_df$cause_var),]
s3_2_coefs_df <- s3_2_coefs_df[match(tax_order3, s3_2_coefs_df$cause_var),]
s3_3_coefs_df <- s3_3_coefs_df[match(tax_order3, s3_3_coefs_df$cause_var),]
s3_coefs_df <- rbind(s3_1_coefs_df, s3_2_coefs_df, s3_3_coefs_df)


# ------------------------------------------------ #
# Save compiled data
# ------------------------------------------------ #
# Save csv
write.csv(s1_coefs_df, "data_tbl/lme_smap_result_height.csv", row.names = FALSE)
write.csv(s2_coefs_df, "data_tbl/lme_smap_result_spad.csv", row.names = FALSE)
write.csv(s3_coefs_df, "data_tbl/lme_smap_result_nstem.csv", row.names = FALSE)

# Save workspace and session information
macam::save_session_info()

