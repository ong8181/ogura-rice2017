####
#### Ogura Rice 2017
#### Load all data
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.08.03
library(phyloseq); packageVersion("phyloseq") # 1.46.0, 2024.11.12

# Set random seeds (for reproduction)
set.seed(1234)

# Create output directory
(outdir <- macam::outdir_create())
dir.create("00_SessionInfo")


# ------------------------------------------------ #
# Load all data
# ------------------------------------------------ #
# Kyoto climate data
d_clim <- read.csv("../01_RiceData/data/data_2017-climate-kyoto.csv")
## Add cumulative values
d_clim$cum_mean_temp <- 0; d_clim$cum_rainfall <- 0; d_clim$cum_sunlight_hours <- 0
for (i in 1:nrow(d_clim)) d_clim$cum_mean_temp[i] <- sum(d_clim$mean_temp[1:i])
for (i in 1:nrow(d_clim)) d_clim$cum_rainfall[i] <- sum(d_clim$daily_rainfall[1:i])
for (i in 1:nrow(d_clim)) d_clim$cum_sunlight_hours[i] <- sum(d_clim$sunlight_hours[1:i])
## Add 7-days (1-week) cumulative values
d_clim$cum_7days_temp <- 0; d_clim$cum_7days_rainfall <- 0; d_clim$cum_7days_sunlihg_hours <- 0
for (i in 1:nrow(d_clim)) {
  if ((i-6) < 0) {j <- 1} else {j <- i-6}
  d_clim$cum_7days_temp[i] <- sum(d_clim$mean_temp[j:i])
  d_clim$cum_7days_rainfall[i] <- sum(d_clim$daily_rainfall[j:i])
  d_clim$cum_7days_sunlihg_hours[i] <- sum(d_clim$sunlight_hours[j:i])
}
# Rice data
d_rice <- read.csv("../01_RiceData/data/data_2017-ogura_rice-monitoring.csv") %>% filter(treatment != "NC")
## Add mean SPAD, n_stem, and height, and 1st differences
d_rice$mean_n_stem <- apply(d_rice[,c("n_stem_1","n_stem_2","n_stem_3")], 1, function(x) mean(x, na.rm = TRUE))
d_rice$mean_spad <- apply(d_rice[,c("spad_1","spad_2","spad_3")], 1, function(x) mean(x, na.rm = TRUE))
d_rice$mean_height <- apply(d_rice[,c("height_1","height_2","height_3")], 1, function(x) mean(x, na.rm = TRUE))
d_rice$mean_n_stem_diff <- d_rice$mean_n_stem - lag(d_rice$mean_n_stem, n = 10) # Change rate
d_rice$mean_spad_diff <- d_rice$mean_spad - lag(d_rice$mean_spad, n = 10) # Change rate
d_rice$mean_height_diff <- d_rice$mean_height - lag(d_rice$mean_height, n = 10) # Change rate
# Volatile and phytoalexin data
d_volt <- read.csv("../02_Volitile/data/data_volatile.csv") %>% filter(treatment != "NC")
d_phyt <- read.csv("../02_Volitile/data/data_phytoalexin.csv") %>% filter(treatment != "NC")
# Plant cover data
d_grss <- read.csv("../04_GrassData/data/data_grass.csv") %>% filter(treatment != "NC")
# Insect data
d_insc <- read.csv("../05_InsectData/data/data_insect.csv") %>% filter(treatment != "NC")
# eDNA data
ps_ecol <- readRDS("../03_eDNA/02_AllDataCompile/04_RemoveSameSpOut/ps_comb_filt.obj")
d_ecol <- data.frame(sample_data(ps_ecol))
d_ecol$date <- ymd(d_ecol$date)
d_ecol$sample_code <- "TEMP"
d_ecol$sample_code[d_ecol$treatment == "CN"] <- "_Conv_"
d_ecol$sample_code[d_ecol$treatment == "NF"] <- "_NoFert_"
d_ecol$sample_code <- paste0("17_", str_sub(d_ecol$date, start = 6, end = 7), str_sub(d_ecol$date, start = 9, end = 10),
                             d_ecol$sample_code, as.numeric(str_sub(d_ecol$plot, start = 2, end = 3)))
sample_data(ps_ecol) <- sample_data(d_ecol)

# Check data structure
dim(d_clim); dim(d_rice); dim(d_volt); dim(d_phyt); dim(d_grss); dim(d_insc); dim(d_ecol)

# Combine all the data
## Prepare a new object
d_all <- d_rice
## Prepare a subset of climate data
d_clim_sub <- d_clim[match(d_all$date, d_clim$date),] %>%
  select(mean_temp, max_temp, min_temp, daily_rainfall, sunlight_hours,
         cum_mean_temp, cum_rainfall, cum_sunlight_hours,
         cum_7days_temp, cum_7days_rainfall, cum_7days_sunlihg_hours) %>% 
  rename_with(~ paste0("clim_", .x))
## Prepare a subset of volatile data
d_volt_sub <- d_volt[match(d_all$sample_code, d_volt$sample_code),] %>% 
  select(Z_3_Hexenal, E_2_Hexenal, Z_3_Hexen_1_ol, X2_4_Heptadienal, linalool, b_caryophyllene, MeSA) %>% 
  rename_with(~ paste0("psmb_", .x))
### Add diterpenoids
d_phyt <- d_phyt %>% mutate(diterpenoids_ug_g_FW = momilactone_A_ug_g_FW + momilactone_B_ug_g_FW + phytocassane_B_ug_g_FW + phytocassane_B3_ug_g_FW)
d_volt_sub$psmb_diterpenoids <- d_phyt[match(d_all$sample_code, d_phyt$sample_code),] %>% select(diterpenoids_ug_g_FW)
## Prepare a subset of grass data
d_grss_sub <- d_grss[match(d_all$sample_code, d_grss$sample_code),] %>% 
  select(cyper_1, cyper_2, cyper_3, cyper_4, poa_1, poa_2, poa_3, poa_4, poa_5,
         portulaca_1, euphorbia_1, brassica_1, aster_1, aster_2) %>% 
  rename_with(~ paste0("grss_", .x))
## Prepare a subset of insect data
d_insc_sub <- d_insc[match(d_all$sample_code, d_insc$sample_code),] %>% 
  select(other_planthopper, cicadellidae_spp, thysanoptera_spp) %>% 
  rename_with(~ paste0("insc_", .x))
## Prepare a subset of eDNA data
d_ecol_sub <- d_ecol[match(d_all$sample_code, d_ecol$sample_code),] %>% 
  select(PRO_dna_copy_sum, ITS_dna_copy_sum, COI_dna_copy_sum, EUK_dna_copy_sum) %>% 
  rename_with(~ paste0("edna_", .x))
## Check dimension
dim(d_clim_sub); dim(d_volt_sub); dim(d_grss_sub); dim(d_insc_sub); dim(d_ecol_sub)

## Combine
d_all <- cbind(d_all, d_clim_sub, d_volt_sub, d_grss_sub, d_insc_sub, d_ecol_sub)

## Change date format
d_all$date <- ymd(d_all$date)
d_all$time <- hm(d_all$time)

# Add phyloseq SampleName2 to d_all
d_all$Sample_Name2 <- d_ecol[match(d_all$sample_code, d_ecol$sample_code),] %>% select(Sample_Name2)

# Rename rows
rownames(d_all) <- 1:nrow(d_all)


# ------------------------------------------------ #
# Compile data for UIC
# ------------------------------------------------ #
# Sort d_all
d_all <- d_all %>%
  group_by(sampling_point) %>%
  arrange(date, .by_group = TRUE)

# Prepare for UIC
d_all <- d_all %>% mutate(uic_lib = TRUE, uic_group = as.factor(sampling_point))
d_all$uic_lib[d_all$date > "2017-08-30" | d_all$date < "2017-06-07"] <- FALSE
d_all$uic_group[d_all$date > "2017-08-30" | d_all$date < "2017-06-07"] <- NA


# ------------------------------------------------ #
# Save compiled data
# ------------------------------------------------ #
saveRDS(d_all, paste0(outdir, "/d_all.obj"))
saveRDS(ps_ecol, paste0(outdir, "/ps_ecol.obj"))

# Save workspace and session information
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))
macam::save_session_info()

