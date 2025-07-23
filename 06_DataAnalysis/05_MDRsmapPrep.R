####
#### Ogura Rice 2017
#### Compile data for MDR S-map
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
# All data
d_all <- readRDS("01_LoadAllDataOut/d_all.obj") %>% data.frame
ps_ecol <- readRDS("01_LoadAllDataOut/ps_ecol.obj")
d_ecol <- otu_table(ps_ecol) %>% data.frame

# Compile sample names
d_all[,"Sample_Name2"][is.na(d_all$Sample_Name2)] <- sprintf("NoDNA%02s", 1:sum(is.na(d_all$Sample_Name2)))
rownames(d_all) <- d_all$Sample_Name2


# ------------------------------------------------ #
# Compile UIC data
# ------------------------------------------------ #
# Load UIC results
# (Object renames for simplicity)
## Climate variables v.s. rice
uic_hgt_clim <- readRDS("02_UIC_GrEnvOut/uic_height_res.obj")
uic_spd_clim <- readRDS("02_UIC_GrEnvOut/uic_spad_res.obj")
uic_stm_clim <- readRDS("02_UIC_GrEnvOut/uic_n_stem_res.obj")
## No causal effects from rice to climate

## Other biological variables v.s. rice
uic_hgt_others <- readRDS("03_UIC_GrOthersOut/uic_others_height.obj")
uic_spd_others <- readRDS("03_UIC_GrOthersOut/uic_others_spad.obj")
uic_stm_others <- readRDS("03_UIC_GrOthersOut/uic_others_n_stem.obj")
uic_hgt_others2 <- readRDS("03_UIC_GrOthersOut/uic_others_height_rev.obj")
uic_spd_others2 <- readRDS("03_UIC_GrOthersOut/uic_others_spad_rev.obj")
uic_stm_others2 <- readRDS("03_UIC_GrOthersOut/uic_others_n_stem_rev.obj")

## eDNA v.s. rice
uic_hgt_edna <- readRDS("04_UIC_GrDNAOut/uic_edna_height.obj")
uic_spd_edna <- readRDS("04_UIC_GrDNAOut/uic_edna_spad.obj")
uic_stm_edna <- readRDS("04_UIC_GrDNAOut/uic_edna_n_stem.obj")
uic_hgt_edna2 <- readRDS("04_UIC_GrDNAOut/uic_edna_height_rev.obj")
uic_spd_edna2 <- readRDS("04_UIC_GrDNAOut/uic_edna_spad_rev.obj")
uic_stm_edna2 <- readRDS("04_UIC_GrDNAOut/uic_edna_n_stem_rev.obj")


# ------------------------------------------------ #
# Choose the strongest causality for each pair
# ------------------------------------------------ #
# Define FDR threshold
fdr_th <- 0.05

# Choose the strongest causality with FDR < 0.05
## Rice Height (1st difference)
uic_hgt_clim <- uic_hgt_clim %>% filter(fdr <= fdr_th) %>% group_by(cause_var) %>% filter(ete == max(ete)) %>% data.frame
uic_hgt_others <- uic_hgt_others %>% filter(fdr <= fdr_th) %>% group_by(cause_var) %>% filter(ete == max(ete)) %>% data.frame
uic_hgt_others2 <- uic_hgt_others2 %>% filter(fdr <= fdr_th) %>% group_by(cause_var) %>% filter(ete == max(ete)) %>% data.frame
uic_hgt_edna <- uic_hgt_edna %>% filter(fdr <= fdr_th) %>% group_by(cause_var) %>% filter(ete == max(ete)) %>% data.frame
uic_hgt_edna2 <- uic_hgt_edna2 %>% filter(fdr <= fdr_th) %>% group_by(cause_var) %>% filter(ete == max(ete)) %>% data.frame
dim(uic_hgt_others); dim(uic_hgt_others2); dim(uic_hgt_edna); dim(uic_hgt_edna2)

## SPAD (1st difference)
uic_spd_clim <- uic_spd_clim %>% filter(fdr <= fdr_th) %>% group_by(cause_var) %>% filter(ete == max(ete)) %>% data.frame
#uic_spd_others <- uic_spd_others %>% filter(fdr <= fdr_th) %>% group_by(cause_var) %>% filter(ete == max(ete)) %>% data.frame # No causality
#uic_spd_others2 <- uic_spd_others2 %>% filter(fdr <= fdr_th) %>% group_by(cause_var) %>% filter(ete == max(ete)) %>% data.frame # No causality
uic_spd_edna <- uic_spd_edna %>% filter(fdr <= fdr_th) %>% group_by(cause_var) %>% filter(ete == max(ete)) %>% data.frame
uic_spd_edna2 <- uic_spd_edna2 %>% filter(fdr <= fdr_th) %>% group_by(cause_var) %>% filter(ete == max(ete)) %>% data.frame
dim(uic_spd_edna); dim(uic_spd_edna2)

## N of stems (1st difference)
uic_stm_clim <- uic_stm_clim %>% filter(fdr <= fdr_th) %>% group_by(cause_var) %>% filter(ete == max(ete)) %>% data.frame
#uic_stm_others <- uic_stm_others %>% filter(fdr <= fdr_th) %>% group_by(cause_var) %>% filter(ete == max(ete)) %>% data.frame # No causality
#uic_stm_others2 <- uic_stm_others2 %>% filter(fdr <= fdr_th) %>% group_by(cause_var) %>% filter(ete == max(ete)) %>% data.frame # No causality
## No causality between SPAD v.s. Others
uic_stm_edna <- uic_stm_edna %>% filter(fdr <= fdr_th) %>% group_by(cause_var) %>% filter(ete == max(ete)) %>% data.frame
uic_stm_edna2 <- uic_stm_edna2 %>% filter(fdr <= fdr_th) %>% group_by(cause_var) %>% filter(ete == max(ete)) %>% data.frame
dim(uic_stm_edna); dim(uic_stm_edna2)


# ------------------------------------------------ #
# Combine all the causal factors
# ------------------------------------------------ #
# Height
uic_hgt <- rbind(uic_hgt_clim, uic_hgt_others, uic_hgt_edna)
uic_hgt2 <- rbind(uic_hgt_others2, uic_hgt_edna2)
# SPAD
uic_spd <- rbind(uic_spd_clim, uic_spd_edna)
uic_spd2 <- rbind(uic_spd_edna2)
# N of stems
uic_stm <- rbind(uic_stm_clim, uic_stm_edna)
uic_stm2 <- rbind(uic_stm_edna2)


# ------------------------------------------------ #
# Save compiled data
# ------------------------------------------------ #
# R objects
saveRDS(uic_hgt, sprintf("%s/uic_hgt.obj", outdir))
saveRDS(uic_hgt2, sprintf("%s/uic_hgt2.obj", outdir))
saveRDS(uic_spd, sprintf("%s/uic_spd.obj", outdir))
saveRDS(uic_spd2, sprintf("%s/uic_spd2.obj", outdir))
saveRDS(uic_stm, sprintf("%s/uic_stm.obj", outdir))
saveRDS(uic_stm2, sprintf("%s/uic_stm2.obj", outdir))

# CSV files
write.csv(uic_hgt, sprintf("%s/CSV_uic_hgt.csv", outdir), row.names = FALSE)
write.csv(uic_hgt2, sprintf("%s/CSV_uic_hgt2.csv", outdir), row.names = FALSE)
write.csv(uic_spd, sprintf("%s/CSV_uic_spd.csv", outdir), row.names = FALSE)
write.csv(uic_spd2, sprintf("%s/CSV_uic_spd2.csv", outdir), row.names = FALSE)
write.csv(uic_stm, sprintf("%s/CSV_uic_stm.csv", outdir), row.names = FALSE)
write.csv(uic_stm2, sprintf("%s/CSV_uic_stm2.csv", outdir), row.names = FALSE)

# Save workspace and session information
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))
macam::save_session_info()

