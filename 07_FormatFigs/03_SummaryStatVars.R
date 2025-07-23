####
#### Ogura Rice 2017
#### Summary stats
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.08.03
library(mgcv); packageVersion("mgcv") # 1.9.1

# Create output directory
dir.create("00_SessionInfo")
set.seed(1234)


# ------------------------------------------------ #
# Load all figure data
# ------------------------------------------------ #
# Climate
c1 <- readRDS("data_robj/Climate_all.obj")
# Rice
r1 <- readRDS("data_robj/Rice_Growth_rawdata.obj")
r2 <- readRDS("data_robj/Rice_Growth_diff.obj")
r3 <- readRDS("data_robj/Rice_Yield.obj")
r4 <- readRDS("data_robj/Rice_Yield_separate.obj")
# Specialized metabolites
m1 <- readRDS("data_robj/SpecializedMetabolites.obj")
# Grass
g1 <- readRDS("data_robj/Grass_temporal-pattern.obj")
# Insect
i1 <- readRDS("data_robj/Insect_temporal-pattern.obj")
# eDNA
e1 <- readRDS("data_robj/eDNA_temporal-pattern.obj")


# ------------------------------------------------ #
# Simple statistical test for the summary table
# ------------------------------------------------ #
# Parameter names
parm_names <- c("variable","int_est","int_t-val","int_p-val","no-fert_est","no-fert_t-val","no-fert_p-val",
                "date_edf_Conv","date_df_Conv","date_f-val_Conv","date_p-val_Conv",
                "date_edf_NonFert","date_df_NonFert","date_f-val_NonFert","date_p-val_NonFert")


# ------------------------------------------------ #
# Rice growth
# ------------------------------------------------ #
# Rice growth
r1_1_df <- r1[[1]]$data %>% na.omit() %>% data.frame
r1_1_df$date_num <- r1_1_df$date %>% as.numeric
r1_1_df$date_num <- r1_1_df$date_num - min(r1_1_df$date_num) + 1
r1_1_df$ind_id_plot <- paste0(r1_1_df$ind_id, "-P", sprintf("%02d", r1_1_df$sampling_point)) %>% factor
r1_1_df$sampling_point <- r1_1_df$sampling_point %>% factor
r1_1_df$treatment <- r1_1_df$treatment %>% factor
r1_1_gam <- gamm(height ~ s(date_num, by = treatment) + treatment, random = list(sampling_point = ~1, ind_id_plot=~1), data = r1_1_df, family = Gamma(link = "log"))
r1_1_sumdf <- summary(r1_1_gam$gam)$p.table
r1_1_coefs <- c("Height", r1_1_sumdf[1,c(1,3,4)], r1_1_sumdf[2,c(1,3,4)], 
                summary(r1_1_gam$gam)$s.table[1,], summary(r1_1_gam$gam)$s.table[2,])
names(r1_1_coefs) <- parm_names

# Rice SPAD
r1_2_df <- r1[[2]]$data %>% na.omit() %>% data.frame
r1_2_df$date_num <- r1_2_df$date %>% as.numeric
r1_2_df$date_num <- r1_2_df$date_num - min(r1_2_df$date_num) + 1
r1_2_df$ind_id_plot <- paste0(r1_2_df$ind_id, "-P", sprintf("%02d", r1_2_df$sampling_point)) %>% factor
r1_2_df$sampling_point <- r1_2_df$sampling_point %>% factor
r1_2_df$treatment <- r1_2_df$treatment %>% factor
r1_2_gam <- gamm(spad ~ s(date_num, by = treatment) + treatment, random = list(sampling_point = ~1, ind_id_plot=~1), data = r1_2_df, family = Gamma(link = "log"))
r1_2_sumdf <- summary(r1_2_gam$gam)$p.table
r1_2_coefs <- c("SPAD", r1_2_sumdf[1,c(1,3,4)], r1_2_sumdf[2,c(1,3,4)], 
                summary(r1_2_gam$gam)$s.table[1,], summary(r1_2_gam$gam)$s.table[2,])
names(r1_2_coefs) <- parm_names

# Rice nstem
r1_3_df <- r1[[3]]$data %>% na.omit() %>% data.frame
r1_3_df$date_num <- r1_3_df$date %>% as.numeric
r1_3_df$date_num <- r1_3_df$date_num - min(r1_3_df$date_num) + 1
r1_3_df$ind_id_plot <- paste0(r1_3_df$ind_id, "-P", sprintf("%02d", r1_3_df$sampling_point)) %>% factor
r1_3_df$sampling_point <- r1_3_df$sampling_point %>% factor
r1_3_df$treatment <- r1_3_df$treatment %>% factor
r1_3_gam <- gamm(n_stem ~ s(date_num, by = treatment) + treatment, random = list(sampling_point = ~1, ind_id_plot=~1), data = r1_3_df, family = Gamma(link = "log"))
r1_3_sumdf <- summary(r1_3_gam$gam)$p.table
r1_3_coefs <- c("n_stems", r1_3_sumdf[1,c(1,3,4)], r1_3_sumdf[2,c(1,3,4)], 
                summary(r1_3_gam$gam)$s.table[1,], summary(r1_3_gam$gam)$s.table[2,])
names(r1_3_coefs) <- parm_names


# ------------------------------------------------ #
# Rice growth (diff)
# ------------------------------------------------ #
# Rice growth
r2_1_df <- r2[[1]]$data %>% na.omit() %>% data.frame
r2_1_df$date_num <- r2_1_df$date %>% as.numeric
r2_1_df$date_num <- r2_1_df$date_num - min(r2_1_df$date_num) + 1
r2_1_df$ind_id_plot <- paste0(r2_1_df$ind_id, "-P", sprintf("%02d", r2_1_df$sampling_point)) %>% factor
r2_1_df$sampling_point <- r2_1_df$sampling_point %>% factor
r2_1_df$treatment <- r2_1_df$treatment %>% factor
r2_1_gam <- gamm(height_diff ~ s(date_num, by = treatment) + treatment, random = list(sampling_point = ~1, ind_id_plot=~1), data = r2_1_df)
r2_1_sumdf <- summary(r2_1_gam$gam)$p.table
r2_1_coefs <- c("Height_diff", r2_1_sumdf[1,c(1,3,4)], r2_1_sumdf[2,c(1,3,4)], 
                summary(r2_1_gam$gam)$s.table[1,], summary(r2_1_gam$gam)$s.table[2,])
names(r2_1_coefs) <- parm_names

# Rice SPAD
r2_2_df <- r2[[2]]$data %>% na.omit() %>% data.frame
r2_2_df$date_num <- r2_2_df$date %>% as.numeric
r2_2_df$date_num <- r2_2_df$date_num - min(r2_2_df$date_num) + 1
r2_2_df$ind_id_plot <- paste0(r2_2_df$ind_id, "-P", sprintf("%02d", r2_2_df$sampling_point)) %>% factor
r2_2_df$sampling_point <- r2_2_df$sampling_point %>% factor
r2_2_df$treatment <- r2_2_df$treatment %>% factor
r2_2_gam <- gamm(spad_diff ~ s(date_num, by = treatment) + treatment, random = list(sampling_point = ~1, ind_id_plot=~1), data = r2_2_df)
r2_2_sumdf <- summary(r2_2_gam$gam)$p.table
r2_2_coefs <- c("SPAD_diff", r2_2_sumdf[1,c(1,3,4)], r2_2_sumdf[2,c(1,3,4)], 
                summary(r2_2_gam$gam)$s.table[1,], summary(r2_2_gam$gam)$s.table[2,])
names(r2_2_coefs) <- parm_names

# Rice nstem
r2_3_df <- r2[[3]]$data %>% na.omit() %>% data.frame
r2_3_df$date_num <- r2_3_df$date %>% as.numeric
r2_3_df$date_num <- r2_3_df$date_num - min(r2_3_df$date_num) + 1
r2_3_df$ind_id_plot <- paste0(r2_3_df$ind_id, "-P", sprintf("%02d", r2_3_df$sampling_point)) %>% factor
r2_3_df$sampling_point <- r2_3_df$sampling_point %>% factor
r2_3_df$treatment <- r2_3_df$treatment %>% factor
r2_3_gam <- gamm(n_stem_diff ~ s(date_num, by = treatment) + treatment, random = list(sampling_point = ~1, ind_id_plot=~1), data = r2_3_df)
r2_3_sumdf <- summary(r2_3_gam$gam)$p.table
r2_3_coefs <- c("n_stems_diff", r2_3_sumdf[1,c(1,3,4)], r2_3_sumdf[2,c(1,3,4)], 
                summary(r2_3_gam$gam)$s.table[1,], summary(r2_3_gam$gam)$s.table[2,])
names(r2_3_coefs) <- parm_names


# ------------------------------------------------ #
# Plant specialized metabolites
# ------------------------------------------------ #
## Z_3_Hexenal
m3_1_df <- m1[[3]]$data %>% na.omit() %>% subset(name == "Z_3_Hexenal") %>% data.frame
m3_1_df$date_num <- m3_1_df$date %>% as.numeric
m3_1_df$date_num <- m3_1_df$date_num - min(m3_1_df$date_num) + 1
m3_1_df$ind_id_plot <- paste0("P", sprintf("%02d", m3_1_df$sampling_point)) %>% factor
m3_1_df$treatment <- m3_1_df$treatment %>% factor
m3_1_gam <- gamm(value ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = m3_1_df, family = Gamma(link = "log"))
m3_1_sumdf <- summary(m3_1_gam$gam)$p.table
m3_1_coefs <- c("Z_3_Hexenal", m3_1_sumdf[1,c(1,3,4)], m3_1_sumdf[2,c(1,3,4)], 
                summary(m3_1_gam$gam)$s.table[1,], summary(m3_1_gam$gam)$s.table[2,])
names(m3_1_coefs) <- parm_names

# E_2_Hexenal
m3_2_df <- m1[[3]]$data %>% na.omit() %>% subset(name == "E_2_Hexenal") %>% data.frame
m3_2_df$date_num <- m3_2_df$date %>% as.numeric
m3_2_df$date_num <- m3_2_df$date_num - min(m3_2_df$date_num) + 1
m3_2_df$ind_id_plot <- paste0("P", sprintf("%02d", m3_2_df$sampling_point)) %>% factor
m3_2_df$treatment <- m3_2_df$treatment %>% factor
m3_2_gam <- gamm(value ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = m3_2_df, family = Gamma(link = "log"))
m3_2_sumdf <- summary(m3_2_gam$gam)$p.table
m3_2_coefs <- c("E_2_Hexenal", m3_2_sumdf[1,c(1,3,4)], m3_2_sumdf[2,c(1,3,4)], 
                summary(m3_2_gam$gam)$s.table[1,], summary(m3_2_gam$gam)$s.table[2,])
names(m3_2_coefs) <- parm_names

# Z_3_Hexen_1_ol
m3_3_df <- m1[[3]]$data %>% na.omit() %>% subset(name == "Z_3_Hexen_1_ol") %>% data.frame
m3_3_df$date_num <- m3_3_df$date %>% as.numeric
m3_3_df$date_num <- m3_3_df$date_num - min(m3_3_df$date_num) + 1
m3_3_df$ind_id_plot <- paste0("P", sprintf("%02d", m3_3_df$sampling_point)) %>% factor
m3_3_df$treatment <- m3_3_df$treatment %>% factor
m3_3_gam <- gamm(value ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = m3_3_df, family = Gamma(link = "log"))
m3_3_sumdf <- summary(m3_3_gam$gam)$p.table
m3_3_coefs <- c("Z_3_Hexen_1_ol", m3_3_sumdf[1,c(1,3,4)], m3_3_sumdf[2,c(1,3,4)], 
                summary(m3_3_gam$gam)$s.table[1,], summary(m3_3_gam$gam)$s.table[2,])
names(m3_3_coefs) <- parm_names

# X2_4_Heptadienal
m3_4_df <- m1[[3]]$data %>% na.omit() %>% subset(name == "X2_4_Heptadienal") %>% data.frame
m3_4_df$date_num <- m3_4_df$date %>% as.numeric
m3_4_df$date_num <- m3_4_df$date_num - min(m3_4_df$date_num) + 1
m3_4_df$ind_id_plot <- paste0("P", sprintf("%02d", m3_4_df$sampling_point)) %>% factor
m3_4_df$treatment <- m3_4_df$treatment %>% factor
m3_4_gam <- gamm(value ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = m3_4_df, family = Gamma(link = "log"))
m3_4_sumdf <- summary(m3_4_gam$gam)$p.table
m3_4_coefs <- c("X2_4_Heptadienal", m3_4_sumdf[1,c(1,3,4)], m3_4_sumdf[2,c(1,3,4)], 
                summary(m3_4_gam$gam)$s.table[1,], summary(m3_4_gam$gam)$s.table[2,])
names(m3_4_coefs) <- parm_names

# linalool
m3_5_df <- m1[[3]]$data %>% na.omit() %>% subset(name == "linalool") %>% data.frame
m3_5_df$date_num <- m3_5_df$date %>% as.numeric
m3_5_df$date_num <- m3_5_df$date_num - min(m3_5_df$date_num) + 1
m3_5_df$ind_id_plot <- paste0("P", sprintf("%02d", m3_5_df$sampling_point))
m3_5_df$treatment <- m3_5_df$treatment %>% factor
m3_5_gam <- gamm(value ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = m3_5_df, family = Gamma(link = "log"))
m3_5_sumdf <- summary(m3_5_gam$gam)$p.table
m3_5_coefs <- c("linalool", m3_5_sumdf[1,c(1,3,4)], m3_5_sumdf[2,c(1,3,4)], 
                summary(m3_5_gam$gam)$s.table[1,], summary(m3_5_gam$gam)$s.table[2,])
names(m3_5_coefs) <- parm_names

# b_caryophyllene
m3_6_df <- m1[[3]]$data %>% na.omit() %>% subset(name == "b_caryophyllene") %>% data.frame
m3_6_df$date_num <- m3_6_df$date %>% as.numeric
m3_6_df$date_num <- m3_6_df$date_num - min(m3_6_df$date_num) + 1
m3_6_df$ind_id_plot <- paste0("P", sprintf("%02d", m3_6_df$sampling_point)) %>% factor
m3_6_df$treatment <- m3_6_df$treatment %>% factor
m3_6_gam <- gamm(value ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = m3_6_df, family = Gamma(link = "log"))
m3_6_sumdf <- summary(m3_6_gam$gam)$p.table
m3_6_coefs <- c("b_caryophyllene", m3_6_sumdf[1,c(1,3,4)], m3_6_sumdf[2,c(1,3,4)], 
                summary(m3_6_gam$gam)$s.table[1,], summary(m3_6_gam$gam)$s.table[2,])
names(m3_6_coefs) <- parm_names

# MeSA
m3_7_df <- m1[[3]]$data %>% na.omit() %>% subset(name == "MeSA") %>% data.frame
m3_7_df$date_num <- m3_7_df$date %>% as.numeric
m3_7_df$date_num <- m3_7_df$date_num - min(m3_7_df$date_num) + 1
m3_7_df$ind_id_plot <- paste0("P", sprintf("%02d", m3_7_df$sampling_point)) %>% factor
m3_7_df$treatment <- m3_7_df$treatment %>% factor
m3_7_gam <- gamm(value ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = m3_7_df, family = Gamma(link = "log"))
m3_7_sumdf <- summary(m3_7_gam$gam)$p.table
m3_7_coefs <- c("MeSA", m3_7_sumdf[1,c(1,3,4)], m3_7_sumdf[2,c(1,3,4)], 
                summary(m3_7_gam$gam)$s.table[1,], summary(m3_7_gam$gam)$s.table[2,])
names(m3_7_coefs) <- parm_names

# diterpenoids_ug_g_FW
m3_8_df <- m1[[3]]$data %>% na.omit() %>% subset(name == "diterpenoids_ug_g_FW") %>% data.frame
m3_8_df$date_num <- m3_8_df$date %>% as.numeric
m3_8_df$date_num <- m3_8_df$date_num - min(m3_8_df$date_num) + 1
m3_8_df$ind_id_plot <- paste0("P", sprintf("%02d", m3_8_df$sampling_point)) %>% factor
m3_8_df$treatment <- m3_8_df$treatment %>% factor
m3_8_df$value[m3_8_df$value == 0] <- 1e-03
m3_8_gam <- gamm(value ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = m3_8_df, family = Gamma(link = "log"), niterPQL = 50)
m3_8_sumdf <- summary(m3_8_gam$gam)$p.table
m3_8_coefs <- c("diterpenoids_ug_g_FW", m3_8_sumdf[1,c(1,3,4)], m3_8_sumdf[2,c(1,3,4)], 
                summary(m3_8_gam$gam)$s.table[1,], summary(m3_8_gam$gam)$s.table[2,])
names(m3_8_coefs) <- parm_names


# ------------------------------------------------ #
# Plant cover (%)
# ------------------------------------------------ #
# cyper_all
g1_1_df <- g1[[2]]$data %>% na.omit() %>% subset(name == "cyper_all") %>% data.frame
g1_1_df$date_num <- g1_1_df$date %>% as.numeric
g1_1_df$date_num <- g1_1_df$date_num - min(g1_1_df$date_num) + 1
g1_1_df$ind_id_plot <- paste0("P", sprintf("%02d", g1_1_df$sampling_point)) %>% factor
g1_1_df$treatment <- g1_1_df$treatment %>% factor
g1_1_df$value[g1_1_df$value == 0] <- 1e-03
g1_1_gam <- gamm(value ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = g1_1_df, family = Gamma(link = "log"), niterPQL = 50)
g1_1_sumdf <- summary(g1_1_gam$gam)$p.table
g1_1_coefs <- c("cyper_all", g1_1_sumdf[1,c(1,3,4)], g1_1_sumdf[2,c(1,3,4)], 
                summary(g1_1_gam$gam)$s.table[1,], summary(g1_1_gam$gam)$s.table[2,])
names(g1_1_coefs) <- parm_names

# poa_all
g1_2_df <- g1[[2]]$data %>% na.omit() %>% subset(name == "poa_all") %>% data.frame
g1_2_df$date_num <- g1_2_df$date %>% as.numeric
g1_2_df$date_num <- g1_2_df$date_num - min(g1_2_df$date_num) + 1
g1_2_df$ind_id_plot <- paste0("P", sprintf("%02d", g1_2_df$sampling_point)) %>% factor
g1_2_df$treatment <- g1_2_df$treatment %>% factor
g1_2_df$value[g1_2_df$value == 0] <- 1e-03
g1_2_gam <- gamm(value ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = g1_2_df, family = Gamma(link = "log"), niterPQL = 50)
g1_2_sumdf <- summary(g1_2_gam$gam)$p.table
g1_2_coefs <- c("poa_all", g1_2_sumdf[1,c(1,3,4)], g1_2_sumdf[2,c(1,3,4)], 
                summary(g1_2_gam$gam)$s.table[1,], summary(g1_2_gam$gam)$s.table[2,])
names(g1_2_coefs) <- parm_names

# others_all
g1_3_df <- g1[[2]]$data %>% na.omit() %>% subset(name == "others_all") %>% data.frame
g1_3_df$date_num <- g1_3_df$date %>% as.numeric
g1_3_df$date_num <- g1_3_df$date_num - min(g1_3_df$date_num) + 1
g1_3_df$ind_id_plot <- paste0("P", sprintf("%02d", g1_3_df$sampling_point)) %>% factor
g1_3_df$treatment <- g1_3_df$treatment %>% factor
g1_3_df$value[g1_3_df$value == 0] <- 1e-03
g1_3_gam <- gamm(value ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = g1_3_df, family = Gamma(link = "log"), niterPQL = 50)
g1_3_sumdf <- summary(g1_3_gam$gam)$p.table
g1_3_coefs <- c("others_all", g1_3_sumdf[1,c(1,3,4)], g1_3_sumdf[2,c(1,3,4)], 
                summary(g1_3_gam$gam)$s.table[1,], summary(g1_3_gam$gam)$s.table[2,])
names(g1_3_coefs) <- parm_names


# ------------------------------------------------ #
# Insect
# ------------------------------------------------ #
# cicadellidae_spp
i1_1_df <- i1$data %>% na.omit() %>% subset(name == "cicadellidae_spp") %>% data.frame
i1_1_df$date_num <- i1_1_df$date %>% as.numeric
i1_1_df$date_num <- i1_1_df$date_num - min(i1_1_df$date_num) + 1
i1_1_df$ind_id_plot <- paste0("P", sprintf("%02d", i1_1_df$sampling_point)) %>% factor
i1_1_df$treatment <- i1_1_df$treatment %>% factor
i1_1_df$value[i1_1_df$value == 0] <- 1e-03
i1_1_gam <- gamm(value ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = i1_1_df, family = Gamma(link = "log"), niterPQL = 50)
i1_1_sumdf <- summary(i1_1_gam$gam)$p.table
i1_1_coefs <- c("cicadellidae_spp", i1_1_sumdf[1,c(1,3,4)], i1_1_sumdf[2,c(1,3,4)],
                summary(i1_1_gam$gam)$s.table[1,], summary(i1_1_gam$gam)$s.table[2,])
names(i1_1_coefs) <- parm_names

# thysanoptera_spp
i1_2_df <- i1$data %>% na.omit() %>% subset(name == "thysanoptera_spp") %>% data.frame
i1_2_df$date_num <- i1_2_df$date %>% as.numeric
i1_2_df$date_num <- i1_2_df$date_num - min(i1_2_df$date_num) + 1
i1_2_df$ind_id_plot <- paste0("P", sprintf("%02d", i1_2_df$sampling_point)) %>% factor
i1_2_df$treatment <- i1_2_df$treatment %>% factor
i1_2_df$value[i1_2_df$value == 0] <- 1e-03
i1_2_gam <- gamm(value ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = i1_2_df, family = Gamma(link = "log"), niterPQL = 50)
i1_2_sumdf <- summary(i1_2_gam$gam)$p.table
i1_2_coefs <- c("thysanoptera_spp", i1_2_sumdf[1,c(1,3,4)], i1_2_sumdf[2,c(1,3,4)],
                summary(i1_2_gam$gam)$s.table[1,], summary(i1_2_gam$gam)$s.table[2,])
names(i1_2_coefs) <- parm_names

# other_planthopper
i1_3_df <- i1$data %>% na.omit() %>% subset(name == "other_planthopper") %>% data.frame
i1_3_df$date_num <- i1_3_df$date %>% as.numeric
i1_3_df$date_num <- i1_3_df$date_num - min(i1_3_df$date_num) + 1
i1_3_df$ind_id_plot <- paste0("P", sprintf("%02d", i1_3_df$sampling_point)) %>% factor
i1_3_df$treatment <- i1_3_df$treatment %>% factor
i1_3_df$value[i1_3_df$value == 0] <- 1e-03
i1_3_gam <- gamm(value ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = i1_3_df, family = Gamma(link = "log"), niterPQL = 100)
i1_3_sumdf <- summary(i1_3_gam$gam)$p.table
i1_3_coefs <- c("other_planthopper", i1_3_sumdf[1,c(1,3,4)], i1_3_sumdf[2,c(1,3,4)],
                summary(i1_3_gam$gam)$s.table[1,], summary(i1_3_gam$gam)$s.table[2,])
names(i1_3_coefs) <- parm_names


# ------------------------------------------------ #
# eDNA
# ------------------------------------------------ #
# Total eDNA abundance
e1_1_df <- e1[[1]]$data %>% na.omit() %>% group_by(date, treatment, plot) %>% summarize(value = sum(Abundance)) %>% data.frame
e1_1_df$date_num <- e1_1_df$date %>% as.numeric
e1_1_df$date_num <- e1_1_df$date_num - min(e1_1_df$date_num) + 1
e1_1_df$ind_id_plot <- e1_1_df$plot %>% factor
e1_1_df$treatment <- e1_1_df$treatment %>% factor
e1_1_gam <- gamm(value ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = e1_1_df, family = Gamma(link = "log"))
e1_1_sumdf <- summary(e1_1_gam$gam)$p.table
e1_1_coefs <- c("total_eDNA_abundance", e1_1_sumdf[1,c(1,3,4)], e1_1_sumdf[2,c(1,3,4)],
                summary(e1_1_gam$gam)$s.table[1,], summary(e1_1_gam$gam)$s.table[2,])
names(e1_1_coefs) <- parm_names

# eDNA diversity (N. OTU)
e1_2_df <- e1[[2]]$data %>% na.omit() %>% data.frame
e1_2_df$date_num <- e1_2_df$date %>% as.numeric
e1_2_df$date_num <- e1_2_df$date_num - min(e1_2_df$date_num) + 1
e1_2_df$ind_id_plot <- e1_2_df$plot %>% factor
e1_2_df$treatment <- e1_2_df$treatment %>% factor
e1_2_gam <- gamm(n_asv ~ s(date_num, by = treatment) + treatment, random = list(ind_id_plot=~1), data = e1_2_df, family = Gamma(link = "log"))
e1_2_sumdf <- summary(e1_2_gam$gam)$p.table
e1_2_coefs <- c("eDNA_n_asv", e1_2_sumdf[1,c(1,3,4)], e1_2_sumdf[2,c(1,3,4)],
                summary(e1_2_gam$gam)$s.table[1,], summary(e1_2_gam$gam)$s.table[2,])
names(e1_2_coefs) <- parm_names


# ------------------------------------------------ #
# Rice yields (approximately Gaussian distribution)
# ------------------------------------------------ #
# Prepare output data.frame
rice_parms <- c("variable", "fval", "pval")
yield_stat_df <- data.frame(NULL)

# Total number of rice head
r4_1_df <- r4[[1]]$data %>% data.frame
r4_1_df$treatment <- factor(r4_1_df$treatment)
r4_1_res <- aov(value ~ treatment, data = r4_1_df)
r4_1_coefs <- summary(r4_1_res)[[1]] %>% as.matrix
r4_1_coefs <- c("total_rice_heads", r4_1_coefs[1,c(4,5)])
names(r4_1_coefs) <- rice_parms

# Total number of rice head
r4_2_df <- r4[[2]]$data %>% data.frame
r4_2_df$treatment <- factor(r4_2_df$treatment)
r4_2_res1 <- aov(value ~ treatment, data = r4_2_df %>% filter(name == "Total wet weight"))
r4_2_res2 <- aov(value ~ treatment, data = r4_2_df %>% filter(name == "Total dry weight"))
r4_2_coefs1 <- summary(r4_2_res1)[[1]] %>% as.matrix
r4_2_coefs2 <- summary(r4_2_res2)[[1]] %>% as.matrix
r4_2_coefs1 <- c("total_wet_weight_12ind", r4_2_coefs1[1,c(4,5)])
r4_2_coefs2 <- c("total_dry_weight_12ind", r4_2_coefs2[1,c(4,5)])
names(r4_2_coefs1) <- rice_parms
names(r4_2_coefs2) <- rice_parms

# Grain count (/10 heads)
r4_3_df <- r4[[3]]$data %>% data.frame
r4_3_df$treatment <- factor(r4_3_df$treatment)
r4_3_res1 <- aov(value ~ treatment, data = r4_3_df %>% filter(name == "Total grain"))
r4_3_res2 <- aov(value ~ treatment, data = r4_3_df %>% filter(name == "Fertile grain"))
r4_3_res3 <- aov(value ~ treatment, data = r4_3_df %>% filter(name == "Sterile grain"))
r4_3_coefs1 <- summary(r4_3_res1)[[1]] %>% as.matrix
r4_3_coefs2 <- summary(r4_3_res2)[[1]] %>% as.matrix
r4_3_coefs3 <- summary(r4_3_res3)[[1]] %>% as.matrix
r4_3_coefs1 <- c("total_grain_10heads", r4_3_coefs1[1,c(4,5)])
r4_3_coefs2 <- c("fertile_grain_10heads", r4_3_coefs2[1,c(4,5)])
r4_3_coefs3 <- c("sterile_grain_10heads", r4_3_coefs3[1,c(4,5)])
names(r4_3_coefs1) <- rice_parms
names(r4_3_coefs2) <- rice_parms
names(r4_3_coefs3) <- rice_parms

# Grain dry weight (/10 heads)
r4_4_df <- r4[[4]]$data %>% data.frame
r4_4_df$treatment <- factor(r4_4_df$treatment)
r4_4_res1 <- aov(value ~ treatment, data = r4_4_df %>% filter(name == "Total grain"))
r4_4_res2 <- aov(value ~ treatment, data = r4_4_df %>% filter(name == "Fertile grain"))
r4_4_res3 <- aov(value ~ treatment, data = r4_4_df %>% filter(name == "Sterile grain"))
r4_4_coefs1 <- summary(r4_4_res1)[[1]] %>% as.matrix
r4_4_coefs2 <- summary(r4_4_res2)[[1]] %>% as.matrix
r4_4_coefs3 <- summary(r4_4_res3)[[1]] %>% as.matrix
r4_4_coefs1 <- c("total_dry_wt_10heads", r4_4_coefs1[1,c(4,5)])
r4_4_coefs2 <- c("fertile_dry_wt_10heads", r4_4_coefs2[1,c(4,5)])
r4_4_coefs3 <- c("sterile_dry_wt_10heads", r4_4_coefs3[1,c(4,5)])
names(r4_4_coefs1) <- rice_parms
names(r4_4_coefs2) <- rice_parms
names(r4_4_coefs3) <- rice_parms


# ------------------------------------------------ #
# Combine all the results (GAMM)
# ------------------------------------------------ #
stat_df <- data.frame(NULL)
stat_df <- stat_df %>% bind_rows(r1_1_coefs, r1_2_coefs, r1_3_coefs,
                                 r2_1_coefs, r2_2_coefs, r2_3_coefs,
                                 m3_1_coefs, m3_2_coefs, m3_3_coefs, m3_4_coefs, m3_5_coefs, m3_6_coefs, m3_7_coefs, m3_8_coefs,
                                 g1_1_coefs, g1_2_coefs, g1_3_coefs,
                                 i1_1_coefs, i1_2_coefs, i1_3_coefs,
                                 e1_1_coefs, e1_2_coefs)


# ------------------------------------------------ #
# Combine all the results (ANOVA)
# ------------------------------------------------ #
stat_yield_df <- data.frame(NULL)
stat_yield_df <- stat_yield_df %>% bind_rows(r4_1_coefs,
                                             r4_2_coefs1, r4_2_coefs2,
                                             r4_3_coefs1, r4_3_coefs2, r4_3_coefs3,
                                             r4_4_coefs1, r4_4_coefs2, r4_4_coefs3)


# ------------------------------------------------ #
# Save compiled data
# ------------------------------------------------ #
# Save csv
write.csv(stat_df, "data_tbl/gamm_result_allvars.csv", row.names = FALSE)
write.csv(stat_yield_df, "data_tbl/anova_result_riceyield.csv", row.names = FALSE)

# Save workspace and session information
macam::save_session_info()

