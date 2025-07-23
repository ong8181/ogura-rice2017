####
#### Plant cover data
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.10.31
library(patchwork); packageVersion("patchwork") # 1.3.0, 2024.10.31
library(cowplot); packageVersion("cowplot") # 1.1.3, 2024.10.31
library(ggforce); packageVersion("ggforce") # 0.4.2, 2024.11.01
library(ggtext); packageVersion("ggtext") # 0.1.2, 2025.07.15

# Create output folder
dir.create("00_SessionInfo")
(outdir <- macam::outdir_create())


# ------------------------------------- #
# Load data and compile data
# ------------------------------------- #
d <- read.csv("data/data_grass.csv") %>% filter(treatment != "NC")
# Edit class
d$date <- ymd(d$date)
d$treatment[d$treatment == "Conv"] <- "Conventional"
d$treatment[d$treatment == "NoFert"] <- "No Fertilizer"
# Remove dates when grass pics were not taken
d <- d %>%
  filter(grass_pic_taken == "Y") %>%
  replace(is.na(.), 0)

# Calculate mean
d_mean <- d %>%
  group_by(date, treatment) %>%
  summarize(cyper_1 = mean(cyper_1),
            cyper_2 = mean(cyper_2),
            cyper_3 = mean(cyper_3),
            cyper_4 = mean(cyper_4),
            poa_1 = mean(poa_1),
            poa_2 = mean(poa_2),
            poa_3 = mean(poa_3),
            poa_4 = mean(poa_4),
            poa_5 = mean(poa_5),
            portulaca_1 = mean(portulaca_1),
            euphorbia_1 = mean(euphorbia_1),
            brassica_1 = mean(brassica_1),
            aster_1 = mean(aster_1),
            aster_2 = mean(aster_2))
d_mean_long <- d_mean %>% 
  pivot_longer(cols = -c(date, treatment))

# Make data longer
d_grs <- d %>%
  select(-total_coverage, -grass_pic_taken) %>% 
  select(-sample_code, -n_week) %>% 
  pivot_longer(cols = -c(date, treatment, sampling_point))


# ------------------------------------- #
# Compile the plant names
# ------------------------------------- #
d_grs$name <- d_grs$name %>% str_replace("aster_1", "_Ixeris japonica_ (Aster.)")
d_grs$name <- d_grs$name %>% str_replace("aster_2", "_Taraxacum officinale_ (Aster.)")
d_grs$name <- d_grs$name %>% str_replace("brassica_1", "_Rorippa palustris_ (Brass.)")
d_grs$name <- d_grs$name %>% str_replace("cyper_1", "_Cyperus eragrostis_ (Cyper.)")
d_grs$name <- d_grs$name %>% str_replace("cyper_2", "_Fimbristylis diphylloides_ (Cyper.)")
d_grs$name <- d_grs$name %>% str_replace("cyper_3", "_Fimbristylis littoralis_ (Cyper.)")
d_grs$name <- d_grs$name %>% str_replace("cyper_4", "_Cyperus brevifolia_ (Cyper.)")
d_grs$name <- d_grs$name %>% str_replace("euphorbia_1", "_Euphorbia maculata_ (Euph.)")
d_grs$name <- d_grs$name %>% str_replace("poa_1", "_Digitaria ciliaris_ (Poac.)")
d_grs$name <- d_grs$name %>% str_replace("poa_2", "_Echinochloa crus-galli var. crus-galli_ (Poac.)")
d_grs$name <- d_grs$name %>% str_replace("poa_3", "_Dinebra chinensis_ (Poac.)")
d_grs$name <- d_grs$name %>% str_replace("poa_4", "_Paspalum thunbergii_ (Poac.)")
d_grs$name <- d_grs$name %>% str_replace("poa_5", "_Paspalum dilatatum_ (Poac.)")
d_grs$name <- d_grs$name %>% str_replace("portulaca_1", "_Portulaca oleracea_ (Port.)")
d_grs$name <- factor(d_grs$name,
                     levels = c("_Ixeris japonica_ (Aster.)", "_Taraxacum officinale_ (Aster.)", "_Rorippa palustris_ (Brass.)",
                                "_Cyperus eragrostis_ (Cyper.)", "_Cyperus brevifolia_ (Cyper.)", "_Fimbristylis diphylloides_ (Cyper.)", "_Fimbristylis littoralis_ (Cyper.)",
                                "_Euphorbia maculata_ (Euph.)", "_Digitaria ciliaris_ (Poac.)", "_Echinochloa crus-galli var. crus-galli_ (Poac.)",
                                "_Dinebra chinensis_ (Poac.)", "_Paspalum thunbergii_ (Poac.)", "_Paspalum dilatatum_ (Poac.)", "_Portulaca oleracea_ (Port.)"))

# ------------------------------------- #
# Visualize data (jitter + line)
# ------------------------------------- #
# Plant cover
g1 <- d_grs %>%
  ggplot(aes(x = date, y = value, color = treatment, group = sampling_point, shape = treatment)) +
  geom_line() +
  geom_jitter(height = 0, width = 0.5, alpha = 0.5) +
  facet_wrap(~ name, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  xlab(NULL) + ggtitle("Plant cover") +
  ylab("Plant coverage (%)") +
  theme(strip.text = element_markdown())
g1

# ------------------------------------- #
# Merge and compile data
# ------------------------------------- #
# Merge and compile
d_merge <- d %>%
  mutate(cyper_all = rowSums(select(., starts_with("cyper_"))),
         poa_all =  rowSums(select(., starts_with("poa_"))),
         others_all = rowSums(select(., portulaca_1, euphorbia_1, brassica_1, starts_with("aster_")))) %>% 
  select(date, n_week, treatment, sampling_point, cyper_all, poa_all, others_all)
d_merge_grs <- d_merge %>% select(-n_week) %>% 
  pivot_longer(cols = -c(date, treatment, sampling_point))
d_merge_mean <- d_merge %>%
  group_by(date, treatment) %>%
  summarize(cyper_all = mean(cyper_all),
            poa_all = mean(poa_all),
            others_all = mean(others_all))
d_merge_mean_long <- d_merge_mean %>% 
  pivot_longer(cols = -c(date, treatment))

g2 <- d_merge_grs %>%
  ggplot(aes(x = date, y = value, color = treatment, shape = treatment)) +
  geom_line(data = d_merge_mean_long, aes(x = date, y = value, color = treatment)) +
  geom_jitter(height = 0, width = 0.5, alpha = 0.5) +
  facet_wrap(~ name, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  xlab(NULL) + ggtitle("Plant cover") +
  ylab("Plant coverage (%)")


# ------------------------------------- #
# Combine all the figures
# ------------------------------------- #
ggsave("01_VisualizeGrassOut/Fig_GrassCov.pdf", g1, height = 12, width = 10)
ggsave("01_VisualizeGrassOut/Fig_GrassCov_merge.pdf", g2, height = 7, width = 10)
saveRDS(list(g1,g2), "../07_FormatFigs/data_robj/Grass_temporal-pattern.obj")


# ------------------------------------- #
# Save workspace
# ------------------------------------- #
macam::save_session_info("00_SessionInfo")
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))


