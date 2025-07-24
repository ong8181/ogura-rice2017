####
#### Ogura Rice 2017
#### Compile and format figures
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.08.03
library(patchwork); packageVersion("patchwork") # 1.3.0, 2024.12.04
library(cowplot); packageVersion("cowplot") # 1.1.3, 2024.12.10
library(ggforce); packageVersion("ggforce") # 0.4.2, 2024.12.10
library(ggtext); packageVersion("ggtext") # 0.1.2, 2024.12.13

# Create output directory
dir.create("00_SessionInfo")


# ------------------------------------------------ #
# Load all figure data
# ------------------------------------------------ #
# Climate
c1 <- readRDS("data_robj/Climate_all.obj")
# Rice
r1 <- readRDS("data_robj/Rice_Growth_rawdata.obj")
r2 <- readRDS("data_robj/Rice_Growth_diff.obj")
r3 <- readRDS("data_robj/Rice_Yield.obj")
# Specialized metabolites
m1 <- readRDS("data_robj/SpecializedMetabolites.obj")
m2 <- readRDS("data_robj/SpecializedMetabolites_DimRed.obj")
# Grass
g1 <- readRDS("data_robj/Grass_temporal-pattern.obj")
g2 <- readRDS("data_robj/Grass_DimRed.obj")
# Insect
i1 <- readRDS("data_robj/Insect_temporal-pattern.obj")
i2 <- readRDS("data_robj/Insect_DimRed.obj")
# eDNA
e1 <- readRDS("data_robj/eDNA_temporal-pattern.obj")
e2 <- readRDS("data_robj/eDNA_DimRed.obj")
e3 <- readRDS("data_robj/eDNA_top_temporal-pattern.obj")
e4 <- readRDS("data_robj/eDNA_TopTaxa.obj")
e5 <- readRDS("data_robj/eDNA_major_groups.obj")
e6 <- readRDS("data_robj/eDNA_caual_taxa.obj")
e7 <- readRDS("data_robj/eDNA_relative.obj")
es1 <- readRDS("data_robj/SeqReadsOverview_16S.obj")
es2 <- readRDS("data_robj/SeqReadsOverview_18S.obj")
es3 <- readRDS("data_robj/SeqReadsOverview_COI.obj")
es4 <- readRDS("data_robj/SeqReadsOverview_ITS.obj")
# MDR S-map results
s1 <- readRDS("data_robj/Smap_Height.obj")
s2 <- readRDS("data_robj/Smap_SPAD.obj")
s3 <- readRDS("data_robj/Smap_Nstem.obj")
# MDR S-map results (all)
s1_all <- readRDS("data_robj/Smap_Height-all.obj")
s2_all <- readRDS("data_robj/Smap_SPAD-all.obj")
s3_all <- readRDS("data_robj/Smap_Nstem-all.obj")


# ------------------------------------------------ #
# Compile figures
# ------------------------------------------------ #
# Set figure theme
theme_set(theme_bw())
col_conv <- "darkgreen"
col_nfrt <- "darkgoldenrod2"
p_alpha <- 0.6
p_size <- 2


# ------------------------------------------------ #
# Figure: Climate
# ------------------------------------------------ #
Fig_climate <- (c1[[1]] + theme(plot.tag = element_text(face = "bold")) +
                  labs(y = "Mean daily\ntemperature (\u00B0C)")) /
  (c1[[2]] + theme(plot.tag = element_text(face = "bold"))) /
  (c1[[3]] + theme(plot.tag = element_text(face = "bold"))) +
  plot_annotation(tag_levels = list(c("b","c","d")))


# ------------------------------------------------ #
# Figure: Rice growth pattern
# ------------------------------------------------ #
## Raw data
r1_1 <- r1[[1]] + theme(legend.position = "top", plot.tag = element_text(face = "bold")) +
  scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
  annotate("text", x = ymd("2017-06-01"), y = 110, label = "Management: ****")
r1_2 <- r1[[2]] + theme(legend.position = "none", plot.tag = element_text(face = "bold")) +
  scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
  annotate("text", x = ymd("2017-06-01"), y = 42, label = "Management: ****")
r1_3 <- r1[[3]] + theme(legend.position = "none", plot.tag = element_text(face = "bold")) +
  scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
  annotate("text", x = ymd("2017-06-01"), y = 31, label = "Management: ****")
r1_1$layers[[2]]$aes_params$alpha <- r1_2$layers[[2]]$aes_params$alpha <- r1_3$layers[[2]]$aes_params$alpha <- p_alpha
r1_1$layers[[2]]$aes_params$size <- r1_2$layers[[2]]$aes_params$size <- r1_3$layers[[2]]$aes_params$size <- p_size
## Diff data
r2_1 <- r2[[1]] + theme(legend.position = "top", plot.tag = element_text(face = "bold")) +
  scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
  annotate("text", x = ymd("2017-08-10"), y = 19, label = "Management: ****")
r2_2 <- r2[[2]] + theme(legend.position = "none", plot.tag = element_text(face = "bold")) +
  scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
  annotate("text", x = ymd("2017-08-10"), y = 17, label = "Management: N.S.")
r2_3 <- r2[[3]] + theme(legend.position = "none", plot.tag = element_text(face = "bold")) +
  scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
  annotate("text", x = ymd("2017-08-10"), y = 8.5, label = "Management: ****")
r2_1$layers[[2]]$aes_params$alpha <- r2_2$layers[[2]]$aes_params$alpha <- r2_3$layers[[2]]$aes_params$alpha <- p_alpha
r2_1$layers[[2]]$aes_params$size <- r2_2$layers[[2]]$aes_params$size <- r2_3$layers[[2]]$aes_params$size <- p_size
## Collect figures
Fig_rice_raw <- r1_1 / r1_2 / r1_3 + plot_annotation(tag_levels = list(c("a","b","c")))
Fig_rice_dif <- r2_1 / r2_2 / r2_3 + plot_annotation(tag_levels = list(c("d","e","f")))


# Rice yield
## Annotate text
y1_label <- data.frame(name = c("Total grain", "Fertile grain", "Sterile grain"),
                       treatment = c("Conventional", "Conventional", "Conventional"),
                       x = c(1.5, 1.5, 1.5),
                       y = c(1200, 1200, 700),
                       label = c("**", "**", "N.S."))
y1_label$name <- factor(y1_label$name, levels = c("Total grain", "Fertile grain", "Sterile grain"))
y2_label <- data.frame(name = c("Total grain", "Fertile grain", "Sterile grain"),
                       treatment = c("Conventional", "Conventional", "Conventional"),
                       x = c(1.5, 1.5, 1.5),
                       y = c(29, 29, 10),
                       label = c("**", "**", "N.S."))
y2_label$name <- factor(y2_label$name, levels = c("Total grain", "Fertile grain", "Sterile grain"))
y3_label <- data.frame(name = c("Total wet weight", "Total dry weight"),
                       treatment = c("Conventional", "Conventional"),
                       x = c(1.5, 1.5),
                       y = c(710, 710),
                       label = c("****", "****"))
y3_label$name <- factor(y3_label$name, levels = c("Total wet weight", "Total dry weight"))

## Edit figures
y1 <- r3[[2]] +
  theme(plot.tag = element_text(face = "bold"), legend.position = "top") + ylim(0,1200) +
  scale_x_discrete(labels=c("Conventional", "No Fertilizer")) +
  geom_text(data = y1_label, aes(x = x, y = y, label = label), color = "black")
y2 <- r3[[3]] + theme(plot.tag = element_text(face = "bold"), legend.position = "none") + ylim(0,30) +
  scale_x_discrete(labels=c("Conventional", "No Fertilizer")) +
  geom_text(data = y2_label, aes(x = x, y = y, label = label), color = "black")
y3 <- (r3[[1]][[1]] + theme(plot.tag = element_text(face = "bold"), legend.position = "top") + ylim(0,330)+
         scale_x_discrete(labels=c("Conventional", "No Fertilizer")) + ylab("Total number of rice heads (/12 inds.)") +
         annotate("text", x = 1.5, y = 320, label = "****")) /
  (r3[[1]][[2]] + theme(plot.tag = element_text(face = "bold"), legend.position = "none") + ylim(0,720) +
     facet_wrap(~ name, ncol = 1) + scale_x_discrete(labels=c("Conventional", "No Fertilizer")) +
     geom_text(data = y3_label, aes(x = x, y = y, label = label), color = "black") +
     ylab("Head weight (g/12 inds.)")) +
  plot_annotation(tag_levels = list(c("g","h"))) + plot_layout(heights = c(1,2))
Fig_rice_yld <- y1 / y2 + plot_annotation(tag_levels = "a")

# Combine rice growth and yield
Fig_rice_all <- plot_grid(Fig_rice_raw, Fig_rice_dif, y3, nrow = 1, rel_widths = c(2,2,1.2))


# ------------------------------------------------ #
# Figure: Ecological variables
# ------------------------------------------------ #
## Specialized metabolites
### Prepare new labels
m_labels <- c("Z_3_Hexenal"="(_Z_)-3-Hexenal", "E_2_Hexenal"="(_E_)-2-Hexenal",
              "Z_3_Hexen_1_ol"="(_Z_)-3-Hexen-1-ol", "X2_4_Heptadienal"="2-4-Heptadienal",
              "linalool"="Linalool", "b_caryophyllene"="&beta;-caryophyllene",
              "MeSA"="Methyl-salicylate (MeSA)", "diterpenoids_ug_g_FW"="Diterpenoids (&mu;g/g)")
## Annotate text
m1_names <- c("Z_3_Hexenal","E_2_Hexenal","Z_3_Hexen_1_ol","X2_4_Heptadienal","linalool", "b_caryophyllene","MeSA","diterpenoids_ug_g_FW")
m1_label <- data.frame(name = m1_names,
                       treatment = rep("Conventional",8), x = rep(ymd("2017-07-15"),8),
                       y = c(300, 600, 16, 10, 7.5, 0.7, 3.9, 16),
                       label = c("****", "**", "N.S.", "N.S.", "N.S.", "N.S.", "**", "N.S."))
m1_label$name <- factor(m1_label$name, levels = m1_names)
## Edit figure
Fig_mtbl_raw <- m1[[3]] +
  scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
  ylab("Normalized concentration\n(relatie peak area/leaf mass)") +
  facet_wrap(~ name, labeller = as_labeller(m_labels), nrow = 2, scales = "free_y") +
  theme(legend.position = "none", strip.text = element_markdown()) +
  geom_text(data = m1_label, aes(x = x, y = y, label = label), color = "black") +
  NULL
Fig_mtbl_raw$layers[[2]]$aes_params$alpha <- p_alpha
Fig_mtbl_raw$layers[[2]]$aes_params$size <- p_size

## Grass cover
grs_levels <- c("cyper_all", "poa_all", "others_all")
g1_label <- data.frame(name = grs_levels,
                       treatment = rep("Conventional",3),
                       x = c(ymd("2017-06-15"),ymd("2017-07-10"),ymd("2017-07-10")),
                       y = c(60, 75, 75),
                       label = c("**", "N.S.", "*"))
g1_label$name <- factor(g1_label$name, levels = grs_levels)

g1[[2]]$data$name <- factor(g1[[2]]$data$name, levels = grs_levels)
g1[[2]]$layers[[1]]$data$name <- factor(g1[[2]]$layers[[1]]$data$name, levels = grs_levels)
g_labels <- c("cyper_all"="Cyperaceae", "poa_all"="Poaceae", "others_all"="Others")
Fig_grss_sum <- g1[[2]] + scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
  facet_wrap(~ name, labeller = as_labeller(g_labels), nrow = 1, scales = "free_y") +
  geom_text(data = g1_label, aes(x = x, y = y, label = label), color = "black") +
  theme(legend.position = "bottom") 
Fig_grss_sum$layers[[2]]$aes_params$alpha <- p_alpha
Fig_grss_sum$layers[[2]]$aes_params$size <- p_size
Fig_grss_all <- g1[[1]] + theme(legend.position = "bottom") + scale_color_manual(values = c(col_conv, col_nfrt), name = NULL)

## Insect
insc_levels <- c("cicadellidae_spp","thysanoptera_spp","other_planthopper")
i1_label <- data.frame(name = insc_levels,
                       treatment = rep("Conventional",3),
                       x = rep(ymd("2017-08-05"), 3),
                       y = c(50, 130, 130),
                       label = c("N.S.", "N.S.", "***"))
i1_label$name <- factor(i1_label$name, levels = insc_levels)

i1$data$name <- factor(i1$data$name, levels = insc_levels)
i1$layers[[1]]$data$name <- factor(i1$layers[[1]]$data$name, levels = insc_levels)
i_labels <- c("cicadellidae_spp"="Cicadellidae","thysanoptera_spp"="Thysanoptera","other_planthopper"="Other planthoppers")
Fig_insc_raw <- i1 + scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
  facet_wrap(~ name, labeller = as_labeller(i_labels), nrow = 1, scales = "free_y")  +
  geom_text(data = i1_label, aes(x = x, y = y, label = label), color = "black") +
  theme(legend.position = "none") 
Fig_insc_raw$layers[[2]]$aes_params$alpha <- p_alpha
Fig_insc_raw$layers[[2]]$aes_params$size <- p_size

## Combine metabolites, grass, and insect figures
Fig_mgi_all <- (Fig_mtbl_raw + theme(plot.tag = element_text(face = "bold"))) / 
  (Fig_insc_raw + theme(plot.tag = element_text(face = "bold"))) / 
  (Fig_grss_sum + theme(plot.tag = element_text(face = "bold"))) +
  plot_layout(heights = c(2,1,1)) +
  plot_annotation(tag_levels = "a") 

## Dimensional reduction for plant cover and metabolites
g2_1 <- g2[[2]] +
  scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
  scale_fill_manual(values = c(col_conv, col_nfrt), name = NULL) +
  theme(plot.tag = element_text(face = "bold")) +
  labs(title = "Plant coverage (%) (stress = 0.1899)",
       subtitle = "The number of weeks, farming practice, and their interactions had statistically clear effects") +
  NULL
m2_1 <- m2[[1]] +
  scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
  scale_fill_manual(values = c(col_conv, col_nfrt), name = NULL) +
  theme(plot.tag = element_text(face = "bold")) +
  labs(title = "Plant specialized metabolites (stress = 0.1568)",
       subtitle = "The number of weeks had statistically clear effects") +
  NULL
Fig_dim1 <- m2_1 / g2_1 + plot_annotation(tag_levels = "a")


# ------------------------------------------------ #
## eDNA
# ------------------------------------------------ #
### Overall patterns
e1_1 <- e1[[2]] + scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
  annotate("text", x = ymd("2017-07-15"), y = 120, label = "Management: ****")
e1_1$layers[[1]]$aes_params$alpha <- p_alpha
e1_1$layers[[1]]$aes_params$size <- p_size
#### NMDS
e2_1 <- e2[[1]] +
  scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
  scale_fill_manual(values = c(col_conv, col_nfrt), name = NULL) +
  labs(title = "eDNA copy numbers (stress = 0.1924)",
       subtitle = "The number of weeks, farming practice, and their interactions had statistically clear effects") +
  NULL
#### Combine all
e1_label <- data.frame(treatment = c("Conventional", "No Fertilizer"),
                       rep_tax = c(NA, NA),
                       x = c(ymd("2017-08-05"), ymd("2017-05-20")),
                       y = c(7e+07, 6e+07),
                       label = c("Management:****", NA))
Fig_edna <- (e1[[1]] +
               geom_text(data = e1_label, aes(x = x, y = y, label = label), color = "black") +
               theme(plot.tag = element_text(face = "bold"))) /
  (e1_1 + ggtitle("OTU diversity (Management: ***)") + ylab("Number of OTUs") +
     theme(plot.tag = element_text(face = "bold"))) /
  (e2_1 + theme(plot.tag = element_text(face = "bold"))) +
  plot_layout(heights = c(1.1, 1, 1)) +
  plot_annotation(tag_levels = "a") 

# eDNA supplementary figures
e3_labels <- c("PRO_Taxa00006"="PRO_Taxa00006 (Comamonadaceae)", "PRO_Taxa00008"="PRO_Taxa00008 (Eukaryota)",
               "PRO_Taxa00010"="PRO_Taxa00010 (_Flavobacterium_)", "PRO_Taxa00009"="PRO_Taxa00009 (Bacteroidota)",
               "PRO_Taxa00011"="PRO_Taxa00011 (Bacteria)", "PRO_Taxa00013"="PRO_Taxa00013 (Actinomycetes)",
               "PRO_Taxa00004"="PRO_Taxa00004 (_Paenibacillus_)", "PRO_Taxa00007"="PRO_Taxa00007 (Actinomycetes)",
               "PRO_Taxa00014"="PRO_Taxa00014 (_Flavobacterium_)", "PRO_Taxa00027"="PRO_Taxa00027 (Bacteria)",
               "PRO_Taxa00019"="PRO_Taxa00019 (Betaproteobacteria)", "PRO_Taxa00023"="PRO_Taxa00023 (_Pedobacter_)")
eDNAtop <- e3 + scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
  scale_fill_manual(values = c(col_conv, col_nfrt), name = NULL) +
  facet_wrap(~ OTU, scales = "free_y", labeller = as_labeller(e3_labels), nrow = 3) +
  theme(legend.position = "bottom", strip.text = element_markdown(), axis.title.y = element_markdown(), plot.tag = element_text(face = "bold")) +
  labs(title = "Top eDNA taxa",
       y = "Abundance (eDNA copies/ml water + 0.5)") +
  NULL

e4_labels <- c("PRO_Taxa00038"="PRO_Taxa00038 (Actinomycetes)", "PRO_Taxa00061"="PRO_Taxa00061 (Gammaproteobacteria)",
               "EUK_Taxa00281"="EUK_Taxa00281 (Chlorophyceae)", "PRO_Taxa00262"="PRO_Taxa00262 (Gallionellaceae)",
               "PRO_Taxa00214"="PRO_Taxa00214 (Verrucomicrobiota)", "EUK_Taxa00423"="EUK_Taxa00423 (Chlorophyceae)",
               "EUK_Taxa00211"="EUK_Taxa00211 (Ulotrichales)", "PRO_Taxa00195"="PRO_Taxa00195 (Silvanigrellaceae)",
               "PRO_Taxa00066"="PRO_Taxa00066 (_Pedobacter_)", "EUK_Taxa00234"="EUK_Taxa00234 (Chlamydomonadaceae)",
               "PRO_Taxa00049"="PRO_Taxa00049 (Bacteria)", "EUK_Taxa00170"="EUK_Taxa00170 (Chrysophyceae)")
Fig_eDNAtopSI1 <- eDNAtop

eDNAtopdif <- e4 + scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
  scale_fill_manual(values = c(col_conv, col_nfrt), name = NULL) +
  facet_wrap(~ OTU, scales = "free_y", labeller = as_labeller(e4_labels), nrow = 3) +
  labs(title = "Top eDNA taxa contributing to the difference between the two paddy fields",
       y = "Abundance (eDNA copies/ml water + 0.5)") +
  theme(legend.position = "bottom", strip.text = element_markdown(), axis.title.y = element_markdown(), plot.tag = element_text(face = "bold")) +
  NULL

# Major microbial groups (to compare them with Kamata et al. 1991)
Fig_eDNA_bac <- (e5[[2]] +
                   scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
                   scale_fill_manual(values = c(col_conv, col_nfrt), name = NULL) +
                   theme(legend.position = "none",
                         plot.tag = element_text(face = "bold")))
Fig_eDNA_fun <- (e5[[1]] +
                   scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
                   scale_fill_manual(values = c(col_conv, col_nfrt), name = NULL) +
                   theme(legend.position = "bottom",
                         plot.tag = element_text(face = "bold")))
Fig_eDNA_rel1 <- (e7[[1]] + theme(legend.position = "right", 
                                  plot.tag = element_text(face = "bold")) +
                    guides(fill = guide_legend(ncol = 2)))
Fig_eDNA_rel2 <- (e7[[2]] + theme(legend.position = "right",
                                  plot.tag = element_text(face = "bold")))
Fig_eDNAtopSI2 <- 
  (Fig_eDNA_rel1 + plot_spacer() + plot_layout(widths = c(1,0.5))) / Fig_eDNA_bac / 
  (Fig_eDNA_rel2 + plot_spacer() + plot_layout(widths = c(1,0.7))) / Fig_eDNA_fun +
  plot_layout(heights = c(0.35, 1, 0.35, 2/3)) +
  plot_annotation(tag_levels = "a") 

Fig_eDNA_rel3 <- (e7[[3]] + theme(legend.position = "right", 
                                  plot.tag = element_text(face = "bold")) +
                    guides(fill = guide_legend(ncol = 4)))
Fig_eDNA_euk <- (e5[[3]] +
                   scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
                   scale_fill_manual(values = c(col_conv, col_nfrt), name = NULL) +
                   theme(legend.position = "bottom",
                         plot.tag = element_text(face = "bold")))
Fig_eDNAtopSI3 <- (Fig_eDNA_rel3 + plot_spacer() + plot_layout(widths = c(1,0.2))) / Fig_eDNA_euk +
  plot_layout(heights = c(0.2, 1)) +
  plot_annotation(tag_levels = "a") 

# Major groups differentiate the two paddy fields
Fig_eDNAtopSI4 <- eDNAtopdif


# ------------------------------------------------ #
# eDNA causal taxa (by UIC and S-map)
# ------------------------------------------------ #
Fig_eDNA_cause1 <- (e6[[1]] +
                      scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
                      scale_fill_manual(values = c(col_conv, col_nfrt), name = NULL) +
                      theme(legend.position = "none",
                            plot.tag = element_text(face = "bold")))
Fig_eDNA_cause2 <- (e6[[2]] +
                      scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
                      scale_fill_manual(values = c(col_conv, col_nfrt), name = NULL) +
                      theme(legend.position = "none",
                            plot.tag = element_text(face = "bold")))
Fig_eDNA_cause3 <- (e6[[3]] +
                      scale_color_manual(values = c(col_conv, col_nfrt), name = NULL) +
                      scale_fill_manual(values = c(col_conv, col_nfrt), name = NULL) +
                      theme(legend.position = "bottom",
                            plot.tag = element_text(face = "bold")))
Fig_eDNA_Causal <- 
  Fig_eDNA_cause1 / Fig_eDNA_cause2 / (Fig_eDNA_cause3 + plot_spacer() + plot_layout(widths = c(1,1.07))) +
  plot_layout(heights = c(1, 2/4, 1/4)) +
  plot_annotation(tag_levels = "a") 


# ------------------------------------------------ #
# Sequence summary (SI)
# ------------------------------------------------ #
SeqRead_16S <- es1 + ggtitle("16S (515F/806R)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none", plot.tag = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))
SeqRead_18S <- es2 + ggtitle("18S (Euk_1391f/EukBr)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "right", plot.tag = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))
SeqRead_COI <- es3 + ggtitle("COI (mlCOIintF/HCO2198)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none", plot.tag = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))
SeqRead_ITS <- es4 + ggtitle("ITS (ITS1-F-KYO1/ITS2-KYO2)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "right", plot.tag = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))
FigSeqRead <- (SeqRead_16S + SeqRead_18S)/(SeqRead_COI + SeqRead_ITS) +
  plot_annotation(tag_levels = "a")


# ------------------------------------------------ #
# Figure: UIC and S-map analysis
# ------------------------------------------------ #
tag_position1 <- c(0.48,1.025)
tag_position2 <- c(-0.03,1.025)
top_margin <- 20

smap_label1 <- c("ITS_Taxa00138_tp0"="Undetermined (ITS_Taxa00138, tp=0)",
                 "PRO_Taxa00019_tp-2"="Betaproteobacteria (PRO_Taxa00019, tp=-2)",
                 "PRO_Taxa00332_tp-2"="Bacteria (PRO_Taxa00332, tp=-2)",
                 "ITS_Taxa00007_tp-2"="Undetermined (ITS_Taxa00007, tp=-2)",
                 "PRO_Taxa00168_tp0"="Bacteria (PRO_Taxa00168, tp=0)", 
                 "PRO_Taxa00186_tp-2"="Pseudomonadota (PRO_Taxa00186, tp=-2)",
                 "PRO_Taxa00037_tp-2"="_Flavobacterium_ (PRO_Taxa00037, tp=-2)", 
                 "PRO_Taxa00024_tp-2"="_Polynucleobacter cosmopolitanus_ (PRO_Taxa00024, tp=-2)",
                 "PRO_Taxa00040_tp-1"="Bacteria (PRO_Taxa00040, tp=-1)",
                 "PRO_Taxa00096_tp0"="Bacteria (PRO_Taxa00096, tp=0)",
                 "COI_Taxa00221_tp0"="_Flavobacterium_ (COI_Taxa00221, tp=0)",
                 "grss_poa_1_tp-1"="_Digitaria ciliaris_ (Poac.) (Plant on the ridge, tp=-1)",
                 "EUK_Taxa00213_tp0"="Undetermined (EUK_Taxa00213, tp=0)",
                 "PRO_Taxa00046_tp-2"="Bacteroidota (PRO_Taxa00046, tp=-2)",
                 "EUK_Taxa00443_tp-2"="Chlamydomonadales (EUK_Taxa00443, tp=-2)",
                 "PRO_Taxa00017_tp-1"="Burkholderiales (PRO_Taxa00017, tp=-1)")
s1_all2 <-
  (s1[[1]] + scale_y_discrete(labels = smap_label1) +
     labs(tag = "a") +
     theme(legend.position = "none",
           axis.text.y = element_markdown(),
           plot.margin = margin(t=top_margin),
           plot.tag.position = tag_position1, plot.tag = element_text(face = "bold"))) +
  (s1[[2]] + labs(tag = "b") +
     theme(axis.text.y = element_blank(),
           plot.tag.position = tag_position2,
           plot.margin = margin(t=top_margin),
           legend.position = "none", plot.tag = element_text(face = "bold"))) +
  (s1[[3]] + labs(tag = "c") +
     theme(plot.tag.position = tag_position2,
           plot.margin = margin(t=top_margin),
           axis.text.y = element_blank(), plot.tag = element_text(face = "bold")))
Fig_smap_hgt <- s1_all2

# Supplementary figures (S-map SPAD)
smap_label2 <- c("COI_Taxa00127_tp0"="_Cochliopodium_ (COI_Taxa00127, tp=0)",
                 "EUK_Taxa00659_tp-1"="Chlorophyceae (EUK_Taxa00659, tp=-1)",
                 "PRO_Taxa00050_tp-2"="Bacteria (PRO_Taxa00050, tp=-2)",
                 "COI_Taxa00134_tp-1"="Ploima (COI_Taxa00134, tp=-1)",
                 "EUK_Taxa00306_tp-2"="Eukaryota (EUK_Taxa00306, tp=-2)", 
                 "PRO_Taxa00013_tp0"="Actinomycetes (PRO_Taxa00013, tp=0)",
                 "PRO_Taxa00011_tp0"="Bacteria (PRO_Taxa00011, tp=0)")
s2_all2 <-
  (s2[[1]] + scale_y_discrete(labels = smap_label2) +
     labs(tag = "a") +
     theme(legend.position = "none",
           axis.text.y = element_markdown(),
           plot.margin = margin(t=top_margin),
           plot.tag.position = c(0.35, 1.03), plot.tag = element_text(face = "bold"))) +
  (s2[[2]] + labs(tag = "b") +
     theme(axis.text.y = element_blank(),
           plot.tag.position = tag_position2,
           plot.margin = margin(t=top_margin),
           legend.position = "none", plot.tag = element_text(face = "bold"))) +
  (s2[[3]] + labs(tag = "c") +
     theme(plot.tag.position = tag_position2,
           plot.margin = margin(t=top_margin),
           axis.text.y = element_blank(), plot.tag = element_text(face = "bold")))
Fig_smap_spd <- s2_all2

# Supplementary figures (S-map SPAD)
smap_label3 <- c("PRO_Taxa00195_tp-1"="Silvanigrellaceae (PRO_Taxa00195, tp=-1)",
                 "EUK_Taxa00423_tp-1"="Chlorophyceae (EUK_Taxa00423, tp=-1)")
s3_all2 <-
  (s3[[1]] + scale_y_discrete(labels = smap_label3) +
     labs(tag = "a") +
     theme(legend.position = "none",
           axis.text.y = element_markdown(),
           plot.margin = margin(t=top_margin),
           plot.tag.position = c(0.38, 1.03), plot.tag = element_text(face = "bold"))) +
  (s3[[2]] + labs(tag = "b") +
     theme(axis.text.y = element_blank(),
           plot.tag.position = tag_position2,
           plot.margin = margin(t=top_margin),
           legend.position = "none", plot.tag = element_text(face = "bold"))) +
  (s3[[3]] + labs(tag = "c") +
     theme(plot.tag.position = tag_position2,
           plot.margin = margin(t=top_margin),
           axis.text.y = element_blank(), plot.tag = element_text(face = "bold")))
Fig_smap_stm <- s3_all2



# ------------------------------------------------ #
# Save compiled figures
# ------------------------------------------------ #
# Main figures
ggsave("formatted_figs/Fig1_Climate.pdf", Fig_climate, width = 6, height = 6)
ggsave("formatted_figs/Fig2_RiceGrowth.pdf", Fig_rice_all, width = 14, height = 10)
#ggsave("formatted_figs/Fig3_MetbGrInsc.pdf", Fig_mgi_all, width = 12, height = 10)
Cairo::CairoPDF(file = "formatted_figs/Fig3_MetbGrInsc.pdf", width = 12, height = 10); print(Fig_mgi_all); dev.off() 
ggsave("formatted_figs/Fig3_MetbGrInsc.jpg", Fig_mgi_all, width = 12, height = 10)
ggsave("formatted_figs/Fig4_NMDSMetbGr.pdf", Fig_dim1, width = 10, height = 9)
ggsave("formatted_figs/Fig5_eDNApattern.pdf", Fig_edna, width = 10, height = 12)
ggsave("formatted_figs/Fig5_eDNApattern.jpg", Fig_edna, width = 10, height = 12)
ggsave("formatted_figs/Fig6_SmapHeight.pdf", Fig_smap_hgt, width = 14, height = 9)

# Supplementary figures
ggsave("formatted_figs/FigS02_RiceGrains.pdf", Fig_rice_yld, width = 7, height = 9)
ggsave("formatted_figs/FigS03_GrassAll.pdf", Fig_grss_all, width = 10, height = 10)
ggsave("formatted_figs/FigS04_SeqReads.pdf", FigSeqRead, width = 12, height = 10)
ggsave("formatted_figs/FigS05_TopeDNA.pdf", Fig_eDNAtopSI1, width = 14, height = 9)
ggsave("formatted_figs/FigS06_MajorMicGroups1.pdf", Fig_eDNAtopSI2, width = 14, height = 20)
ggsave("formatted_figs/FigS07_MajorMicGroups2.pdf", Fig_eDNAtopSI3, width = 14, height = 16)
ggsave("formatted_figs/FigS08_TopDifeDNA.pdf", Fig_eDNAtopSI4, width = 14, height = 9)
ggsave("formatted_figs/FigS09_SmapSPAD.pdf", Fig_smap_spd, width = 14, height = 7)
ggsave("formatted_figs/FigS10_SmapNstem.pdf", Fig_smap_stm, width = 14, height = 3)
ggsave("formatted_figs/FigS11_CausalSppDynamics.pdf", Fig_eDNA_Causal, width = 14, height = 16)


# ------------------------------------------------ #
# Save compiled data
# ------------------------------------------------ #
# Save workspace and session information
macam::save_session_info()

