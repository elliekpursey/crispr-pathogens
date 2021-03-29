library(tidyverse)
library(ggeffects) 
library(ghibli)
library(MuMIn)
library(DHARMa)
library(ggpubr)
library(cowplot)

setwd("results")

# read in dataframe
all_df <- read_csv("all_df.csv")

# remove prediction for binomial models, so only one binomial entry per genome

model_df <- all_df %>% 
  dplyr::select(-Prediction, -taxid, -n_crisprcas_systems, -total_spacers) %>%
  group_by(id, amr_binomial, length, crisprcas_binomial, species, n_plasmid_replicons, n_ICEs, n_amr_genes, n_intI1) %>%
  summarise(across(everything(), sum))
  
######## CRISPR-Cas binomial models ########

# model1 - Is genome length predictive of crispr-cas presence/absence?

# density plot of length - fairly normally distributed within species although not overall
ggplot(model_df, aes(x=length)) + 
  geom_histogram() +
  facet_grid(~species)

# linear model with interaction term
lm1 <- lm(length ~ crisprcas_binomial + species + crisprcas_binomial*species, data=model_df, na.action = "na.fail")
lm_summary <- summary(lm1)

# model selection
lm_dredge <- dredge(lm1)
write.csv(lm_dredge, "length_AICs.csv")

lm_tidy <- lm_summary$coefficients
write.csv(lm_tidy, file = "length_effects.csv")

# diagnostics - qqplot looks bad but using model anyway as only predicting per species, so overall distribution less important
plot(lm1)

# plot of residuals - they are normally distributed
hist(residuals(lm1), breaks = 100)

# model2 - What best predicts CRISPR-Cas presence/absence?

# remove NG as it never has CRISPR-Cas
glm_df <- model_df %>%
  filter(species != "NG")

# glm
glm <- glm(crisprcas_binomial ~ species*n_amr_genes + species*n_plasmid_replicons 
           + species*n_ICEs + species*n_intI1, 
           data=glm_df, family="binomial",  na.action = "na.fail")


glm_dredge <- dredge(glm)
write.csv(glm_dredge, "binom_AICs.csv")

glm_summary <- summary(glm)
glm_tidy <- glm_summary$coefficients
write.csv(glm_tidy, file = "binom_effects.csv")

# test fit of model 
testDispersion(glm)
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)

# prediction plot - lengths vs crispr

# create prediction dataframe, with modified + and - labels, removing NG (485)
pred_lm1 <- ggpredict(lm1, terms = c("crisprcas_binomial", "species"), type='fixed', ci.lvl = 0.95) %>%
  mutate(x, x = ifelse(x == 1, "+", "-")) %>%
  filter(group != "NG") %>%
  mutate(group = as.character(group))


# make prediction plot of probability of having CRISPR-Cas system vs genome length, faceted by species
length_plot <- ggplot(pred_lm1, aes(x = x, y = predicted), colour=x)+ 
  theme_classic() +
  geom_point(aes(colour=x), size=2) +
  geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high, colour=x), size=0.7) +  
  labs(x = "CRISPR-Cas", y = "Genome length") + 
  facet_grid(~ group) +
  theme(legend.position="none", plot.margin=unit(c(0.1,0.1,1,0.1),"cm")) +
  scale_colour_ghibli_d("PonyoMedium", direction = -1)

ggsave(plot=length_plot, "pred_plot_length_crispr.pdf", width=10, dpi=300, height=10)

####### prediction plots #########

# get max no. of AMR genes observed in all genomes
amr_range <- model_df %>%
  group_by(species) %>%
  mutate(max_amr = max(n_amr_genes)) %>%
  ungroup() %>%
  dplyr::select(species, max_amr) %>%
  distinct()

# filter for genomes that have acquired AMR genes (FT and MT only have intrinsic ones that the pipeline identifies, NM has max 1 AMR gene)
pred_glm_amr <- ggpredict(glm, terms = c("n_amr_genes", "species"), type='fixed', ci.lvl = 0.95) %>%
  filter(group != "FT" & group != "MT" & group != "NM")  %>%
  mutate(group = as.character(group))

# make prediction plot of probability of having CRISPR-Cas system vs AMR count, faceted by species
# 6 chosen as max 
AMR_pred_plot <- ggplot(pred_glm_amr, aes(x = x, y = predicted)) + 
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "ABR gene count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none") +
  scale_x_continuous(limits = c(0, 10), breaks=c(0, 5, 10))

# by species

# AB

pred_amr_AB <- pred_glm_amr %>%
  filter(group == "AB")

ABR_plot_AB <- ggplot(pred_amr_AB, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "ABR gene count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none", plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x = element_blank()) +
  scale_x_continuous(limits = c(0,36))

# EFm

# subset prediction dataframe
pred_amr_EFm <- pred_glm_amr %>%
  filter(group == "EFm")

ABR_plot_EFm <- ggplot(pred_amr_EFm, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "ABR gene count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,35))

# EFs

# subset prediction dataframe
pred_amr_EFs <- pred_glm_amr %>%
  filter(group == "EFs")

ABR_plot_EFs <- ggplot(pred_amr_EFs, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "ABR gene count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,26))

# PA

# subset prediction dataframe
pred_amr_PA <- pred_glm_amr %>%
  filter(group == "PA")

ABR_plot_PA <- ggplot(pred_amr_PA, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "ABR gene count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,32))

# SA

# subset prediction dataframe
pred_amr_SA <- pred_glm_amr %>%
  filter(group == "SA")

ABR_plot_SA <- ggplot(pred_amr_SA, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "ABR gene count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(), 
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,32))

# SE

# subset prediction dataframe
pred_amr_SE <- pred_glm_amr %>%
  filter(group == "SE")

ABR_plot_SE <- ggplot(pred_amr_SE, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "ABR gene count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,22))

# SP

# subset prediction dataframe
pred_amr_SP <- pred_glm_amr %>%
  filter(group == "SP")

ABR_plot_SP <- ggplot(pred_amr_SP, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "ABR gene count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,6))

# arrange all ABR plots
all_ABR <- ggarrange(ABR_plot_AB, ABR_plot_EFm, ABR_plot_EFs, ABR_plot_PA, ABR_plot_SA, ABR_plot_SE, ABR_plot_SP,
                     ncol = 7, nrow = 1)


# get max no. of IntI1 genes observed in all genomes
intI1_range <- model_df %>%
  group_by(species) %>%
  mutate(max_intI1 = max(n_intI1)) %>%
  ungroup() %>%
  dplyr::select(species, max_intI1) %>%
  distinct()

# intI1 genes - filter for genomes that have intI1
pred_glm_intI1 <- ggpredict(glm, terms = c("n_intI1", "species"), type='fixed', ci.lvl = 0.95) %>%
  filter(group == "AB" | group == "PA") %>%
  mutate(group = as.character(group))

intI1_pred_plot <- ggplot(pred_glm_intI1, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "intI1 copies", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0, 3), breaks=c(1, 2, 3))

# by species

# AB
pred_inti1_AB <- pred_glm_intI1 %>%
  filter(group == "AB")

inti1_plot_AB <- ggplot(pred_inti1_AB, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "IntI1 copy count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,3))

# PA

# subset prediction dataframe
pred_inti1_PA <- pred_glm_intI1 %>%
  filter(group == "PA")

inti1_plot_PA <- ggplot(pred_inti1_PA, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "IntI1 copy count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,6))


all_inti1 <- ggarrange(inti1_plot_AB, inti1_plot_PA,
                     ncol = 2, nrow = 1)


# ICEs 

# get max no. of ICEs observed in all genomes - 3
ICE_range <- model_df %>%
  group_by(species) %>%
  mutate(max_ICEs = max(n_ICEs)) %>%
  ungroup() %>%
  dplyr::select(species, max_ICEs) %>%
  distinct()

# prediction df - filter for genomes that have ICEs
pred_glm_ICEs <- ggpredict(glm, terms = c("n_ICEs", "species"), type='fixed', ci.lvl = 0.95) %>%
  mutate(group = as.character(group)) %>%
  filter(group == "PA" | group == "SP" | group == "EFm" | group == "EFs" | group == "SA")
  
ICES_pred_plot <- ggplot(pred_glm_ICEs, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +  
  labs(x = "ICE count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0, 3), breaks=c(1, 2, 3))

# by species

# EFm

# subset prediction dataframe
pred_ICE_EFm <- pred_glm_ICEs %>%
  filter(group == "EFm")

ICE_plot_EFm <- ggplot(pred_ICE_EFm, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "ICE count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank()) +
  scale_x_continuous(limits = c(0,3))

# EFs

# subset prediction dataframe
pred_ICE_EFs <- pred_glm_ICEs %>%
  filter(group == "EFs")

ICE_plot_EFs <- ggplot(pred_ICE_EFs, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "ICE count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,5))

# PA

# subset prediction dataframe
pred_ICE_PA <- pred_glm_ICEs %>%
  filter(group == "PA")

ICE_plot_PA <- ggplot(pred_ICE_PA, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "ICE count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,4))

# SA

# subset prediction dataframe
pred_ICE_SA <- pred_glm_ICEs %>%
  filter(group == "SA")

ICE_plot_SA <- ggplot(pred_ICE_SA, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "ICE count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(), 
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,11))

# SP

# subset prediction dataframe
pred_ICE_SP <- pred_glm_ICEs %>%
  filter(group == "SP")

ICE_plot_SP <- ggplot(pred_ICE_SP, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "ICE count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,3))

# arrange all ABR plots
all_ICE <- ggarrange(ICE_plot_EFm, ICE_plot_EFs, ICE_plot_PA, ICE_plot_SA, ICE_plot_SP,
                     ncol = 5, nrow = 1)


# plasmids 

# get max no. of ICEs observed in all genomes - 3
plasmid_range <- model_df %>%
  group_by(species) %>%
  mutate(max_plasmids = max(n_plasmid_replicons)) %>%
  ungroup() %>%
  dplyr::select(species, max_plasmids) %>%
  distinct()

# prediction df - filter for genomes that ever have > 2 plasmids
pred_glm_plasmids <- ggpredict(glm, terms = c("n_plasmid_replicons", "species"), type='fixed', ci.lvl = 0.95) %>%
  mutate(group = as.character(group)) %>%
  filter(group != "MT" & group != "NM" & group != "NG" & group != "MT" & group != "FT")

plasmids_pred_plot <- ggplot(pred_glm_plasmids, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +  
  labs(x = "Plasmid replicon count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0, 4), breaks=c(1, 2, 3,4))

# by species

# AB

pred_plasmid_AB <- pred_glm_plasmids %>%
  filter(group == "AB")

plasmid_plot_AB <- ggplot(pred_plasmid_AB, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "Plasmid replicon count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank()) +
  scale_x_continuous(limits = c(0,8))

# EFm

# subset prediction dataframe
pred_plasmid_EFm <- pred_glm_plasmids %>%
  filter(group == "EFm")

plasmid_plot_EFm <- ggplot(pred_plasmid_EFm, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "Plasmid replicon count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,18))

# EFs

# subset prediction dataframe
pred_plasmid_EFs <- pred_glm_plasmids %>%
  filter(group == "EFs")

plasmid_plot_EFs <- ggplot(pred_plasmid_EFs, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "Plasmid replicon count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,12), breaks=c(3,6,9))

# PA

# subset prediction dataframe
pred_plasmid_PA <- pred_glm_plasmids %>%
  filter(group == "PA")

plasmid_plot_PA <- ggplot(pred_plasmid_PA, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "Plasmid replicon count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,6))

# SA

# subset prediction dataframe
pred_plasmid_SA <- pred_glm_plasmids %>%
  filter(group == "SA")

plasmid_plot_SA <- ggplot(pred_plasmid_SA, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "Plasmid replicon count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(), 
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,8))

# SE

# subset prediction dataframe
pred_plasmid_SE <- pred_glm_plasmids %>%
  filter(group == "SE")

plasmid_plot_SE <- ggplot(pred_plasmid_SE, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "Plasmid replicon count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none",    plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,13))

# SP

# subset prediction dataframe
pred_plasmid_SP <- pred_glm_plasmids %>%
  filter(group == "SP")

plasmid_plot_SP <- ggplot(pred_plasmid_SP, aes(x = x, y = predicted)) + 
  geom_line() +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha=0.2) +  
  labs(x = "Plasmid replicon count", y = "CRISPR-Cas prob.") + 
  facet_grid(~ group) +
  theme_classic() + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5,0.75,1)) +
  theme(legend.position="none", plot.margin=unit(c(0.1,0.1,1,0.1),"cm"), axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits = c(0,4))

# arrange all plasmid plots
all_plasmids <- ggarrange(plasmid_plot_AB, plasmid_plot_EFm, plasmid_plot_EFs, plasmid_plot_PA, plasmid_plot_SA, plasmid_plot_SE, plasmid_plot_SP,
                     ncol = 7, nrow = 1)


# arrange all plot panels into one figure
middle_row <- plot_grid(all_ICE, all_inti1, labels = c('B', "C"), label_size = 15,
                        rel_widths = c(2.5,1))
bottom_row <- plot_grid(all_plasmids, labels = c('D'), label_size = 15)

x_labels <- c("No. of ABR genes", "No. of ICEs", "IntI1 copies", "No. of plasmid replicons")

plot_arrange <- plot_grid(all_ABR,
                          middle_row, 
                          bottom_row, 
                          labels = c('A', ''), 
                          rel_heights = c(2,1,1), label_size = 15, ncol = 1)

final_plots <- plot_arrange + 
  draw_text(x_labels, x = c(0.45,0.35,0.82,0.45), y = c(.53, .28,.28, 0.03), hjust = 0, size=12)
                                
ggsave(plot=final_plots, "CRISPR_binom_models.pdf", width = 27, dpi=300, height = 22, units = "cm")

