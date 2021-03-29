library(tidyverse)
library(ggeffects)
library(lme4)
library(ghibli)
library(MuMIn)
library(bbmle)
library(DHARMa)
library(ggpubr)

# for testing only
setwd("results")

arrays_cols <- c("ignore", "Contig","CRISPR","Start","End","Consensus_repeat", "N_repeats", "Repeat_len", "Spacer_len_avg", "Repeat_identity",
                 "Spacer_identity", "Spacer_len_sem", "Trusted", "Prediction", "Subtype", "Subtype_probability", "id")
amr_cols <- c("ignore", "FILE","SEQUENCE","START","END","STRAND","GENE","COVERAGE", "COVERAGE_MAP", "GAPS", "%COVERAGE",
              "%IDENTITY", "DATABASE", "ACCESSION", "PRODUCT", "RESISTANCE", "id")

# load dataframes

# load crispr arrays, select only needed columns, amalgamate type IV predictions as in main full df
crispr <- read_csv("crispr/crispr_arrays.csv", col_names = arrays_cols)

# tidy crispr array data and add spacer counts per array
arrays <- crispr %>%
  separate(col=id, into=c("results", "crispr", "taxid", "id"), sep="/") %>% # split filename column into usable data ie taxid and id
  dplyr::select(id, Contig, CRISPR, Prediction, taxid, N_repeats) %>%
  mutate(Prediction = str_replace(Prediction, "IV-A1", "IV")) %>% # convert all type IV subtypes to type IV (little data)
  mutate(Prediction = str_replace(Prediction, "IV-A2", "IV")) %>%
  mutate(Prediction = str_replace(Prediction, "IV-A3", "IV")) %>%
  mutate(total_spacers=N_repeats-1) %>% # get total_spacers from repeats - 1
  dplyr::select(-N_repeats)

# read in MGEs with AMR data frame 
MGEs_with_AMR <- read_csv("MGEs_with_AMR.csv") %>%
  filter(pident > 90) %>% 
  dplyr::select(qseqid, sseqid, target) %>%
  distinct(sseqid, target) %>% # remove multiple hits per MGE
  dplyr::select(-target) %>%
  mutate(MGE_w_AMR_targeting_spacers=1) # add column for summing in next step

# read in spacer target dataframe and tidy
spacer_targets <- read_csv("spacer_targets.csv") 

spacers_tidy <- spacer_targets %>%
  separate(col=qseqid, into=c("CRISPR", "spacer_id"), sep=":") # split column to get CRISPR array id

# join spacer df with MGEs with AMR dataframe, by matching sseqid (ie the MGE name (spacer sequence target or AMR hit))
spacers_against_MGEs_w_AMR_convert_binom <- spacers_tidy %>% 
  left_join(., MGEs_with_AMR) %>%
  filter(pident == 100) %>% # get spacers only with 100% matches
  dplyr::select(CRISPR, spacer_id, target, MGE_w_AMR_targeting_spacers) %>%
  replace(is.na(.), 0) %>%
  group_by(CRISPR, spacer_id, target) %>% # account for there being multiple hits per spacer by summing hits per spacer  
  summarise(across(everything(), sum)) %>% # sum all AMR hits together for each spacer
  ungroup() %>%
  mutate(MGE_w_AMR_targeting_spacers, MGE_w_AMR_targeting_spacers = ifelse(MGE_w_AMR_targeting_spacers ==0, 0, 1)) # convert to binomial i.e. if multiple spacers match a target with amr, that counts as blocking 1 amr gene
  
spacers_against_MGEs_w_AMR <- spacers_against_MGEs_w_AMR_convert_binom %>%
  dplyr::select(-spacer_id, -target) %>%
  group_by(CRISPR) %>% # group again to get counts for each spacer target type per CRISPR array
  summarise(across(everything(), sum)) %>%
  ungroup() 

# get total spacers targeting amr across whole dataset for results
total_spacers_targeting_AMR <- spacers_against_MGEs_w_AMR %>%
  mutate(sum_spacers = sum(MGE_w_AMR_targeting_spacers))

# filter for only spacers with 100% match, count total MGE-targeting spacers for each array
spacers_against_all_MGEs <- spacers_tidy %>%
  filter(pident == 100) %>%
  dplyr::select(CRISPR, spacer_id, target) %>%
  distinct() %>% # remove duplicate spacer hits so only one hit per spacer
  dplyr::select(-spacer_id) %>%
  count(CRISPR, name="MGE_targeting_spacers")

total_spacers_targeting_MGEs <- spacers_against_all_MGEs %>%
  mutate(sum_spacers = sum(MGE_targeting_spacers))

# join dataframes for both spacer types (target any known MGE and MGE with AMR)
spacers_all <- left_join(spacers_against_all_MGEs, spacers_against_MGEs_w_AMR)

# join with all CRISPR array data   
spacers_join <- spacers_all %>%
  right_join(., arrays, by="CRISPR") %>% # join with arrays 
  filter(Prediction != "Unknown") %>% # remove unknown predictions for CRISPR-Cas
  replace(is.na(.), 0)

# sum per genome and Prediction
spacer_types_genomes <- spacers_join %>%
  dplyr::select(id, MGE_targeting_spacers, Prediction, MGE_w_AMR_targeting_spacers, total_spacers) %>%
  group_by(id, Prediction) %>%
  summarise(across(everything(), sum)) %>%
  mutate(MGE_wo_AMR_targeting_spacers = MGE_targeting_spacers - MGE_w_AMR_targeting_spacers) %>%
  ungroup()

# join with main dataframe

all_df <- read_csv("all_df.csv") %>%
  dplyr::select(id, n_amr_genes, species)

final_df <- spacer_types_genomes %>% 
  left_join(., all_df, by=c('id')) %>%
  mutate(n_amr_genes, n_amr_genes = as.numeric(n_amr_genes)) 

# total number of spacers targeting vectors of AMR

p <- ggplot(final_df, aes(x=MGE_wo_AMR_targeting_spacers, y=n_amr_genes, colour='MGE targeting')) + 
  geom_point()
            
p +  geom_point(aes(x=MGE_w_AMR_targeting_spacers, y=n_amr_genes, colour='MGE with AMR-targeting'))

### models proportion

model_df <- final_df %>%
  filter(species != "NG" & species != "MT" & species != "FT" & species != "NM") 

ggplot(model_df, aes(x=n_amr_genes)) +
  geom_density()

glm_mge_amr <- MCMCglmm(n_amr_genes ~ MGE_w_AMR_targeting_spacers + species, family="poisson", data=as.data.frame(model_df), burnin = 10000,  nitt = 20000, thin = 5)
glm_mge_no_amr <- MCMCglmm(n_amr_genes ~ MGE_wo_AMR_targeting_spacers + species, family="poisson", data=as.data.frame(model_df), burnin = 10000,  nitt = 20000, thin = 5)
summary_w_amr <- summary(glm_mge_amr)
summary_wo_amr <- summary(glm_mge_no_amr)

plot(glm_mge_amr)
plot(glm_mge_no_amr)


# Prediction dataframe for MGEs with AMR
preddf_amr<- expand.grid(species = c("AB", "EFm", "EFs", "PA", "SA", "SE", "SP"),
                              MGE_w_AMR_targeting_spacers = c(0:22),
                              n_amr_genes = 0) 

pred_obj_amr <- predict(glm_mge_amr, newdata=preddf_amr, type="response", interval="confidence")

preddf_amr$n_amr_genes <- pred_obj_amr[,1]
preddf_amr$lower <- pred_obj_amr[,2]
preddf_amr$upper <- pred_obj_amr[,3]

preddf_amr$facet = 'MGE with ABR'

# Prediction dataframe for MGEs without AMR
preddf_no_amr<- expand.grid(species = c("AB", "EFm", "EFs", "PA", "SA", "SE", "SP"),
                         MGE_wo_AMR_targeting_spacers = c(0:22),
                         n_amr_genes = 0) 

pred_obj_no_amr <- predict(glm_mge_no_amr, newdata=preddf_no_amr, type="response", interval="confidence")

preddf_no_amr$n_amr_genes <- pred_obj_no_amr[,1]
preddf_no_amr$lower <- pred_obj_no_amr[,2]
preddf_no_amr$upper <- pred_obj_no_amr[,3]

preddf_no_amr$facet = 'MGE without ABR'

# plot

plot1 <- ggplot(preddf_amr, aes(x=MGE_w_AMR_targeting_spacers, y=n_amr_genes, colour=facet, fill=facet)) +
  geom_line() +
  geom_ribbon(aes(x = MGE_w_AMR_targeting_spacers, ymin = lower, ymax = upper), alpha=0.5) +
  labs(x = "No. of spacers with known targets", y = "Predicted no. of ABR genes", colour="Spacer target", fill="Spacer target") +
  theme_classic() +
  facet_grid(~species) +
  scale_colour_manual(values=c("#F4ADB3FF", "#AE93BEFF")) +
  scale_fill_manual(values=c("#F4ADB3FF", "#AE93BEFF"))

combined_plot <- plot1 +
  geom_line(data=preddf_no_amr, aes(x = MGE_wo_AMR_targeting_spacers, y=n_amr_genes, colour=facet)) +
  geom_ribbon(data=preddf_no_amr, aes(x = MGE_wo_AMR_targeting_spacers, ymin = lower, ymax = upper, fill=facet), alpha=0.5)

ggsave(plot=combined_plot, "spacer_type_models.pdf", width = 25, dpi = 300, height = 20, units = "cm")

# write model output data frames

write.csv(summary_w_amr$solutions, file = "targets_w_amr.csv")
write.csv(summary_wo_amr$solutions, file = "targets_wo_amr.csv")