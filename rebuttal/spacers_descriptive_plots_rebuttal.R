library(tidyverse)
library(ggeffects)
library(lme4)
library(ghibli)
library(MuMIn)
library(bbmle)
library(DHARMa)
library(ggpubr)
library(glmmTMB)

# for testing only
setwd("/Volumes/Rosalind/genomic_analysis/phil_trans_21/results")

# read in datafiles 

# set column names 

arrays_cols <- c("ignore", "Contig","CRISPR","Start","End","Consensus_repeat", "N_repeats", "Repeat_len", "Spacer_len_avg", "Repeat_identity",
                 "Spacer_identity", "Spacer_len_sem", "Trusted", "Prediction", "Subtype", "Subtype_probability", "id")
amr_cols <- c("ignore", "FILE","SEQUENCE","START","END","STRAND","GENE","COVERAGE", "COVERAGE_MAP", "GAPS", "%COVERAGE",
              "%IDENTITY", "DATABASE", "ACCESSION", "PRODUCT", "RESISTANCE", "id")

# load dataframes

# load crispr arrays, select only needed columns, amalgamate type IV predictions as in main full df
crispr <- read_csv("final_dfs/crispr_arrays.csv", col_names = arrays_cols)

arrays <- crispr %>%
  separate(col=id, into=c("results", "crispr", "taxid", "id"), sep="/") %>%
  select(id, Contig, CRISPR, Prediction, taxid, N_repeats) %>%
  mutate(Prediction = str_replace(Prediction, "IV-A1", "IV")) %>%
  mutate(Prediction = str_replace(Prediction, "IV-A2", "IV")) %>%
  mutate(Prediction = str_replace(Prediction, "IV-A3", "IV")) %>%
  mutate(total_spacers=N_repeats-1) %>% # get total_spacers from repeats - 1
  select(-N_repeats)

# read in spacer df
spacers <- read_csv("final_dfs/spacer_targets.csv") %>%
  distinct(qseqid, .keep_all=TRUE) %>%    # remove duplicate spacer hits so only one hit per spacer
  separate(col=qseqid, into=c("CRISPR"), sep=":") # split column to get CRISPR array id

# filter for only spacers with 100% match, count spacers for each target type 
spacers_join <- spacers %>%
  filter(pident == 100) %>%
  count(target, CRISPR, name="n_targeting_spacers") %>%
  right_join(., arrays, by="CRISPR") %>% # join with arrays 
  filter(Prediction != "Unknown") %>% # remove unknown predictions for CRISPR-Cas
  mutate(target = replace_na(target, 'unknown')) %>% # replace NA for target column to ensure all arrays remain in df
  mutate(n_targeting_spacers = replace_na(n_targeting_spacers, 0)) # ditto for n_targeting spacers

# calculate number of each spacer type for each array
spacer_types_arrays <- spacers_join %>%
  pivot_wider(names_from = target, values_from = n_targeting_spacers) %>%
  replace(is.na(.), 0)
  
# sum per genome and Prediction
spacer_types_genomes <- spacer_types_arrays %>%
  group_by(Prediction, id) %>%
  summarise(across(c(ICE, plasmid, virus, unknown, total_spacers), sum)) %>%
  ungroup()

# join with main dataframe

all_df <- read_csv("final_dfs/all_df.csv") 

join_all_df <- spacer_types_genomes %>% 
  left_join(., all_df, by=c('id', 'Prediction', 'total_spacers')) %>%
  drop_na() # small number of cases (n=4) where NAs due to a Type III CRISPR having Type II spacers

# calculate number of each spacer type for each CRISPR-Cas type for each genome
spacers_all <- join_all_df %>%
  mutate(unknown = total_spacers - ICE - plasmid - virus) %>%
  pivot_longer(cols = c(virus, ICE, plasmid, unknown), names_to='target', values_to="spacer_count_per_target") %>% # convert to long format for spacer targets
  mutate(proportion = spacer_count_per_target/total_spacers) # add proportions for each genome

# counts of unknown spacer targets in dataset

unknown_counts <- spacers_all %>%
  filter(target=="unknown") %>%
  mutate(sum = sum(spacer_count_per_target)) 

ICE_counts <- spacers_all %>%
  filter(target=="ICE") %>%
  mutate(sum = sum(spacer_count_per_target)) %>%
  mutate(db_length_kb = 28377422/1000) %>%
  mutate(spacers_per_kb = sum/db_length_kb)

plasmid_counts <- spacers_all %>%
  filter(target=="plasmid") %>%
  mutate(sum = sum(spacer_count_per_target)) %>%
  mutate(db_length_kb = 856388404/1000) %>%
  mutate(spacers_per_kb = sum/db_length_kb)

virus_counts <- spacers_all %>%
  filter(target=="virus")  %>%
  mutate(sum = sum(spacer_count_per_target)) %>%
  mutate(db_length_kb = 352271305/1000) %>%
  mutate(spacers_per_kb = sum/db_length_kb)

total_no_spacers <- spacers_all %>%
  filter(target=="unknown") %>% # get one row per genome to sum total spacers
  mutate(total_spacers_in_dataset = sum(total_spacers))

# write to file
write.csv(spacers_all, "final_dfs/spacers_all_95id.csv")

# descriptive plots
plot_df <- spacers_all
  
plot_df$target <- factor(plot_df$target, levels = c("virus", "plasmid", "ICE", "unknown"))

# by species
species <- ggplot(plot_df, aes(x=species, y=proportion, fill=target)) +
  geom_boxplot(position=position_dodge()) +
  theme_classic() +
  labs(x="CRISPR-Cas system type", y="Proportion of spacers", fill="Spacer target") +
  scale_fill_ghibli_d("PonyoMedium", direction = -1) +
  theme(plot.margin=unit(c(0.7,0.7,0.7,0.7),"cm"))

# by system type
type <- ggplot(plot_df, aes(x=Prediction, y=proportion, fill=target)) +
  geom_boxplot(position=position_dodge()) +
  theme_classic() +
  labs(x="CRISPR-Cas system type", y="Proportion of spacers", fill="Spacer target") +
  scale_fill_ghibli_d("PonyoMedium", direction = -1) +
  theme(plot.margin=unit(c(0.7,0.7,0.7,0.7),"cm"))

descript_plots_overview <- ggarrange(species, type,
                            ncol = 2, nrow = 1,
                            labels=c('A', 'B'),
                            common.legend = TRUE, legend = "bottom")

#ggsave(plot=descript_plots_overview, "desc_plots_spacers_type_95pid.pdf", width = 30, height = 15, units = "cm")
  

# plot by species + system type - filter df to create separate dataframe for each species
# (Cannot facet_grid as the CRISPR-Cas system types do not match up and its a mess!)

PA_df <- plot_df %>%
  filter(species == "PA")

PA_d <- ggplot(PA_df, aes(x=Prediction, y=proportion, colour=target)) +
  geom_boxplot(position=position_dodge()) +
  theme_classic() +
  facet_grid(~species) +
  labs(x="CRISPR-Cas system type", y="Proportion of spacers", colour="Spacer target") +
  scale_colour_ghibli_d("PonyoMedium", direction = 1) 

AB_df <- plot_df %>%
  filter(species == "AB")

AB_d <- ggplot(AB_df, aes(x=Prediction, y=proportion, colour=target)) +
  geom_boxplot(position=position_dodge(width=1.1)) +
  theme_classic() +
  facet_grid(~species) +
  labs(x="CRISPR-Cas system type", y="Proportion of spacers", colour="Spacer target") +
  scale_colour_ghibli_d("PonyoMedium", direction = 1) 

SE_df <- plot_df %>%
  filter(species == "SE")

SE_d <- ggplot(SE_df, aes(x=Prediction, y=proportion, colour=target)) +
  geom_boxplot(position=position_dodge(width=1.1)) +
  theme_classic() +
  facet_grid(~species) +
  labs(x="CRISPR-Cas system type", y="Proportion of spacers", colour="Spacer target") +
  scale_colour_ghibli_d("PonyoMedium", direction = 1) 

SP_df <- plot_df %>%
  filter(species == "SP")

SP_d <- ggplot(SP_df, aes(x=Prediction, y=proportion, colour=target)) +
  geom_boxplot(position=position_dodge(width=1.1)) +
  theme_classic() +
  facet_grid(~species) +
  labs(x="CRISPR-Cas system type", y="Proportion of spacers", colour="Spacer target") +
  scale_colour_ghibli_d("PonyoMedium", direction = 1) 

NM_df <- plot_df %>%
  filter(species == "NM")

NM_d <- ggplot(NM_df, aes(x=Prediction, y=proportion, colour=target)) +
  geom_boxplot(position=position_dodge(width=1.1)) +
  theme_classic() +
  facet_grid(~species) +
  labs(x="CRISPR-Cas system type", y="Proportion of spacers", colour="Spacer target") +
  scale_colour_ghibli_d("PonyoMedium", direction = 1) 

Efm_df <- plot_df %>%
  filter(species == "EFm")

Efm_d <- ggplot(Efm_df, aes(x=Prediction, y=proportion, colour=target)) +
  geom_boxplot(position=position_dodge(width=1.1)) +
  theme_classic() +
  facet_grid(~species) +
  labs(x="CRISPR-Cas system type", y="Proportion of spacers", colour="Spacer target") +
  scale_colour_ghibli_d("PonyoMedium", direction = 1) 

Efs_df <- plot_df %>%
  filter(species == "EFs")

Efs_d <- ggplot(Efs_df, aes(x=Prediction, y=proportion, colour=target)) +
  geom_boxplot(position=position_dodge(width=1.1)) +
  theme_classic() +
  facet_grid(~species) +
  labs(x="CRISPR-Cas system type", y="Proportion of spacers", colour="Spacer target") +
  scale_colour_ghibli_d("PonyoMedium", direction = 1) 

SA_df <- plot_df %>%
  filter(species == "SA")

SA_d <- ggplot(SA_df, aes(x=Prediction, y=proportion, colour=target)) +
  geom_boxplot(position=position_dodge(width=1.1)) +
  theme_classic() +
  facet_grid(~species) +
  labs(x="CRISPR-Cas system type", y="Proportion of spacers", colour="Spacer target") +
  scale_colour_ghibli_d("PonyoMedium", direction = 1) 


descript_plots <- ggarrange(AB_d, PA_d, SE_d, SA_d, NM_d, Efs_d, Efm_d, SP_d,
                            ncol = 4, nrow = 2,
                            common.legend = TRUE, legend = "bottom")

ggsave(plot=descript_plots, "desc_plots_spacers_95_pid.pdf", dpi=300, width = 40, height = 20, units = "cm")

