library(tidyverse)
library(ghibli)
library(MuMIn)
library(bbmle)
library(DHARMa)
library(ggpubr)
library(pscl)
library(MCMCglmm)
library(ape)
library(phangorn)
library(picante)

######### CRISPR-Cas type/spacer models ###########

# read in organised spacer dataframe (created in 'spacer_analysis.R')
spacers_all_df <- read_csv("results/spacers_all.csv")

# filter spacer df: remove species with no CRISPR or no acquired AMR
spacers <- spacers_all_df %>%
  dplyr::select(-X1, -proportion, -target, -spacer_count_per_target) %>% # this will leave duplicate rows
  distinct() %>% # remove duplicates 
  filter(species != "NG" & species != "FT" & species != "MT")

# get counts of rows for each level of Prediction
count_spacers <- spacers %>%
  count(Prediction)

# rejoin with all_df to get data for species with no CRISPR-Cas system
all_df <- read_csv("results/all_df.csv")

spacers_type <- full_join(all_df, spacers)

# get counts of rows for each level of Prediction
count_type <- spacers_type %>%
  count(Prediction)

# get counts of rows for each level of ABR per species
count_amr <- spacers_type %>%
  count(n_amr_genes, species)

# specify names for model dataframes to be used throughout
model_df_type <- spacers_type
model_df_spacers <- spacers

#### split by species

######## PA #########

## type

# load original mash sketch and remove outliers
outliers_PA_tree <- read_tsv("final_trees/287_dist.tab", col_names = c("Ref_ID", "Query_ID", "mash_dist", "p_value", "matching_hashes")) %>%
  filter(mash_dist > 0.1) %>%
  separate(col=Query_ID, into=c("taxid", "id"), sep="/") %>%
  separate(col=id, into=c("id", "ext"), sep=".fna") %>%
  dplyr::select(id)

PA_tree <- read.tree('final_trees/287_tree_1d.tre')

# subset model dataframe to just include PA
# remove any genomes for which more than one CRISPR-Cas system is present
model_df_type_PA <- model_df_type %>%
  filter(species == "PA") %>%
  filter(id %in% (outliers_PA_tree$id) == FALSE) %>% # remove outliers
  rename(animal = id) %>% # rename to work as MCMCglmm pedigree
  dplyr::select(Prediction, animal, n_amr_genes, total_spacers, n_ICEs, n_plasmid_replicons)  %>%
  filter(animal %in% PA_tree$tip.label) %>% # make sure all genomes are in tree
  group_by(animal) %>%
  filter(n() == 1) %>% # only include genomes with one CRISPR-Cas system type
  ungroup() %>%
  filter(Prediction != "II-C") %>%
  mutate(Prediction= fct_relevel(Prediction, "None"))

model_df_spacers_PA <- model_df_type_PA %>% 
  filter(Prediction != "None") 

# histogram of response variable
ggplot(model_df_type_PA, aes(x=n_amr_genes)) +
  geom_histogram()

ggplot(model_df_spacers_PA, aes(x=n_amr_genes)) +
  geom_histogram()

## models 
phylo_PA_type <- MCMCglmm(n_amr_genes ~ Prediction, family="poisson", data=as.data.frame(model_df_type_PA), pedigree=PA_tree, burnin = 20000,  nitt = 30000, thin = 5)
phylo_PA_spacers <- MCMCglmm(n_amr_genes ~ total_spacers + Prediction + total_spacers*Prediction, family="poisson", data=as.data.frame(model_df_spacers_PA), pedigree=PA_tree, burnin = 20000,  nitt = 30000, thin = 5)

phylo_PA_type_summary <- summary(phylo_PA_type)
phylo_type_plot_PA <- plot(phylo_PA_type)

phylo_PA_spacers_summary <- summary(phylo_PA_spacers)
phylo_spacers_plot_PA <- (phylo_PA_spacers)

## type plus spacers

# plots

# prediction df for type
preddf_type_PA <- expand.grid(Prediction = c("None", "I-C", "I-E", "I-F", "IV"),
                              n_amr_genes = 0) 

pred_obj_type_PA <- predict(phylo_PA_type, newdata=preddf_type_PA, type="response", interval="confidence")

preddf_type_PA$n_amr_genes <- pred_obj_type_PA[,1]
preddf_type_PA$lower <- pred_obj_type_PA[,2]
preddf_type_PA$upper <- pred_obj_type_PA[,3]

# prediction df for spacers*type

preddf_spacers_PA <- expand.grid(Prediction = c("IV","I-C", "I-E", "I-F"),
                                 total_spacers = c(1:72),
                                 n_amr_genes = 0)

pred_obj_spacers_PA <- predict(phylo_PA_spacers, newdata=preddf_spacers_PA, type="response", interval="confidence")

preddf_spacers_PA$n_amr_genes <- pred_obj_spacers_PA[,1]
preddf_spacers_PA$lower <- pred_obj_spacers_PA[,2]
preddf_spacers_PA$upper <- pred_obj_spacers_PA[,3]

plot_PA_spacers <- preddf_spacers_PA %>%
  filter(Prediction != "IV")

PA_type_plot <- ggplot(preddf_type_PA, aes(x=Prediction, y=n_amr_genes)) +
  geom_point() +
  geom_errorbar(aes(x = Prediction, ymin = lower, ymax = upper), alpha=0.5, width=0.2) +
  labs(x = "CRISPR-Cas system type", y = "Predicted no. of ABR genes") +
  theme_classic() +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(0.7,0.2,0.2,0.2),"cm"))

PA_spacers_plot <- ggplot(plot_PA_spacers, aes(x=total_spacers, y=n_amr_genes, colour=Prediction, fill=Prediction)) +
  geom_line() +
  geom_ribbon(aes(x = total_spacers, ymin = lower, ymax = upper), alpha=0.5) +
  labs(x = "No. of spacers", y = "Predicted no. of ABR genes", fill="CRISPR-Cas system type", colour= "CRISPR-Cas system type") +
  theme_classic() +
  theme(legend.position="top", legend.title = element_blank()) +
  scale_fill_manual(values = c("#D8AF39FF","#AE93BEFF", "#E7A79BFF")) +
  scale_colour_manual(values = c("#D8AF39FF","#AE93BEFF", "#E7A79BFF")) +
  theme(plot.margin=unit(c(0.7,0.2,0.2,0.2),"cm"))

######## Efm #########

## type

# load original mash sketch and remove outliers
outliers_EFm_tree <- read_tsv("final_trees/1352_dist.tab", col_names = c("Ref_ID", "Query_ID", "mash_dist", "p_value", "matching_hashes")) %>%
  filter(mash_dist > 0.1) %>%
  separate(col=Query_ID, into=c("1", "taxid", "id"), sep="/") %>% # remove 1 later
  separate(col=id, into=c("id", "ext"), sep=".fna") %>%
  dplyr::select(id)

EFm_tree <- read.tree('final_trees/1352_tree_1d.tre')

# subset model dataframe to just include PA
# remove any genomes for which more than one CRISPR-Cas system is present
model_df_type_EFm <- model_df_type %>%
  filter(species == "EFm") %>%
  filter(id %in% (outliers_EFm_tree$id) == FALSE) %>% # remove outliers
  rename(animal = id) %>% # rename to work as MCMCglmm pedigree
  dplyr::select(Prediction, animal, n_amr_genes, total_spacers)  %>%
  filter(animal %in% EFm_tree$tip.label) %>% # make sure all genomes are in tree
  group_by(animal) %>%
  filter(n() == 1) %>% # only include genomes with one CRISPR-Cas system type
  ungroup() %>%
  mutate(Prediction= fct_relevel(Prediction, "None"))

model_df_spacers_EFm <- model_df_type_EFm %>% 
  filter(Prediction != "None")

# histogram of response variable
ggplot(model_df_type_EFm, aes(x=n_amr_genes)) +
  geom_histogram()

ggplot(model_df_spacers_EFm, aes(x=n_amr_genes)) +
  geom_histogram()

## models 
phylo_EFm_type <- MCMCglmm(n_amr_genes ~ Prediction, family="poisson", data=as.data.frame(model_df_type_EFm), pedigree=EFm_tree, burnin = 20000,  nitt = 30000, thin = 5)
phylo_EFm_spacers <- MCMCglmm(n_amr_genes ~ total_spacers, family="poisson", data=as.data.frame(model_df_spacers_EFm), pedigree=EFm_tree, burnin = 20000,  nitt = 30000, thin = 5)

phylo_EFm_type_summary <- summary(phylo_EFm_type)
phylo_type_plot_EFm <- plot(phylo_EFm_type)

phylo_EFm_spacers_summary <- summary(phylo_EFm_spacers)
phylo_spacers_plot_EFm <- plot(phylo_EFm_spacers)

## type plus spacers

# plots

# prediction df for type
preddf_type_EFm <- expand.grid(Prediction = c("None","II-A"),
                              n_amr_genes = 0) 

pred_obj_type_EFm <- predict(phylo_EFm_type, newdata=preddf_type_EFm, type="response", interval="confidence")

preddf_type_EFm$n_amr_genes <- pred_obj_type_EFm[,1]
preddf_type_EFm$lower <- pred_obj_type_EFm[,2]
preddf_type_EFm$upper <- pred_obj_type_EFm[,3]

# prediction df for spacers*type

preddf_spacers_EFm <- expand.grid(Prediction = c("II-A"),
                                 total_spacers = c(1:20),
                                 n_amr_genes = 0)

pred_obj_spacers_EFm <- predict(phylo_EFm_spacers, newdata=preddf_spacers_EFm, type="response", interval="confidence")

preddf_spacers_EFm$n_amr_genes <- pred_obj_spacers_EFm[,1]
preddf_spacers_EFm$lower <- pred_obj_spacers_EFm[,2]
preddf_spacers_EFm$upper <- pred_obj_spacers_EFm[,3]

EFm_type_plot <- ggplot(preddf_type_EFm, aes(x=Prediction, y=n_amr_genes)) +
  geom_point() +
  geom_errorbar(aes(x = Prediction, ymin = lower, ymax = upper), alpha=0.5, width=0.2) +
  labs(x = "CRISPR-Cas system type", y = "Predicted no. of ABR genes") +
  theme_classic() +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(0.7,0.2,0.2,0.2),"cm"))

EFm_spacers_plot <- ggplot(preddf_spacers_EFm, aes(x=total_spacers, y=n_amr_genes, colour=Prediction, fill=Prediction)) +
  geom_line() +
  geom_ribbon(aes(x = total_spacers, ymin = lower, ymax = upper), alpha=0.5) +
  labs(x = "No. of spacers", y = "Predicted no. of ABR genes", fill="CRISPR-Cas type", colour= "CRISPR-Cas type") +
  theme_classic() +
  theme(legend.position="top", legend.title = element_blank()) +
  scale_fill_manual(values = c("#E75B64FF")) +
  scale_colour_manual(values = c("#E75B64FF")) +
  theme(plot.margin=unit(c(0.7,0.2,0.2,0.2),"cm"))

###### EFs ########

## type

# load original mash sketch and remove outliers
outliers_EFs_tree <- read_tsv("final_trees/1351_dist.tab", col_names = c("Ref_ID", "Query_ID", "mash_dist", "p_value", "matching_hashes")) %>%
  filter(mash_dist > 0.1) %>%
  separate(col=Query_ID, into=c("1", "taxid", "id"), sep="/") %>% # remove 1 later
  separate(col=id, into=c("id", "ext"), sep=".fna") %>%
  dplyr::select(id)

EFs_tree <- read.tree('final_trees/1351_tree_1d.tre')

# subset model dataframe to just include PA
# remove any genomes for which more than one CRISPR-Cas system is present
model_df_type_EFs <- model_df_type %>%
  filter(species == "EFs") %>%
  filter(id %in% (outliers_EFs_tree$id) == FALSE) %>% # remove outliers
  rename(animal = id) %>% # rename to work as MCMCglmm pedigree
  dplyr::select(Prediction, animal, n_amr_genes, total_spacers)  %>%
  filter(animal %in% EFs_tree$tip.label) %>% # make sure all genomes are in tree
  group_by(animal) %>%
  filter(n() == 1) %>% # only include genomes with one CRISPR-Cas system type
  ungroup()  %>%
  mutate(Prediction= fct_relevel(Prediction, "None"))

model_df_spacers_EFs <- model_df_type_EFs %>% 
  filter(Prediction != "None")

# histogram of response variable
ggplot(model_df_type_EFs, aes(x=n_amr_genes)) +
  geom_histogram()

ggplot(model_df_spacers_EFs, aes(x=n_amr_genes)) +
  geom_histogram()

## models 
phylo_EFs_type <- MCMCglmm(n_amr_genes ~ Prediction, family="poisson", data=as.data.frame(model_df_type_EFs), pedigree=EFs_tree, burnin = 20000,  nitt = 30000, thin = 5)
phylo_EFs_spacers <- MCMCglmm(n_amr_genes ~ total_spacers, family="poisson", data=as.data.frame(model_df_spacers_EFs), pedigree=EFs_tree, burnin = 20000,  nitt = 30000, thin = 5)

phylo_EFs_type_summary <- summary(phylo_EFs_type)
plot(phylo_EFs_type)

phylo_EFs_spacers_summary <- summary(phylo_EFs_spacers)
plot(phylo_EFs_spacers)

## type plus spacers

# plots

# prediction df for type
preddf_type_EFs <- expand.grid(Prediction = c("None","II-A"),
                              n_amr_genes = 0) 

pred_obj_type_EFs <- predict(phylo_EFs_type, newdata=preddf_type_EFs, type="response", interval="confidence")

preddf_type_EFs$n_amr_genes <- pred_obj_type_EFs[,1]
preddf_type_EFs$lower <- pred_obj_type_EFs[,2]
preddf_type_EFs$upper <- pred_obj_type_EFs[,3]

# prediction df for spacers*type

preddf_spacers_EFs <- expand.grid(Prediction = c("II-A"),
                                 total_spacers = c(1:36),
                                 n_amr_genes = 0)

pred_obj_spacers_EFs <- predict(phylo_EFs_spacers, newdata=preddf_spacers_EFs, type="response", interval="confidence")

preddf_spacers_EFs$n_amr_genes <- pred_obj_spacers_EFs[,1]
preddf_spacers_EFs$lower <- pred_obj_spacers_EFs[,2]
preddf_spacers_EFs$upper <- pred_obj_spacers_EFs[,3]

EFs_type_plot <- ggplot(preddf_type_EFs, aes(x=Prediction, y=n_amr_genes)) +
  geom_point() +
  geom_errorbar(aes(x = Prediction, ymin = lower, ymax = upper), alpha=0.5, width=0.2) +
  labs(x = "CRISPR-Cas system type", y = "Predicted no. of ABR genes") +
  theme_classic() +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(0.7,0.2,0.2,0.2),"cm"))

EFs_spacers_plot <- ggplot(preddf_spacers_EFs, aes(x=total_spacers, y=n_amr_genes, colour=Prediction, fill=Prediction)) +
  geom_line() +
  geom_ribbon(aes(x = total_spacers, ymin = lower, ymax = upper), alpha=0.5) +
  labs(x = "No. of spacers", y = "Predicted no. of ABR genes", fill="CRISPR-Cas type", colour= "CRISPR-Cas type") +
  theme_classic() +
  theme(legend.position="top", legend.title = element_blank()) +
  scale_fill_manual(values = c("#E75B64FF")) +
  scale_colour_manual(values = c("#E75B64FF")) +
  theme(plot.margin=unit(c(0.7,0.2,0.2,0.2),"cm"))

###### AB ########

## type

# load original mash sketch and remove outliers
outliers_AB_tree <- read_tsv("final_trees/470_dist.tab", col_names = c("Ref_ID", "Query_ID", "mash_dist", "p_value", "matching_hashes")) %>%
  filter(mash_dist > 0.1) %>%
  separate(col=Query_ID, into=c("taxid", "id"), sep="/") %>% # remove 1 later
  separate(col=id, into=c("id", "ext"), sep=".fna") %>%
  dplyr::select(id)

AB_tree <- read.tree('final_trees/470_tree_10.tre')

# subset model dataframe to just include PA
# remove any genomes for which more than one CRISPR-Cas system is present
model_df_type_AB <- model_df_type %>%
  filter(species == "AB") %>%
  filter(id %in% (outliers_AB_tree$id) == FALSE) %>% # remove outliers
  rename(animal = id) %>% # rename to work as MCMCglmm pedigree
  dplyr::select(Prediction, animal, n_amr_genes, total_spacers)  %>%
  filter(animal %in% AB_tree$tip.label) %>% # make sure all genomes are in tree
  group_by(animal) %>%
  filter(n() == 1) %>% # only include genomes with one CRISPR-Cas system type
  ungroup()  %>%
  mutate(Prediction= fct_relevel(Prediction, "None"))

model_df_spacers_AB <- model_df_type_AB %>% 
  filter(Prediction != "None")

# histogram of response variable
ggplot(model_df_type_AB, aes(x=n_amr_genes)) +
  geom_histogram()

ggplot(model_df_spacers_AB, aes(x=n_amr_genes)) +
  geom_histogram()

## models 
phylo_AB_type <- MCMCglmm(n_amr_genes ~ Prediction, family="poisson", data=as.data.frame(model_df_type_AB), pedigree=AB_tree, burnin = 20000,  nitt = 30000, thin = 5)
phylo_AB_spacers <- MCMCglmm(n_amr_genes ~ total_spacers, family="poisson", data=as.data.frame(model_df_spacers_AB), pedigree=AB_tree, burnin = 20000,  nitt = 30000, thin = 5)

phylo_AB_type_summary <- summary(phylo_AB_type)
plot(phylo_AB_type)

phylo_AB_spacers_summary <- summary(phylo_AB_spacers)
plot(phylo_AB_spacers)

## type plus spacers

# plots

# prediction df for type
preddf_type_AB <- expand.grid(Prediction = c("None", "I-F"),
                               n_amr_genes = 0) 

pred_obj_type_AB <- predict(phylo_AB_type, newdata=preddf_type_AB, type="response", interval="confidence")

preddf_type_AB$n_amr_genes <- pred_obj_type_AB[,1]
preddf_type_AB$lower <- pred_obj_type_AB[,2]
preddf_type_AB$upper <- pred_obj_type_AB[,3]

# prediction df for spacers*type

preddf_spacers_AB <- expand.grid(Prediction = c("I-F"),
                                  total_spacers = c(1:223),
                                  n_amr_genes = 0)

pred_obj_spacers_AB <- predict(phylo_AB_spacers, newdata=preddf_spacers_AB, type="response", interval="confidence")

preddf_spacers_AB$n_amr_genes <- pred_obj_spacers_AB[,1]
preddf_spacers_AB$lower <- pred_obj_spacers_AB[,2]
preddf_spacers_AB$upper <- pred_obj_spacers_AB[,3]

AB_type_plot <- ggplot(preddf_type_AB, aes(x=Prediction, y=n_amr_genes)) +
  geom_point() +
  geom_errorbar(aes(x = Prediction, ymin = lower, ymax = upper), alpha=0.5, width=0.2) +
  labs(x = "CRISPR-Cas system type", y = "Predicted no. of ABR genes") +
  theme_classic() +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(0.7,0.2,0.2,0.2),"cm"))

AB_spacers_plot <- ggplot(preddf_spacers_AB, aes(x=total_spacers, y=n_amr_genes, colour=Prediction, fill=Prediction)) +
  geom_line() +
  geom_ribbon(aes(x = total_spacers, ymin = lower, ymax = upper), alpha=0.5) +
  labs(x = "No. of spacers", y = "Predicted no. of ABR genes", fill="CRISPR-Cas type", colour= "CRISPR-Cas type") +
  theme_classic() +
  theme(legend.position="top", legend.title = element_blank()) +
  scale_fill_manual(values = c("#E7A79BFF")) +
  scale_colour_manual(values = c("#E7A79BFF")) +
  theme(plot.margin=unit(c(0.7,0.2,0.2,0.2),"cm"))

###### SP ########

## type

# load original mash sketch and remove outliers
outliers_SP_tree <- read_tsv("final_trees/1314_dist.tab", col_names = c("Ref_ID", "Query_ID", "mash_dist", "p_value", "matching_hashes")) %>%
  filter(mash_dist > 0.1) %>%
  separate(col=Query_ID, into=c("taxid", "id"), sep="/") %>% # remove 1 later
  separate(col=id, into=c("id", "ext"), sep=".fna") %>%
  dplyr::select(id)

SP_tree <- read.tree('final_trees/1314_tree_1.tre')

# subset model dataframe to just include PA
# remove any genomes for which more than one CRISPR-Cas system is present
model_df_type_SP <- model_df_type %>%
  filter(species == "SP") %>%
  filter(id %in% (outliers_SP_tree$id) == FALSE) %>% # remove outliers
  rename(animal = id) %>% # rename to work as MCMCglmm pedigree
  dplyr::select(Prediction, animal, n_amr_genes, total_spacers)  %>%
  filter(animal %in% SP_tree$tip.label) %>% # make sure all genomes are in tree
  group_by(animal) %>%
  filter(n() == 1) %>% # only include genomes with one CRISPR-Cas system type
  ungroup()  %>%
  mutate(Prediction= fct_relevel(Prediction, "None"))

model_df_spacers_SP <- model_df_type_SP %>% 
  filter(Prediction != "None")

# histogram of response variable
ggplot(model_df_type_SP, aes(x=n_amr_genes)) +
  geom_histogram()

ggplot(model_df_spacers_SP, aes(x=n_amr_genes)) +
  geom_histogram()

## models 
phylo_SP_type <- MCMCglmm(n_amr_genes ~ Prediction, family="poisson", data=as.data.frame(model_df_type_SP), pedigree=SP_tree, burnin = 20000,  nitt = 30000, thin = 5)
phylo_SP_spacers <- MCMCglmm(n_amr_genes ~ total_spacers + Prediction + total_spacers*Prediction, family="poisson", data=as.data.frame(model_df_spacers_SP), pedigree=SP_tree, burnin = 20000,  nitt = 30000, thin = 5)

phylo_SP_type_summary <- summary(phylo_SP_type)
plot(phylo_SP_type)

phylo_SP_spacers_summary <- summary(phylo_SP_spacers)
plot(phylo_SP_spacers)

## type plus spacers

# plots

# prediction df for type
preddf_type_SP <- expand.grid(Prediction = c("None", "I-C", "II-A"),
                              n_amr_genes = 0) 

pred_obj_type_SP <- predict(phylo_SP_type, newdata=preddf_type_SP, type="response", interval="confidence")

preddf_type_SP$n_amr_genes <- pred_obj_type_SP[,1]
preddf_type_SP$lower <- pred_obj_type_SP[,2]
preddf_type_SP$upper <- pred_obj_type_SP[,3]

# prediction df for spacers*type

preddf_spacers_SP <- expand.grid(Prediction = c("II-A", "I-C"),
                                 total_spacers = c(1:17),
                                 n_amr_genes = 0)

pred_obj_spacers_SP <- predict(phylo_SP_spacers, newdata=preddf_spacers_SP, type="response", interval="confidence")

preddf_spacers_SP$n_amr_genes <- pred_obj_spacers_SP[,1]
preddf_spacers_SP$lower <- pred_obj_spacers_SP[,2]
preddf_spacers_SP$upper <- pred_obj_spacers_SP[,3]

SP_type_plot <- ggplot(preddf_type_SP, aes(x=Prediction, y=n_amr_genes)) +
  geom_point() +
  geom_errorbar(aes(x = Prediction, ymin = lower, ymax = upper), alpha=0.5, width=0.2) +
  labs(x = "CRISPR-Cas system type", y = "Predicted no. of ABR genes") +
  theme_classic() +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(0.7,0.2,0.2,0.2),"cm"))

SP_spacers_plot <- ggplot(preddf_spacers_SP, aes(x=total_spacers, y=n_amr_genes, colour=Prediction, fill=Prediction)) +
  geom_line() +
  geom_ribbon(aes(x = total_spacers, ymin = lower, ymax = upper), alpha=0.5) +
  labs(x = "No. of spacers", y = "Predicted no. of ABR genes", fill="CRISPR-Cas type", colour= "CRISPR-Cas type") +
  theme_classic() +
  theme(legend.position="top", legend.title = element_blank()) +
  scale_fill_manual(values = c("#E75B64FF", "#D8AF39FF")) +
  scale_colour_manual(values = c("#E75B64FF", "#D8AF39FF")) +
  theme(plot.margin=unit(c(0.7,0.2,0.2,0.2),"cm"))


###### SE ####### 

## type

# load original mash sketch and remove outliers
outliers_SE_tree <- read_tsv("final_trees/1282_dist.tab", col_names = c("Ref_ID", "Query_ID", "mash_dist", "p_value", "matching_hashes")) %>%
  filter(mash_dist > 0.1) %>%
  separate(col=Query_ID, into=c("taxid", "id"), sep="/") %>% # remove 1 later
  separate(col=id, into=c("id", "ext"), sep=".fna") %>%
  dplyr::select(id)

SE_tree <- read.tree('final_trees/1282_tree_1.tre')

# subset model dataframe to just include PA
# remove any genomes for which more than one CRISPR-Cas system is present
model_df_type_SE <- model_df_type %>%
  filter(species == "SE") %>%
  filter(id %in% (outliers_SE_tree$id) == FALSE) %>% # remove outliers
  rename(animal = id) %>% # rename to work as MCMCglmm pedigree
  dplyr::select(Prediction, animal, n_amr_genes, total_spacers)  %>%
  filter(animal %in% SE_tree$tip.label) %>% # make sure all genomes are in tree
  group_by(animal) %>%
  filter(n() == 1) %>% # only include genomes with one CRISPR-Cas system type
  ungroup() %>%
  filter(Prediction != "II-A") %>% #only 18 samples
  mutate(Prediction= fct_relevel(Prediction, "None"))

model_df_spacers_SE <- model_df_type_SE %>% 
  filter(Prediction != "None")

# histogram of response variable
ggplot(model_df_type_SE, aes(x=n_amr_genes)) +
  geom_histogram()

ggplot(model_df_spacers_SE, aes(x=n_amr_genes)) +
  geom_histogram()

# convert "zero" length edge 
SE_tree$edge.length[which(SE_tree$edge.length < 1e-16)] <- 1e-16

prior1 <- list(R = list(V = diag(1), nu = 0.002), B = list(mu = diag(1,2), V = diag(1e+08/2, 2)))
prior2 <- list(R = list(V = 0, nu = -2))

## models 
phylo_SE_type <- MCMCglmm(n_amr_genes ~ Prediction, family="poisson", data=as.data.frame(model_df_type_SE), pedigree=SE_tree, burnin = 20000,  nitt = 30000, thin = 5)
#phylo_SE_spacers <- MCMCglmm(n_amr_genes ~ total_spacers, family="poisson", data=as.data.frame(model_df_spacers_SE), prior=prior2, pedigree=SE_tree, burnin = 20000,  nitt = 30000, thin = 5)

phylo_SE_type_summary <- summary(phylo_SE_type)
plot(phylo_SE_type)

# phylo_SE_spacers_summary <- summary(phylo_SE_spacers)
# plot(phylo_SE_spacers)

## type plus spacers

# plots

# prediction df for type
preddf_type_SE <- expand.grid(Prediction = c("None","III-A"),
                               n_amr_genes = 0) 

pred_obj_type_SE <- predict(phylo_SE_type, newdata=preddf_type_SE, type="response", interval="confidence")

preddf_type_SE$n_amr_genes <- pred_obj_type_SE[,1]
preddf_type_SE$lower <- pred_obj_type_SE[,2]
preddf_type_SE$upper <- pred_obj_type_SE[,3]

SE_type_plot <- ggplot(preddf_type_SE, aes(x=Prediction, y=n_amr_genes)) +
  geom_point() +
  geom_errorbar(aes(x = Prediction, ymin = lower, ymax = upper), alpha=0.5, width=0.2) +
  labs(x = "CRISPR-Cas system type", y = "Predicted no. of ABR genes") +
  theme_classic() +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(0.7,0.2,0.2,0.2),"cm"))

all_plots <- ggarrange(AB_type_plot, EFm_type_plot, EFs_type_plot,  PA_type_plot, SP_type_plot, SE_type_plot, 
                       AB_spacers_plot, EFm_spacers_plot, EFs_spacers_plot, PA_spacers_plot, SP_spacers_plot,
                       nrow=2, ncol=6, 
                       hjust=c(-1.8, -1,-1.2, -1.7, -1.9, -1.7),
                       labels=c('AB', 'EFm', 'EFs', 'PA', 'SP', 'SE')) +
  annotate("text", x = 0.01, y = 0.98, label = "A", size=8) + 
  annotate("text", x = 0.01, y = 0.44, label = "B", size=8)
  
ggsave(plot=all_plots, "plots/phylo_plots.pdf", width = 40, height = 17, dpi=300, units = "cm")


### get model estimates

# type
PA_type_est <- phylo_PA_type_summary$solutions
PA_type_write <- cbind(PA_type_est, c("PA")) 

EFm_type_est <- phylo_EFm_type_summary$solutions
EFm_type_write <- cbind(EFm_type_est, c("EFm")) 

EFs_type_est <- phylo_EFs_type_summary$solutions
EFs_type_write <- cbind(EFs_type_est, c("EFs")) 

AB_type_est <- phylo_AB_type_summary$solutions
AB_type_write <- cbind(AB_type_est, c("AB")) 

SP_type_est <- phylo_SP_type_summary$solutions
SP_type_write <- cbind(SP_type_est, c("SP")) 

SE_type_est <- phylo_SE_type_summary$solutions
SE_type_write <- cbind(SE_type_est, c("SE")) 

write_df_type <- rbind(PA_type_write, EFm_type_write, EFs_type_write, AB_type_write,
            SP_type_write, SE_type_write)

write.csv(write_df_type, file = "results/phylo_type.csv")

# spacers
PA_spacers_est <- phylo_PA_spacers_summary$solutions
PA_spacers_write <- cbind(PA_spacers_est, c("PA")) 

EFm_spacers_est <- phylo_EFm_spacers_summary$solutions
EFm_spacers_write <- cbind(EFm_spacers_est, c("EFm")) 

EFs_spacers_est <- phylo_EFs_spacers_summary$solutions
EFs_spacers_write <- cbind(EFs_spacers_est, c("EFs")) 

AB_spacers_est <- phylo_AB_spacers_summary$solutions
AB_spacers_write <- cbind(AB_spacers_est, c("AB")) 

SP_spacers_est <- phylo_SP_spacers_summary$solutions
SP_spacers_write <- cbind(SP_spacers_est, c("SP")) 

write_df_spacers <- rbind(PA_spacers_write, EFm_spacers_write, EFs_spacers_write, AB_spacers_write,
                       SP_spacers_write)

write.csv(write_df_spacers, file = "results/phylo_spacers.csv")



