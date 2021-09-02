library(tidyverse)
library(ghibli)
library(ggpubr)

# creates tidy dataframe and descriptive plots  

setwd("/Volumes/Rosalind/genomic_analysis/phil_trans_21")

# read in datafiles 

# set column names 

crisprcas_cols <- c("ignore", "Contig","Operon","Operon_Pos","Prediction","CRISPRs","Distances","Prediction_Cas","Prediction_CRISPRs", "id")
arrays_cols <- c("ignore", "Contig","CRISPR","Start","End","Consensus_repeat", "N_repeats", "Repeat_len", "Spacer_len_avg", "Repeat_identity",
                 "Spacer_identity", "Spacer_len_sem", "Trusted", "Prediction", "Subtype", "Subtype_probability", "id")
amr_cols <- c("ignore", "FILE","SEQUENCE","START","END","STRAND","GENE","COVERAGE", "COVERAGE_MAP", "GAPS", "%COVERAGE",
              "%IDENTITY", "DATABASE", "ACCESSION", "PRODUCT", "RESISTANCE", "id")
acr_cols<- c('index', 'qseqid','seqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qseq','id')

# read in dataframes and generate same id column for all for merging

genome_lengths <- read_csv("results/final_dfs/genome_lengths.csv", col_names = c("id", "length")) %>%
  separate(col=id, into=c("resources", "genomes", "taxid", "id"), sep="/") %>%
  dplyr::select(-resources, -genomes) %>%
  transform(id=str_replace(id,".fna","")) %>%
  drop_na()
  
crisprcas_raw <- read_csv("results/final_dfs/crisprcas.csv", col_names = crisprcas_cols) %>%
  separate(col=id, into=c("results", "crispr", "taxid", "id"), sep="/") %>%
  dplyr::select(-ignore, -results, -crispr)

arrays_raw <- read_csv("results/final_dfs/crispr_arrays.csv", col_names = arrays_cols) %>%
  separate(col=id, into=c("results", "crispr", "taxid", "id"), sep="/") %>%
  dplyr::select(-ignore, -results, -crispr)

amr_raw <- read_csv("results/final_dfs/amr.csv", col_names=amr_cols) %>%
  separate(col=id, into=c("results", "amr", "taxid", "id"), sep="/") %>%
  dplyr::select(-ignore, -results, -amr) %>%
  transform(id=str_replace(id,".tab",""))

plasmids_raw <- read_csv("results/final_dfs/plasmids.csv", col_names=amr_cols) %>%
  separate(col=id, into=c("results", "plasmid", "taxid", "id"), sep="/") %>%
  dplyr::select(-ignore, -results, -plasmid) %>%
  transform(id=str_replace(id,".tab",""))

ICEs_raw <- read_csv("results/final_dfs/ICEs.csv", col_names=amr_cols) %>%
  separate(col=id, into=c("results", "plasmid", "taxid", "id"), sep="/") %>%
  dplyr::select(-ignore, -results, -plasmid) %>%
  transform(id=str_replace(id,".tab",""))

intI1_raw <- read_csv("results/final_dfs/intI1.csv", col_names=amr_cols) %>%
  separate(col=id, into=c("results", "intI1", "taxid", "id"), sep="/") %>%
  dplyr::select(-ignore, -results, -intI1) %>%
  transform(id=str_replace(id,".tab",""))

acrs_raw <- read_csv("rebuttal/all_acrs.csv", skip=1, col_names=acr_cols) %>%
  separate(col=id, into=c("results", "rebuttal", "acrs", 'id'), sep="/") %>%
  dplyr::select( -results, -rebuttal, -acrs, -index, -qseq) %>%
  transform(id=str_replace(id,"_acrs.tab","")) %>%
  mutate(taxid = "287")


# create dataframe of pathogen lifestyles and taxids

# type <- c("accidental/opportunistic", "obligate", "accidental/opportunistic", "accidental/opportunistic",
#           "obligate", "accidental/opportunistic", "obligate", "accidental/opportunistic","accidental/opportunistic",
#           "accidental/opportunistic", "accidental/opportunistic")
taxid <- c("1282", "1280", "263", "287", "470", "485", "487", "1773", "1314", "1351", "1352")
species <- c("SE", "SA", "FT", "PA", "AB", "NG", "NM", "MT", "SP", "EFs", "EFm")
  
pathogen <- data.frame(taxid, species)

# tidy data

# Table 1 - counts of genomes for each taxid
 genomes_totals <- genome_lengths %>%
  count(taxid, name="total")

# make df with all plasmids, virulence, genome lengths and crispr-cas presence/absence

amr_tidy <- amr_raw %>%
  dplyr::select(GENE, taxid, id) %>%
  count(taxid, id, name="n_amr_genes") 

plasmids_tidy <- plasmids_raw %>%
  dplyr::select(GENE, taxid, id) %>%
  count(taxid, id, name="n_plasmid_replicons") 

ICEs_tidy <- ICEs_raw %>%
  dplyr::select(taxid, id) %>%
  count(taxid, id, name="n_ICEs")

intI1_tidy <- intI1_raw %>%
  dplyr::select(GENE, taxid, id) %>%
  count(taxid, id, name="n_intI1") 

crisprcas <- crisprcas_raw %>%
  dplyr::select(Prediction, taxid, id) %>%
  count(Prediction,id,taxid, name="n_crisprcas_systems") 

arrays <- arrays_raw %>%
  dplyr::select(Prediction, taxid, id, N_repeats) %>%
  mutate(total_spacers = N_repeats - 1) %>% # get spacers per array
  group_by(id, taxid, Prediction) %>%
  mutate(total_spacers = sum(total_spacers)) %>%# get total spacers per genome for each CRISPR-Cas type
  ungroup() %>%
  select(-N_repeats) %>%
  distinct()

acrs <- acrs_raw %>%
  filter(pident > 90) %>%
  dplyr::select(taxid, id) %>%
  count(taxid, id, name="n_acrs")

acr_types <- acrs_raw %>%
  filter(pident > 90) %>%
  dplyr::select(seqid) %>%
  count(seqid, name="n_acrs")


# join together all dataframes
join_all <- full_join(crisprcas, genome_lengths) %>%
  left_join(., amr_tidy) %>%
  left_join(., plasmids_tidy) %>%
  left_join(., ICEs_tidy) %>%
  left_join(., intI1_tidy) %>%
  left_join(., pathogen) %>%
  left_join(., arrays) %>%
  left_join(., acrs)
  
# merge type IV systems into one category, add binomial presence/absence columns 
all_df <- join_all %>%
  mutate(Prediction = str_replace(Prediction, "IV-A1", "IV")) %>%
  mutate(Prediction = str_replace(Prediction, "IV-A2", "IV")) %>%
  mutate(Prediction = str_replace(Prediction, "IV-A3", "IV")) %>%
  mutate(Prediction = replace_na(Prediction, "None")) %>%
  replace(is.na(.), 0) %>%
  mutate(n_crisprcas_systems, crisprcas_binomial = ifelse(n_crisprcas_systems > 0, 1, 0)) %>%
  mutate(n_amr_genes, amr_binomial = ifelse(n_amr_genes > 0, 1, 0)) %>%
  mutate(n_plasmid_replicons, plasmids_binomial = ifelse(n_plasmid_replicons > 0, 1, 0)) %>%
  mutate(n_ICEs, ICEs_binomial = ifelse(n_ICEs > 0, 1, 0)) %>%
  mutate(n_intI1, intI1_binomial = ifelse(n_intI1 > 0, 1, 0)) %>%
  mutate(n_acrs, acrs_binomial = ifelse(n_acrs > 0, 1, 0)) %>%
  mutate(taxid=as.factor(taxid)) %>%
  mutate(species=as.character(species)) %>%
  mutate(length=as.numeric(length))

write_csv(all_df, 'rebuttal/all_df_rebuttal.csv')

# count dfs and plot of acrs

acr_count <-  all_df %>%
  select(-Prediction, -total_spacers) %>%
  distinct() %>%
  count(taxid, species, acrs_binomial, name="count") %>%
  mutate(acrs_binomial, acrs_binomial = ifelse(acrs_binomial == 1, "present", "absent")) %>% 
  left_join(., genomes_totals) %>%
  rowwise() %>%
  mutate(prop = count/total)

# bayesian phylogenetic models for PA with acrs vs. versus ICEs and each CRISPR-Cas system type

spacers_all_df <- read_csv("results/final_dfs/spacers_all.csv")

# filter spacer df: remove species with no CRISPR or no acquired AMR
spacers <- spacers_all_df %>%
  dplyr::select(-X1, -proportion, -target, -spacer_count_per_target) %>% # this will leave duplicate rows
  distinct() %>% # remove duplicates 
  filter(species == "PA") %>%
  mutate(taxid, taxid = (as.factor(taxid)))

# get counts of rows for each level of Prediction
count_spacers <- spacers %>%
  count(Prediction)

# rejoin with all_df to get data for species with no CRISPR-Cas system
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
outliers_PA_tree <- read_tsv("results/final_trees/287_dist.tab", col_names = c("Ref_ID", "Query_ID", "mash_dist", "p_value", "matching_hashes")) %>%
  filter(mash_dist > 0.1) %>%
  separate(col=Query_ID, into=c("taxid", "id"), sep="/") %>%
  separate(col=id, into=c("id", "ext"), sep=".fna") %>%
  dplyr::select(id)

PA_tree <- read.tree('results/final_trees/287_tree_1d.tre')

# subset model dataframe to just include PA
# remove any genomes for which more than one CRISPR-Cas system is present
model_df_type_PA <- model_df_type %>%
  filter(species == "PA") %>%
  filter(id %in% (outliers_PA_tree$id) == FALSE) %>% # remove outliers
  rename(animal = id) %>% # rename to work as MCMCglmm pedigree
  dplyr::select(Prediction, animal, n_amr_genes, total_spacers, n_ICEs, n_plasmid_replicons, n_acrs) %>%
  filter(animal %in% PA_tree$tip.label) %>% # make sure all genomes are in tree
  group_by(animal) %>%
  filter(n() == 1) %>% # only include genomes with one CRISPR-Cas system type
  ungroup() %>%
  filter(Prediction != "II-C" & Prediction != "IV") %>%
  mutate(Prediction= fct_relevel(Prediction, "None"))

model_df_spacers_PA <- model_df_type_PA %>% 
  filter(Prediction != "None") 

# histogram of response variable
ggplot(model_df_type_PA, aes(x=n_amr_genes)) +
  geom_histogram()

ggplot(model_df_spacers_PA, aes(x=n_amr_genes)) +
  geom_histogram()

## models 

# acrs
phylo_PA_type_acrs <- MCMCglmm(n_amr_genes ~ Prediction*n_acrs, family="poisson", data=as.data.frame(model_df_type_PA), pedigree=PA_tree, burnin = 20000,  nitt = 30000, thin = 5)

phylo_PA_type_acrs_summary <- summary(phylo_PA_type_acrs)
plot(phylo_PA_type_acrs)


# ICEs
phylo_PA_type_ICEs <- MCMCglmm(n_amr_genes ~ Prediction*n_ICEs, family="poisson", data=as.data.frame(model_df_type_PA), pedigree=PA_tree, burnin = 20000,  nitt = 30000, thin = 5)

phylo_PA_type_ICEs_summary <- summary(phylo_PA_type_ICEs)
plot(phylo_PA_type_ICEs)


# plots

# acrs

# prediction df for type
preddf_type_PA <- expand.grid(Prediction = c("None", "I-C", "I-E", "I-F"),
                              n_acrs = c(0, 1, 2),
                              n_amr_genes = 0)

pred_obj_type_PA <- predict(phylo_PA_type_acrs, newdata=preddf_type_PA, type="response", interval="confidence")

preddf_type_PA$n_amr_genes <- pred_obj_type_PA[,1]
preddf_type_PA$lower <- pred_obj_type_PA[,2]
preddf_type_PA$upper <- pred_obj_type_PA[,3]

plot_PA_type <- preddf_type_PA %>%
  filter(Prediction != "IV") 

PA_type_plot_acrs <- ggplot(plot_PA_type, aes(x=n_acrs, y=n_amr_genes, colour=Prediction, fill=Prediction)) +
  geom_line() +
  geom_ribbon(aes(x = n_acrs, ymin = lower, ymax = upper), width=0.2, alpha=0.5) +
  labs(x = "No. of type I-C acrs", y = "Predicted no. of ABR genes", colour="CRISPR-Cas type", fill="CRISPR-Cas type") +
  theme_classic() +
  scale_fill_manual(values = c("#D8AF39FF","#AE93BEFF", "#E7A79BFF", "#AD8152FF")) +
  scale_colour_manual(values = c("#D8AF39FF","#AE93BEFF", "#E7A79BFF", "#AD8152FF")) +
  scale_x_continuous(breaks=seq(0,2,1)) +
  theme(plot.margin=unit(c(0.7,0.2,0.2,0.2),"cm"))

# ICEs

# prediction df for type
preddf_type_PA_ICE <- expand.grid(Prediction = c("None", "I-C", "I-E", "I-F"),
                              n_amr_genes = 0,
                              n_ICEs = c(0,1, 2, 3, 4)) 

pred_obj_type_PA_ICE <- predict(phylo_PA_type_ICEs, newdata=preddf_type_PA_ICE, type="response", interval="confidence")

preddf_type_PA_ICE$n_amr_genes <- pred_obj_type_PA_ICE[,1]
preddf_type_PA_ICE$lower <- pred_obj_type_PA_ICE[,2]
preddf_type_PA_ICE$upper <- pred_obj_type_PA_ICE[,3]

PA_type_plot_ICEs <- ggplot(preddf_type_PA_ICE, aes(x=n_ICEs, y=n_amr_genes, colour=Prediction, fill=Prediction)) +
  geom_line() +
  geom_ribbon(aes(x = n_ICEs, ymin = lower, ymax = upper), alpha=0.5) +
  labs(x = "No. of ICEs", y = "Predicted no. of ABR genes", colour="CRISPR-Cas type", fill="CRISPR-Cas type") +
  theme_classic() +
  scale_fill_manual(values = c("#D8AF39FF","#AE93BEFF", "#E7A79BFF", "#AD8152FF")) +
  scale_colour_manual(values = c("#D8AF39FF","#AE93BEFF", "#E7A79BFF", "#AD8152FF")) +
  theme(plot.margin=unit(c(0.7,0.2,0.2,0.2),"cm"))

arranged <- ggarrange(PA_type_plot_acrs, PA_type_plot_ICEs,
                            ncol = 2, nrow = 1,
                            labels=c('A', 'B'),
                            common.legend = TRUE, legend = "bottom")

ggsave(plot=arranged, "PA_acrs_ICEs.pdf", dpi=300, width = 40, height = 20, units = "cm")

### get model estimates

# type
PA_type_est_acrs <- phylo_PA_type_acrs_summary$solutions
PA_type_est_ICEs <- phylo_PA_type_ICEs_summary$solutions

write.csv(PA_type_est_acrs, file = "results/phylo_type_rebuttal_acrs_PA.csv")

write.csv(PA_type_est_ICEs, file = "results/phylo_type_rebuttal_ICEs_PA.csv")

