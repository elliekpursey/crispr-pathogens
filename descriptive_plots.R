library(tidyverse)
library(ghibli)
library(ggpubr)

# creates tidy dataframe and descriptive plots  

setwd("results")

# read in datafiles 

# set column names 

crisprcas_cols <- c("ignore", "Contig","Operon","Operon_Pos","Prediction","CRISPRs","Distances","Prediction_Cas","Prediction_CRISPRs", "id")
arrays_cols <- c("ignore", "Contig","CRISPR","Start","End","Consensus_repeat", "N_repeats", "Repeat_len", "Spacer_len_avg", "Repeat_identity",
                 "Spacer_identity", "Spacer_len_sem", "Trusted", "Prediction", "Subtype", "Subtype_probability", "id")
amr_cols <- c("ignore", "FILE","SEQUENCE","START","END","STRAND","GENE","COVERAGE", "COVERAGE_MAP", "GAPS", "%COVERAGE",
              "%IDENTITY", "DATABASE", "ACCESSION", "PRODUCT", "RESISTANCE", "id")

# read in dataframes and generate same id column for all for merging

genome_lengths <- read_csv("genome_lengths.csv", col_names = c("id", "length")) %>%
  separate(col=id, into=c("resources", "genomes", "taxid", "id"), sep="/") %>%
  dplyr::select(-resources, -genomes) %>%
  transform(id=str_replace(id,".fna","")) %>%
  drop_na()
  
crisprcas_raw <- read_csv("crispr/crisprcas.csv", col_names = crisprcas_cols) %>%
  separate(col=id, into=c("results", "crispr", "taxid", "id"), sep="/") %>%
  dplyr::select(-ignore, -results, -crispr)

arrays_raw <- read_csv("crispr/crispr_arrays.csv", col_names = arrays_cols) %>%
  separate(col=id, into=c("results", "crispr", "taxid", "id"), sep="/") %>%
  dplyr::select(-ignore, -results, -crispr)

amr_raw <- read_csv("amr/amr.csv", col_names=amr_cols) %>%
  separate(col=id, into=c("results", "amr", "taxid", "id"), sep="/") %>%
  dplyr::select(-ignore, -results, -amr) %>%
  transform(id=str_replace(id,".tab",""))

plasmids_raw <- read_csv("plasmids/plasmids.csv", col_names=amr_cols) %>%
  separate(col=id, into=c("results", "plasmid", "taxid", "id"), sep="/") %>%
  dplyr::select(-ignore, -results, -plasmid) %>%
  transform(id=str_replace(id,".tab",""))

ICEs_raw <- read_csv("ICEs/ICEs.csv", col_names=amr_cols) %>%
  separate(col=id, into=c("results", "plasmid", "taxid", "id"), sep="/") %>%
  dplyr::select(-ignore, -results, -plasmid) %>%
  transform(id=str_replace(id,".tab",""))

intI1_raw <- read_csv("intI1/intI1.csv", col_names=amr_cols) %>%
  separate(col=id, into=c("results", "intI1", "taxid", "id"), sep="/") %>%
  dplyr::select(-ignore, -results, -intI1) %>%
  transform(id=str_replace(id,".tab",""))


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

# join together all dataframes
join_all <- full_join(crisprcas, genome_lengths) %>%
  left_join(., amr_tidy) %>%
  left_join(., plasmids_tidy) %>%
  left_join(., ICEs_tidy) %>%
  left_join(., intI1_tidy) %>%
  left_join(., pathogen) %>%
  left_join(., arrays)
  
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
  mutate(taxid=as.factor(taxid)) %>%
  mutate(species=as.character(species)) %>%
  mutate(length=as.numeric(length))

write_csv(all_df, 'all_df.csv')

# make df of spacers counts, crispr-cas system types and amr genes

# count dfs

crispr_count <- all_df %>%
  select(-Prediction, -total_spacers) %>%
  distinct() %>%
  count(taxid, species, crisprcas_binomial, name="count") %>%
  mutate(crisprcas_binomial, crisprcas_binomial = ifelse(crisprcas_binomial == 1, "present", "absent")) %>%
  left_join(., genomes_totals) %>%
  rowwise() %>%
  mutate(prop = count/total)

amr_count <-  all_df %>%
  select(-Prediction, -total_spacers) %>%
  distinct() %>%
  count(taxid, species, amr_binomial, name="count") %>%
  mutate(amr_binomial, amr_binomial = ifelse(amr_binomial == 1, "present", "absent")) %>% 
  left_join(., genomes_totals) %>%
  rowwise() %>%
  mutate(prop = count/total)

plasmid_count <-  all_df %>%
  select(-Prediction, -total_spacers) %>%
  distinct() %>%
  count(taxid, species, plasmids_binomial, name="count") %>%
  mutate(plasmids_binomial, plasmids_binomial = ifelse(plasmids_binomial == 1, "present", "absent")) %>% 
  left_join(., genomes_totals) %>%
  rowwise() %>%
  mutate(prop = count/total)

ICE_count <-  all_df %>%
  select(-Prediction, -total_spacers) %>%
  distinct() %>%
  count(taxid, species, ICEs_binomial, name="count") %>%
  mutate(ICEs_binomial, ICEs_binomial = ifelse(ICEs_binomial == 1, "present", "absent")) %>% 
  left_join(., genomes_totals) %>%
  rowwise() %>%
  mutate(prop = count/total)

intI1_count <- all_df %>%
  select(-Prediction, -total_spacers) %>%
  distinct() %>%
  count(taxid, species, intI1_binomial, name="count") %>%
  mutate(intI1_binomial, intI1_binomial = ifelse(intI1_binomial == 1, "present", "absent")) %>% 
  left_join(., genomes_totals) %>%
  rowwise() %>%
  mutate(prop = count/total)

prediction_counts <- all_df %>%
  dplyr::select(Prediction, species) %>%
  filter(Prediction != "None" & Prediction != "Unknown") %>%
  count(species, Prediction, name="count") %>%
  replace(is.na(.), "None") %>%
  left_join(., pathogen) %>%
  left_join(., genomes_totals) %>%
  rowwise() %>%
  mutate(prop = count/total) 

# descriptive plots

crispr_d <- ggplot(crispr_count, aes(x = species, y = prop, fill=crisprcas_binomial)) + 
  geom_col() +
  labs(x = "Species", y = "Proportion of genomes") + 
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_fill_ghibli_d("PonyoMedium", direction = -1) +
  theme(plot.margin=unit(c(0.7,0.7,0.7,0.7),"cm"))

amr_d <- ggplot(amr_count, aes(x = species, y = prop, fill = amr_binomial)) + 
  geom_col() +
  labs(x = "Species", y = "Proportion of genomes") + 
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_fill_ghibli_d("PonyoMedium", direction = -1) +
  theme(plot.margin=unit(c(0.7,0.7,0.7,0.7),"cm"))
  
plasmid_d <- ggplot(plasmid_count, aes(x = species, y = prop, fill = plasmids_binomial)) + 
  geom_col() +
  labs(x = "Species", y = "Proportion of genomes") + 
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_fill_ghibli_d("PonyoMedium", direction = -1) +
  theme(plot.margin=unit(c(0.7,0.7,0.7,0.7),"cm"))
 
intI1_d <- ggplot(intI1_count, aes(x = species, y = prop, fill = intI1_binomial)) + 
  geom_col() +
  labs(x = "Species", y = "Proportion of genomes") + 
  theme(legend.title = element_blank()) +
  theme_classic() +
  scale_fill_ghibli_d("PonyoMedium", direction = -1) +
  theme(plot.margin=unit(c(0.7,0.7,0.7,0.7),"cm"))

ICE_d <- ggplot(ICE_count, aes(x = species, y = prop, fill = ICEs_binomial)) + 
  geom_col() +
  labs(x = "Species", y = "Proportion of genomes") + 
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_fill_ghibli_d("PonyoMedium", direction = -1) +
  theme(plot.margin=unit(c(0.7,0.7,0.7,0.7),"cm"))

descript_plots <- ggarrange(crispr_d, amr_d, plasmid_d, intI1_d, ICE_d,
                            labels=c("CRISPR-Cas", "ABR genes", "Plasmid replicons",
                                      "intI1 copies", "ICEs"),
                            label.y=1,
                            align='h',
                            ncol = 3, nrow = 2,
                            common.legend = TRUE, legend = "bottom")

ggsave(plot=descript_plots, "desc_plots.pdf", dpi= 300, width = 30, height = 20, units = "cm")


# number of genomes with each type of CRISPR-Cas system
CRISPR_types <- ggplot(prediction_counts, aes(x=species, y = prop, fill=species)) + 
  theme_classic() +
  geom_col(position = "dodge") +
  labs(x = "CRISPR-Cas system type", y = "Proportion of genomes") + 
  facet_wrap(~ Prediction) +
  scale_fill_ghibli_d("MarnieMedium1", direction = 1) 

# plot crispr types per species

AB_type <- prediction_counts %>%
  filter(species == "AB")

AB_type_plot <- ggplot(AB_type, aes(x=Prediction, y = prop, fill=Prediction)) + 
  theme_classic() +
  geom_col() +
  labs(x = "CRISPR-Cas type", y = "Proportion of genomes") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_fill_manual(values = c("#E7A79BFF","#4C413FFF"))

PA_type <- prediction_counts %>%
  filter(species == "PA")

PA_type_plot <- ggplot(PA_type, aes(x=Prediction, y = prop, fill=Prediction)) + 
  theme_classic() +
  geom_col() +
  labs(x = "CRISPR-Cas type", y = "Proportion of genomes") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_fill_manual(values = c("#D8AF39FF","#AE93BEFF", "#E7A79BFF", "#4C413FFF", "#4C413FFF"))

EFm_type <- prediction_counts %>%
  filter(species == "EFm")

EFm_type_plot <- ggplot(EFm_type, aes(x=Prediction, y = prop, fill=Prediction)) + 
  theme_classic() +
  geom_col() +
  labs(x = "CRISPR-Cas type", y = "Proportion of genomes") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_fill_manual(values = c("#E75B64FF"))

EFs_type <- prediction_counts %>%
  filter(species == "EFs")

EFs_type_plot <- ggplot(EFs_type, aes(x=Prediction, y = prop, fill=Prediction)) + 
  theme_classic() +
  geom_col() +
  labs(x = "CRISPR-Cas type", y = "Proportion of genomes") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_fill_manual(values = c("#E75B64FF"))

SA_type <- prediction_counts %>%
  filter(species == "SA")

SA_type_plot <- ggplot(SA_type, aes(x=Prediction, y = prop, fill=Prediction)) + 
  theme_classic() +
  geom_col() +
  labs(x = "CRISPR-Cas type", y = "Proportion of genomes") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_fill_manual(values = c("#E75B64FF", "#278B9AFF"))

SE_type <- prediction_counts %>%
  filter(species == "SE")

SE_type_plot <- ggplot(SE_type, aes(x=Prediction, y = prop, fill=Prediction)) + 
  theme_classic() +
  geom_col() +
  labs(x = "CRISPR-Cas type", y = "Proportion of genomes") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_fill_manual(values = c("#E75B64FF", "#278B9AFF"))

SP_type <- prediction_counts %>%
  filter(species == "SP")

SP_type_plot <- ggplot(SP_type, aes(x=Prediction, y = prop, fill=Prediction)) + 
  theme_classic() +
  geom_col() +
  labs(x = "CRISPR-Cas type", y = "Proportion of genomes") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_fill_manual(values = c("#D8AF39FF", "#E75B64FF"))

FT_type <- prediction_counts %>%
  filter(species == "FT")

FT_type_plot <- ggplot(FT_type, aes(x=Prediction, y = prop, fill=Prediction)) + 
  theme_classic() +
  geom_col() +
  labs(x = "CRISPR-Cas type", y = "Proportion of genomes") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_fill_manual(values = c("#833437FF", "#C5A387FF"))

MT_type <- prediction_counts %>%
  filter(species == "MT")

MT_type_plot <- ggplot(MT_type, aes(x=Prediction, y = prop, fill=Prediction)) + 
  theme_classic() +
  geom_col() +
  labs(x = "CRISPR-Cas type", y = "Proportion of genomes") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_fill_manual(values = c("#278B9AFF"))

NM_type <- prediction_counts %>%
  filter(species == "NM")

NM_type_plot <- ggplot(NM_type, aes(x=Prediction, y = prop, fill=Prediction)) + 
  theme_classic() +
  geom_col() +
  labs(x = "CRISPR-Cas type", y = "Proportion of genomes") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_fill_manual(values = c("#D8AF39FF", "#819A7AFF"))

descript_plots_type <- ggarrange(AB_type_plot, PA_type_plot, SE_type_plot, SA_type_plot,
                                   NM_type_plot, EFs_type_plot, EFm_type_plot, SP_type_plot,
                                   MT_type_plot, FT_type_plot,
                                   ncol = 4, nrow = 3, legend  = "none" )

ggsave(plot=descript_plots_type, "desc_plots_crispr_type_distrib.pdf", dpi = 300, width = 40, height = 30, units = "cm")


# number of CRISPR repeats per CRISPR type and species
CRISPR_repeats_plot <- all_df %>%
  filter(Prediction != "None")

# plot crispr repeats per species

AB_repeats <- CRISPR_repeats_plot %>%
  filter(species == "AB")

AB_repeats_plot <- ggplot(AB_repeats, aes(x=Prediction, y = total_spacers, fill=Prediction)) + 
  theme_classic() +
  geom_boxplot() +
  labs(x = "CRISPR-Cas type", y = "Spacers per genome") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_fill_manual(values = c("#E7A79BFF","#4C413FFF"))

PA_repeats <- CRISPR_repeats_plot %>%
  filter(species == "PA")

PA_repeats_plot <- ggplot(PA_repeats, aes(x=Prediction, y = total_spacers, fill=Prediction)) + 
  theme_classic() +
  geom_boxplot() +
  labs(x = "CRISPR-Cas type", y = "Spacers per genome") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_fill_manual(values = c("#D8AF39FF","#AE93BEFF", "#E7A79BFF", "#4C413FFF", "#4C413FFF"))


EFm_repeats <- CRISPR_repeats_plot %>%
  filter(species == "EFm")

EFm_repeats_plot <- ggplot(EFm_repeats, aes(x=Prediction, y = total_spacers, fill=Prediction)) + 
  theme_classic() +
  geom_boxplot() +
  labs(x = "CRISPR-Cas type", y = "Spacers per genome") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_fill_manual(values = c("#E75B64FF"))

EFs_repeats <- CRISPR_repeats_plot %>%
  filter(species == "EFs")

EFs_repeats_plot <- ggplot(EFs_repeats, aes(x=Prediction, y = total_spacers, fill=Prediction)) + 
  theme_classic() +
  geom_boxplot() +
  labs(x = "CRISPR-Cas type", y = "Spacers per genome") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_fill_manual(values = c("#E75B64FF"))

SA_repeats <- CRISPR_repeats_plot %>%
  filter(species == "SA")

SA_repeats_plot <- ggplot(SA_repeats, aes(x=Prediction, y = total_spacers, fill=Prediction)) + 
  theme_classic() +
  geom_boxplot() +
  labs(x = "CRISPR-Cas type", y = "Spacers per genome") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_fill_manual(values = c("#E75B64FF", "#278B9AFF"))

SE_repeats <- CRISPR_repeats_plot %>%
  filter(species == "SE")

SE_repeats_plot <- ggplot(SE_repeats, aes(x=Prediction, y = total_spacers, fill=Prediction)) + 
  theme_classic() +
  geom_boxplot() +
  labs(x = "CRISPR-Cas type", y = "Spacers per genome") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_fill_manual(values = c("#E75B64FF", "#278B9AFF"))

SP_repeats <- CRISPR_repeats_plot %>%
  filter(species == "SP")

SP_repeats_plot <- ggplot(SP_repeats, aes(x=Prediction, y = total_spacers, fill=Prediction)) + 
  theme_classic() +
  geom_boxplot() +
  labs(x = "CRISPR-Cas type", y = "Spacers per genome") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_fill_manual(values = c("#D8AF39FF", "#E75B64FF"))

FT_repeats <- CRISPR_repeats_plot %>%
  filter(species == "FT")

FT_repeats_plot <- ggplot(FT_repeats, aes(x=Prediction, y = total_spacers, fill=Prediction)) + 
  theme_classic() +
  geom_boxplot() +
  labs(x = "CRISPR-Cas type", y = "Spacers per genome") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_fill_manual(values = c("#833437FF", "#C5A387FF"))


MT_repeats <- CRISPR_repeats_plot %>%
  filter(species == "MT")

MT_repeats_plot <- ggplot(MT_repeats, aes(x=Prediction, y = total_spacers, fill=Prediction)) + 
  theme_classic() +
  geom_boxplot() +
  labs(x = "CRISPR-Cas type", y = "Spacers per genome") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_fill_manual(values = c("#278B9AFF"))

NM_repeats <- CRISPR_repeats_plot %>%
  filter(species == "NM")

NM_repeats_plot <- ggplot(NM_repeats, aes(x=Prediction, y = total_spacers, fill=Prediction)) + 
  theme_classic() +
  geom_boxplot() +
  labs(x = "CRISPR-Cas type", y = "Spacers per genome") + 
  facet_wrap(~ species) +
  theme(legend.position="none") +
  scale_fill_manual(values = c("#D8AF39FF", "#819A7AFF"))


descript_plots_crispr <- ggarrange(AB_repeats_plot, PA_repeats_plot, SE_repeats_plot, SA_repeats_plot,
                                   NM_repeats_plot, EFs_repeats_plot, EFm_repeats_plot, SP_repeats_plot,
                                   MT_repeats_plot, FT_repeats_plot,
                            ncol = 4, nrow = 3, legend  = "none" )

ggsave(plot=descript_plots_crispr, "desc_plots_spacers_distrib.pdf", dpi = 300, width = 40, height = 30, units = "cm")
