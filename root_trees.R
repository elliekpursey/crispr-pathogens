library(ape)
library(phangorn)
library(tidyverse)

#### SP
SP_tree <- read.tree('1314.dnd')

# load original mash sketch and remove outliers
outliers_SP_tree <- read_tsv("1314_dist.tab", col_names = c("Ref_ID", "Query_ID", "mash_dist", "p_value", "matching_hashes")) %>%
  filter(mash_dist > 0.1) %>%
  separate(col=Query_ID, into=c("taxid", "id"), sep="/") %>%
  separate(col=id, into=c("id", "ext"), sep=".fna") %>%
  dplyr::select(id)

SP_tree_pruned <- drop.tip(SP_tree, outliers_SP_tree$id)

SP_tree_rooted <- root(SP_tree_pruned, outgroup="GCF_900638415.1_57635_F01_genomic", resolve.root = TRUE)

SP_calibration <- makeChronosCalib(SP_tree_rooted, node="root")

print("Strep pyogenes trees")

print("correlated")
SP_tree_1 <- chronos(SP_tree_rooted, lambda = 1, model = "correlated", calibration = SP_calibration, control = chronos.control())
SP_tree_10 <- chronos(SP_tree_rooted, lambda = 10, model = "correlated", calibration = SP_calibration, control = chronos.control())
print("discrete")
SP_tree_1d <- chronos(SP_tree_rooted, lambda = 1, model = "discrete", calibration = SP_calibration, control = chronos.control())
SP_tree_10d <- chronos(SP_tree_rooted, lambda = 10, model = "discrete", calibration = SP_calibration, control = chronos.control())
print("relaxed")
SP_tree_1r <- chronos(SP_tree_rooted, lambda = 1, model = "relaxed", calibration = SP_calibration, control = chronos.control())
SP_tree_10r <- chronos(SP_tree_rooted, lambda = 10, model = "relaxed", calibration = SP_calibration, control = chronos.control())

# save best tree
write.tree(SP_tree_1 , file="1314_tree_1.tre")
write.tree(SP_tree_10 , file="1314_tree_10.tre")
write.tree(SP_tree_1d , file="1314_tree_1d.tre")
write.tree(SP_tree_10d , file="1314_tree_10d.tre")
write.tree(SP_tree_1r , file="1314_tree_1r.tre")
write.tree(SP_tree_10r , file="1314_tree_10r.tre")

### SE
SE_tree <- read.tree('1282.dnd')

# load original mash sketch and remove outliers
outliers_SE_tree <- read_tsv("1282_dist.tab", col_names = c("Ref_ID", "Query_ID", "mash_dist", "p_value", "matching_hashes")) %>%
  filter(mash_dist > 0.1) %>%
  separate(col=Query_ID, into=c("taxid", "id"), sep="/") %>%
  separate(col=id, into=c("id", "ext"), sep=".fna") %>%
  dplyr::select(id)

SE_tree_pruned <- drop.tip(SE_tree, outliers_SE_tree$id)

SE_tree_rooted <- root(SE_tree_pruned, outgroup="GCF_000013425.1_ASM1342v1_genomic", resolve.root = TRUE)

SE_calibration <- makeChronosCalib(SE_tree_rooted, node="root")

print("Staph epidermidis trees")

print("correlated")
SE_tree_1 <- chronos(SE_tree_rooted, lambda = 1, model = "correlated", calibration = SE_calibration, control = chronos.control())
SE_tree_10 <- chronos(SE_tree_rooted, lambda = 10, model = "correlated", calibration = SE_calibration, control = chronos.control())
print("discrete")
SE_tree_1d <- chronos(SE_tree_rooted, lambda = 1, model = "discrete", calibration = SE_calibration, control = chronos.control())
SE_tree_10d <- chronos(SE_tree_rooted, lambda = 10, model = "discrete", calibration = SE_calibration, control = chronos.control())
print("relaxed")
SE_tree_1r <- chronos(SE_tree_rooted, lambda = 1, model = "relaxed", calibration = SE_calibration, control = chronos.control())
SE_tree_10r <- chronos(SE_tree_rooted, lambda = 10, model = "relaxed", calibration = SE_calibration, control = chronos.control())

# save best tree
write.tree(SE_tree_1 , file="1282_tree_1.tre")
write.tree(SE_tree_10 , file="1282_tree_10.tre")
write.tree(SE_tree_1d , file="1282_tree_1d.tre")
write.tree(SE_tree_10d , file="1282_tree_10d.tre")
write.tree(SE_tree_1r , file="1282_tree_1r.tre")
write.tree(SE_tree_10r , file="1282_tree_10r.tre")

#### Efm

Efm_tree <- read.tree('1352.dnd')

# load original mash sketch and remove outliers
outliers_Efm_tree <- read_tsv("1352_dist.tab", col_names = c("Ref_ID", "Query_ID", "mash_dist", "p_value", "matching_hashes")) %>%
  filter(mash_dist > 0.1) %>%
  separate(col=Query_ID, into=c("taxid", "id"), sep="/") %>%
  separate(col=id, into=c("id", "ext"), sep=".fna") %>%
  dplyr::select(id)

Efm_tree_pruned <- drop.tip(Efm_tree, outliers_Efm_tree$id)

Efm_tree_rooted <- root(Efm_tree_pruned, outgroup="GCF_000393015.1_Ente_faec_T5_V1_genomic", resolve.root = TRUE)

Efm_calibration <- makeChronosCalib(Efm_tree_rooted, node="root")

print("E. faecium trees")

print("correlated")
Efm_tree_1 <- chronos(Efm_tree_rooted, lambda = 1, model = "correlated", calibration = Efm_calibration, control = chronos.control())
Efm_tree_10 <- chronos(Efm_tree_rooted, lambda = 10, model = "correlated", calibration = Efm_calibration, control = chronos.control())
print("discrete")
Efm_tree_1d <- chronos(Efm_tree_rooted, lambda = 1, model = "discrete", calibration = Efm_calibration, control = chronos.control())
Efm_tree_10d <- chronos(Efm_tree_rooted, lambda = 10, model = "discrete", calibration = Efm_calibration, control = chronos.control())
print("relaxed")
Efm_tree_1r <- chronos(Efm_tree_rooted, lambda = 1, model = "relaxed", calibration = Efm_calibration, control = chronos.control())
Efm_tree_10r <- chronos(Efm_tree_rooted, lambda = 10, model = "relaxed", calibration = Efm_calibration, control = chronos.control())

# save best tree
write.tree(Efm_tree_1 , file="1352_tree_1.tre")
write.tree(Efm_tree_10 , file="1352_tree_10.tre")
write.tree(Efm_tree_1d , file="1352_tree_1d.tre")
write.tree(Efm_tree_1d , file="1352_tree_10d.tre")
write.tree(Efm_tree_1r , file="1352_tree_1r.tre")
write.tree(Efm_tree_10r , file="1352_tree_10r.tre")

# #### Efs

Efs_tree <- read.tree('1351.dnd')

# load original mash sketch and remove outliers
outliers_Efs_tree <- read_tsv("1351_dist.tab", col_names = c("Ref_ID", "Query_ID", "mash_dist", "p_value", "matching_hashes")) %>%
  filter(mash_dist > 0.1) %>%
  separate(col=Query_ID, into=c("taxid", "id"), sep="/") %>%
  separate(col=id, into=c("id", "ext"), sep=".fna") %>%
  dplyr::select(id)

Efs_tree_pruned <- drop.tip(Efs_tree, outliers_Efs_tree$id)

Efs_tree_rooted <- root(Efs_tree_pruned, outgroup="GCF_010120755.1_ASM1012075v1_genomic", resolve.root = TRUE)

Efs_calibration <- makeChronosCalib(Efs_tree_rooted, node="root")

print("E. faecalis trees")

print("correlated")
Efs_tree_1 <- chronos(Efs_tree_rooted, lambda = 1, model = "correlated", calibration = Efs_calibration, control = chronos.control())
Efs_tree_10 <- chronos(Efs_tree_rooted, lambda = 10, model = "correlated", calibration = Efs_calibration, control = chronos.control())
print("discrete")
Efs_tree_1d <- chronos(Efs_tree_rooted, lambda = 1, model = "discrete", calibration = Efs_calibration, control = chronos.control())
Efs_tree_10d <- chronos(Efs_tree_rooted, lambda = 10, model = "discrete", calibration = Efs_calibration, control = chronos.control())
print("relaxed")
Efs_tree_1r <- chronos(Efs_tree_rooted, lambda = 1, model = "relaxed", calibration = Efs_calibration, control = chronos.control())
Efs_tree_10r <- chronos(Efs_tree_rooted, lambda = 10, model = "relaxed", calibration = Efs_calibration, control = chronos.control())

# save best tree
write.tree(Efs_tree_1 , file="1351_tree_1.tre")
write.tree(Efs_tree_10 , file="1351_tree_10.tre")
write.tree(Efs_tree_1d , file="1351_tree_1d.tre")
write.tree(Efs_tree_10d , file="1351_tree_10d.tre")
write.tree(Efs_tree_1r , file="1351_tree_1r.tre")
write.tree(Efs_tree_10r , file="1351_tree_10r.tre")

##### PA

PA_tree <- read.tree('287.dnd')

# load original mash sketch and remove outliers
outliers_PA_tree <- read_tsv("287_dist.tab", col_names = c("Ref_ID", "Query_ID", "mash_dist", "p_value", "matching_hashes")) %>%
  filter(mash_dist > 0.1) %>%
  separate(col=Query_ID, into=c("taxid", "id"), sep="/") %>%
  separate(col=id, into=c("id", "ext"), sep=".fna") %>%
  dplyr::select(id)

PA_tree_pruned <- drop.tip(PA_tree, outliers_PA_tree$id)

PA_tree_rooted <- root(PA_tree_pruned, outgroup="GCF_002220155.1_ASM222015v1_genomic", resolve.root = TRUE)

PA_calibration <- makeChronosCalib(PA_tree_rooted, node="root")

print("P. aeruginosa trees")

print("correlated")
PA_tree_1 <- chronos(PA_tree_rooted, lambda = 1, model = "correlated", calibration = PA_calibration, control = chronos.control())
PA_tree_10 <- chronos(PA_tree_rooted, lambda = 10, model = "correlated", calibration = PA_calibration, control = chronos.control())
print("discrete")
PA_tree_1d <- chronos(PA_tree_rooted, lambda = 1, model = "discrete", calibration = PA_calibration, control = chronos.control())
PA_tree_10d <- chronos(PA_tree_rooted, lambda = 10, model = "discrete", calibration = PA_calibration, control = chronos.control())
print("relaxed")
PA_tree_1r <- chronos(PA_tree_rooted, lambda = 1, model = "relaxed", calibration = PA_calibration, control = chronos.control())
PA_tree_10r <- chronos(PA_tree_rooted, lambda = 10, model = "relaxed", calibration = PA_calibration, control = chronos.control())

# save best tree
write.tree(PA_tree_1 , file="287_tree_1.tre")
write.tree(PA_tree_10 , file="287_tree_10.tre")
write.tree(PA_tree_1d , file="287_tree_1d.tre")
write.tree(PA_tree_10d , file="287_tree_10d.tre")
write.tree(PA_tree_1r , file="287_tree_1r.tre")
write.tree(PA_tree_10r , file="287_tree_10r.tre")

#### AB

AB_tree <- read.tree('470.dnd')

# load original mash sketch and remove outliers
outliers_AB_tree <- read_tsv("470_dist.tab", col_names = c("Ref_ID", "Query_ID", "mash_dist", "p_value", "matching_hashes")) %>%
  filter(mash_dist > 0.1) %>%
  separate(col=Query_ID, into=c("taxid", "id"), sep="/") %>%
  separate(col=id, into=c("id", "ext"), sep=".fna") %>%
  dplyr::select(id)

AB_tree_pruned <- drop.tip(AB_tree, outliers_AB_tree$id)

AB_tree_rooted <- root(AB_tree_pruned, outgroup="GCF_000092265.1_ASM9226v1_genomic", resolve.root = TRUE)

AB_calibration <- makeChronosCalib(AB_tree_rooted, node="root")

print("A. baumannii trees")

print("correlated")
AB_tree_1 <- chronos(AB_tree_rooted, lambda = 1, model = "correlated", calibration = AB_calibration, control = chronos.control())
AB_tree_10 <- chronos(AB_tree_rooted, lambda = 10, model = "correlated", calibration = AB_calibration, control = chronos.control())
print("discrete")
AB_tree_1d <- chronos(AB_tree_rooted, lambda = 1, model = "discrete", calibration = AB_calibration, control = chronos.control())
AB_tree_10d <- chronos(AB_tree_rooted, lambda = 10, model = "discrete", calibration = AB_calibration, control = chronos.control())
print("relaxed")
AB_tree_1r <- chronos(AB_tree_rooted, lambda = 1, model = "relaxed", calibration = AB_calibration, control = chronos.control())
AB_tree_10r <- chronos(AB_tree_rooted, lambda = 10, model = "relaxed", calibration = AB_calibration, control = chronos.control())

# save best tree
write.tree(AB_tree_1 , file="470_tree_1.tre")
write.tree(AB_tree_10 , file="470_tree_10.tre")
write.tree(AB_tree_1d , file="470_tree_1d.tre")
write.tree(AB_tree_10d , file="470_tree_10d.tre")
write.tree(AB_tree_1r , file="470_tree_1r.tre")
write.tree(AB_tree_10r , file="470_tree_10r.tre")

#### SA

SA_tree <- read.tree('1280.dnd')

# load original mash sketch and remove outliers
outliers_SA_tree <- read_tsv("1280_dist.tab", col_names = c("Ref_ID", "Query_ID", "mash_dist", "p_value", "matching_hashes")) %>%
  filter(mash_dist > 0.1) %>%
  separate(col=Query_ID, into=c("taxid", "id"), sep="/") %>%
  separate(col=id, into=c("id", "ext"), sep=".fna") %>%
  dplyr::select(id)

print(outliers_SA_tree)

print(outliers_SA_tree)

SA_tree_pruned <- drop.tip(SA_tree, outliers_SA_tree$id)

SA_tree_rooted <- root(SA_tree_pruned, outgroup="GCF_016028795.1_ASM1602879v1_genomic", resolve.root = TRUE)

SA_calibration <- makeChronosCalib(SA_tree_rooted, node="root")

print("Staph aureus trees")

print("correlated")
SA_tree_1 <- chronos(SA_tree_rooted, lambda = 1, model = "correlated", calibration = SA_calibration, control = chronos.control())
SA_tree_10 <- chronos(SA_tree_rooted, lambda = 10, model = "correlated", calibration = SA_calibration, control = chronos.control())
print("discrete")
SA_tree_1d <- chronos(SA_tree_rooted, lambda = 1, model = "discrete", calibration = SA_calibration, control = chronos.control())
SA_tree_10d <- chronos(SA_tree_rooted, lambda = 10, model = "discrete", calibration = SA_calibration, control = chronos.control())
print("relaxed")
SA_tree_1r <- chronos(SA_tree_rooted, lambda = 1, model = "relaxed", calibration = SA_calibration, control = chronos.control())
SA_tree_10r <- chronos(SA_tree_rooted, lambda = 10, model = "relaxed", calibration = SA_calibration, control = chronos.control())

# save best tree
write.tree(SA_tree_1 , file="1280_tree_1.tre")
write.tree(SA_tree_10 , file="1280_tree_10.tre")
write.tree(SA_tree_1d , file="1280_tree_1d.tre")
write.tree(SA_tree_10d , file="1280_tree_10d.tre")
write.tree(SA_tree_1r , file="1280_tree_1r.tre")
write.tree(SA_tree_10r , file="1280_tree_10r.tre")