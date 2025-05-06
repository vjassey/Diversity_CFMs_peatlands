# Diversity, abundance, and biogeography of CO2-fixing microorganisms in peatlands
# Marie Le Geay, Kyle Mayers, Anna Sytiuk, Ellen Dorrepaal, Martin Küttim, 
# Mariusz Lamentowicz, Eeva-Stiina Tuittila, Béatrice Lauga, and Vincent E.J. Jassey

# Author: Marie Le Geay and Vincent Jassey

################################################################################
# 4. Exploration of ASVs

# Main text:
## Figure 4

# Supplementary Information:
## Figures S1, S11, S12, S13 and S14

################################################################################
# |   |  Load packages ----
library(ggplot2); library(plyr); library(reshape2); library(grid); library(gridExtra); library(ape); library(picante)
library(seqinr); library(phytools); library(scales); library(nlme); library(vegan); library(RColorBrewer)
library("BiocManager"); library("Biobase"); library("rhdf5"); library("devtools"); library(phyloseq); library(msa)
library(ggtree); library(tidyjson); library(RColorBrewer);library(tidyverse); library(viridis); library(randomForest); library(caret); library(cluster)
library(tidyverse); library(cluster);  library(parallelDist); library(boral);  library(corrplot); library(igraph); library(car); library(fossil); library(ROCR); library(doMC)
library(effectsize); library(phyloseq); library(Biostrings); library(dplyr); library(tidyr); library(bioseq)
library(vegan); library(ade4); library(BAT); library(patchwork); library(hrbrthemes); library(gcookbook)
library(htmltools); library(HH) ;library(ggupset) ;library(UpSetR)

# |   |  Load data ----
setwd("C:/Users/mlegeay/OneDrive/Documents/chapter_2/Redaction/Révision/RStudio-to-publish/Data")# set up a new path

# Metadata ----
metadata_16S <- read.table("meta_16S.csv", header = T, sep = ";")
metadata_23S <- read.table("meta_23S.csv", header = T, sep = ";")
metadata_cbbL <- read.table("meta_cbbL.csv", header = T, sep = ";")
metadata_bchY <- read.table("meta_bchY.csv", header = T, sep = ";")

# ASV matrices ----
ASV_23S <- read.table("t_ASV_23S_HP.csv", header = T, sep = ",")
rownames(ASV_23S) <- ASV_23S$X
ASV_23S <- ASV_23S[,-1]

ASV_cbbL <- read.table("t_ASV_cbbL.csv", header = T, sep = ",")
rownames(ASV_cbbL) <- ASV_cbbL$X
ASV_cbbL <- ASV_cbbL[,-1]

ASV_bchY <- read.table("t_ASV_bchY.csv", header = T, sep = ",")
rownames(ASV_bchY) <- ASV_bchY$X
ASV_bchY <- ASV_bchY[,-1]


# Taxonomy ----
Taxo_23S <- read.table("Taxo_23S.csv", header = T, sep = ",")
rownames(Taxo_23S) <- Taxo_23S$X
Taxo_23S <- Taxo_23S[,-1]

Taxo_cbbL <- read.table("Taxo_cbbL.csv", header = T, sep = ",")
rownames(Taxo_cbbL) <- Taxo_cbbL$X
Taxo_cbbL <- Taxo_cbbL[,-1]

Taxo_bchY <- read.table("Taxo_bchY.csv", header = T, sep = ",")
rownames(Taxo_bchY) <- Taxo_bchY$X
Taxo_bchY <- Taxo_bchY[,-1]

# 1. Rarefaction curves ----
# |  | 16S ----
data_rar_16S <- read.table("ASV_16S_RarCurve.csv", header = T, sep = ",")
rownames(data_rar_16S) <- data_rar_16S$X
data_rar_16S <- data_rar_16S[,-1]
names(metadata_16S)[1] <- "Site"
rar_16S <- rarecurve(data_rar_16S,
                     step = 50,
                     tidy = TRUE) %>%
  left_join(metadata_16S)

plot_16S <- ggplot(rar_16S) +
  geom_line(aes(x = Sample, y = Species, group = Site, 
                colour = as.factor(Location), linetype = as.factor(Depth))) +
  theme_classic()   
plot_16S

# |  | 23S ----
data_rar_23S <- read.table("ASV_23S_RarCurve.csv", header = T, sep = ",")
rownames(data_rar_23S) <- data_rar_23S$X
data_rar_23S <- data_rar_23S[,-1]
names(metadata_23S)[1] <- "Site"
rar_23S <- rarecurve(data_rar_23S,
                     step = 50,
                     tidy = TRUE) %>%
  left_join(metadata_23S)

plot_23S <- ggplot(rar_23S) +
  geom_line(aes(x = Sample, y = Species, group = Site, 
                colour = as.factor(Location), linetype = as.factor(Depth))) +
  theme_classic()   
plot_23S  

# |  | cbbL ----
data_rar_cbbL <- read.table("ASV_cbbL_RarCurve.csv", header = T, sep = ",")
rownames(data_rar_cbbL) <- data_rar_cbbL$X
data_rar_cbbL <- data_rar_cbbL[,-1]
names(metadata_cbbL)[1] <- "Site"
rar_cbbL <- rarecurve(data_rar_cbbL,
                      step = 50,
                      tidy = TRUE) %>%
  left_join(metadata_cbbL)

plot_cbbL <- ggplot(rar_cbbL) +
  geom_line(aes(x = Sample, y = Species, group = Site, 
                colour = Location, linetype = as.factor(Depth))) +
  theme_classic()   
plot_cbbL

# |  | bchY ----
data_rar_bchY <- read.table("ASV_bchY_RarCurve.csv", header = T, sep = ",")
rownames(data_rar_bchY) <- data_rar_bchY$X
data_rar_bchY <- data_rar_bchY[,-1]
names(metadata_bchY)[1] <- "Site"
rar_bchY <- rarecurve(data_rar_bchY,
                      step = 50,
                      tidy = TRUE) %>%
  left_join(metadata_bchY)

plot_bchY <- ggplot(rar_bchY) +
  geom_line(aes(x = Sample, y = Species, group = Site, 
                colour = as.factor(Location), linetype = as.factor(Depth))) +
  theme_classic()   
plot_bchY

# All graph together ----
p_all_rarefy <-  (plot_16S | plot_23S) / (plot_cbbL | plot_bchY)
p_all_rarefy



# 2. Heatmaps at class level ----
# | | 23S ---- 
asv_23S_loc_depth <- aggregate(ASV_23S, by = list(as.character(metadata_23S$Location_Depth)), FUN = sum)
rownames(asv_23S_loc_depth) <- asv_23S_loc_depth$Group.1
asv_23S_loc_depth <- asv_23S_loc_depth[,-1]

## 2 .aggregation by taxononic rank : 
asv_23S_loc_depth2 <- aggregate(as.data.frame(t(asv_23S_loc_depth)), by = list(as.character(Taxo_23S$Class)), FUN = sum)

# Relative abundance :
asv_23S_loc_depth2 <- subset(asv_23S_loc_depth2, 
                             asv_23S_loc_depth2$Group.1 != "Multi-affiliation")

asv_23S_loc_depth2$Abisko_1 <- asv_23S_loc_depth2$Abisko_1 *100 / sum (asv_23S_loc_depth2$Abisko_1)
asv_23S_loc_depth2$Abisko_2 <- asv_23S_loc_depth2$Abisko_2 *100 / sum (asv_23S_loc_depth2$Abisko_2)
asv_23S_loc_depth2$Abisko_3 <- asv_23S_loc_depth2$Abisko_3 *100 / sum (asv_23S_loc_depth2$Abisko_3)

asv_23S_loc_depth2$Counozouls_1 <- asv_23S_loc_depth2$Counozouls_1  *100 / sum (asv_23S_loc_depth2$Counozouls_1)
asv_23S_loc_depth2$Counozouls_2 <- asv_23S_loc_depth2$Counozouls_2  *100 / sum (asv_23S_loc_depth2$Counozouls_2)
asv_23S_loc_depth2$Counozouls_3 <- asv_23S_loc_depth2$Counozouls_3  *100 / sum (asv_23S_loc_depth2$Counozouls_3)

asv_23S_loc_depth2$Mann_1 <- asv_23S_loc_depth2$Mann_1 *100 / sum (asv_23S_loc_depth2$Mann_1)
asv_23S_loc_depth2$Mann_2 <- asv_23S_loc_depth2$Mann_2 *100 / sum (asv_23S_loc_depth2$Mann_2)
asv_23S_loc_depth2$Mann_3 <- asv_23S_loc_depth2$Mann_3 *100 / sum (asv_23S_loc_depth2$Mann_3)

asv_23S_loc_depth2$Siikaneva_1 <- asv_23S_loc_depth2$Siikaneva_1 *100 / sum (asv_23S_loc_depth2$Siikaneva_1)
asv_23S_loc_depth2$Siikaneva_2 <- asv_23S_loc_depth2$Siikaneva_2 *100 / sum (asv_23S_loc_depth2$Siikaneva_2)
asv_23S_loc_depth2$Siikaneva_3 <- asv_23S_loc_depth2$Siikaneva_3 *100 / sum (asv_23S_loc_depth2$Siikaneva_3)

asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Abisko_1 <= 0.05] <- "Others" # if relative abundance is too low, go on other
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Abisko_2 <= 0.05]<- "Others" # if relative abundance is too low, go on other
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Abisko_3 <= 0.05] <- "Others" # if relative abundance is too low, go on other

asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Counozouls_1 <= 0.05] <- "Others"
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Counozouls_2 <= 0.05] <- "Others"
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Counozouls_3 <= 0.05] <- "Others"

asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Mann_1 <= 0.05] <- "Others"
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Mann_2 <= 0.05] <- "Others"
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Mann_3 <= 0.05] <- "Others"

asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Siikaneva_1 <= 0.05] <- "Others"
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Siikaneva_2 <= 0.05] <- "Others"
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Siikaneva_3 <= 0.05] <- "Others"

#asv_23S_loc_depth2$Total_1 <-  rowSums(asv_23S_agg_loc2[,2:ncol(asv_23S_agg_loc2)]) *100 / sum(asv_23S_agg_loc2[,2:ncol(asv_23S_agg_loc2)]) 

# reshape :
asv_23S_loc_depth2_melt <- reshape2::melt(asv_23S_loc_depth2)

# plot :
asv_23S_loc_depth2_melt <- as.data.frame(asv_23S_loc_depth2_melt)
asv_23S_loc_depth2_melt$variable

asv_23S_loc_depth2_melt$variable <- fct_relevel(asv_23S_loc_depth2_melt$variable, 
                                                c("Counozouls_1","Counozouls_2","Counozouls_3",
                                                  "Mann_1","Mann_2","Mann_3",
                                                  "Siikaneva_1","Siikaneva_2","Siikaneva_3",
                                                  "Abisko_1","Abisko_2","Abisko_3"))
S23_heatmap <- ggplot(asv_23S_loc_depth2_melt, 
                      aes(x = variable , y = Group.1, fill = value)) +
  geom_tile(color = "black",
            #lwd = 1.5,
            #linetype = 1
  ) +
  scale_fill_gradient2(low = "white", 
                       high = "#3E938BFF") +
  coord_fixed()
S23_heatmap


# | |  cbbL ----
asv_cbbL_loc_depth <- aggregate(ASV_cbbL, by = list(as.character(metadata_cbbL$Location_Depth)), FUN = sum)
rownames(asv_cbbL_loc_depth) <- asv_cbbL_loc_depth$Group.1
asv_cbbL_loc_depth <- asv_cbbL_loc_depth[,-1]

## 2 .aggregation by taxononic rank : Class
asv_cbbL_loc_depth2 <- aggregate(as.data.frame(t(asv_cbbL_loc_depth)), by = list(as.character(Taxo_cbbL$Class)), FUN = sum)

# Relative abundance :
asv_cbbL_loc_depth2 <- subset(asv_cbbL_loc_depth2, 
                              asv_cbbL_loc_depth2$Group.1 != "Multi-affiliation")

asv_cbbL_loc_depth2$Abisko_1 <- asv_cbbL_loc_depth2$Abisko_1 *100 / sum (asv_cbbL_loc_depth2$Abisko_1)
asv_cbbL_loc_depth2$Abisko_2 <- asv_cbbL_loc_depth2$Abisko_2 *100 / sum (asv_cbbL_loc_depth2$Abisko_2)
asv_cbbL_loc_depth2$Abisko_3 <- asv_cbbL_loc_depth2$Abisko_3 *100 / sum (asv_cbbL_loc_depth2$Abisko_3)

asv_cbbL_loc_depth2$Counozouls_1 <- asv_cbbL_loc_depth2$Counozouls_1  *100 / sum (asv_cbbL_loc_depth2$Counozouls_1)
asv_cbbL_loc_depth2$Counozouls_2 <- asv_cbbL_loc_depth2$Counozouls_2  *100 / sum (asv_cbbL_loc_depth2$Counozouls_2)
asv_cbbL_loc_depth2$Counozouls_3 <- asv_cbbL_loc_depth2$Counozouls_3  *100 / sum (asv_cbbL_loc_depth2$Counozouls_3)

asv_cbbL_loc_depth2$Mann_1 <- asv_cbbL_loc_depth2$Mann_1 *100 / sum (asv_cbbL_loc_depth2$Mann_1)
asv_cbbL_loc_depth2$Mann_2 <- asv_cbbL_loc_depth2$Mann_2 *100 / sum (asv_cbbL_loc_depth2$Mann_2)
asv_cbbL_loc_depth2$Mann_3 <- asv_cbbL_loc_depth2$Mann_3 *100 / sum (asv_cbbL_loc_depth2$Mann_3)

asv_cbbL_loc_depth2$Siikaneva_1 <- asv_cbbL_loc_depth2$Siikaneva_1 *100 / sum (asv_cbbL_loc_depth2$Siikaneva_1)
asv_cbbL_loc_depth2$Siikaneva_2 <- asv_cbbL_loc_depth2$Siikaneva_2 *100 / sum (asv_cbbL_loc_depth2$Siikaneva_2)
asv_cbbL_loc_depth2$Siikaneva_3 <- asv_cbbL_loc_depth2$Siikaneva_3 *100 / sum (asv_cbbL_loc_depth2$Siikaneva_3)

asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Abisko_1 <= 0.05] <- "Others" # if relative abundance is too low, go on other
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Abisko_2 <= 0.05]<- "Others" # if relative abundance is too low, go on other
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Abisko_3 <= 0.05] <- "Others" # if relative abundance is too low, go on other

asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Counozouls_1 <= 0.05] <- "Others"
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Counozouls_2 <= 0.05] <- "Others"
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Counozouls_3 <= 0.05] <- "Others"

asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Mann_1 <= 0.05] <- "Others"
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Mann_2 <= 0.05] <- "Others"
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Mann_3 <= 0.05] <- "Others"

asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Siikaneva_1 <= 0.05] <- "Others"
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Siikaneva_2 <= 0.05] <- "Others"
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Siikaneva_3 <= 0.05] <- "Others"

# reshape :
asv_cbbL_loc_depth2_melt <- reshape2::melt(asv_cbbL_loc_depth2)

# plot :
asv_cbbL_loc_depth2_melt <- as.data.frame(asv_cbbL_loc_depth2_melt)
asv_cbbL_loc_depth2_melt$variable
asv_cbbL_loc_depth2_melt$variable <- fct_relevel(asv_cbbL_loc_depth2_melt$variable, 
                                                 c("Counozouls_1","Counozouls_2","Counozouls_3",
                                                   "Mann_1","Mann_2","Mann_3",
                                                   "Siikaneva_1","Siikaneva_2","Siikaneva_3",
                                                   "Abisko_1","Abisko_2","Abisko_3"))

cbbL_heatmap <- ggplot(asv_cbbL_loc_depth2_melt, 
                       aes(x = variable , y = Group.1, fill = value)) +
  geom_tile(color = "black",
            #lwd = 1.5,
            #linetype = 1
  ) +
  scale_fill_gradient2(low = "white", 
                       high = "#3E938BFF") +
  coord_fixed()
cbbL_heatmap 


#  | | bchY ----
asv_bchY_loc_depth <- aggregate(ASV_bchY, by = list(as.character(metadata_bchY$Location_Depth)), FUN = sum)
rownames(asv_bchY_loc_depth) <- asv_bchY_loc_depth$Group.1
asv_bchY_loc_depth <- asv_bchY_loc_depth[,-1]

## 2 .aggregation by taxononic rank : Class
asv_bchY_loc_depth2 <- aggregate(as.data.frame(t(asv_bchY_loc_depth)), by = list(as.character(Taxo_bchY$Class)), FUN = sum)

# Relative abundance :
asv_bchY_loc_depth2 <- subset(asv_bchY_loc_depth2, 
                              asv_bchY_loc_depth2$Group.1 != "Multi-affiliation")

asv_bchY_loc_depth2$Abisko_1 <- asv_bchY_loc_depth2$Abisko_1 *100 / sum (asv_bchY_loc_depth2$Abisko_1)
asv_bchY_loc_depth2$Abisko_2 <- asv_bchY_loc_depth2$Abisko_2 *100 / sum (asv_bchY_loc_depth2$Abisko_2)
asv_bchY_loc_depth2$Abisko_3 <- asv_bchY_loc_depth2$Abisko_3 *100 / sum (asv_bchY_loc_depth2$Abisko_3)

asv_bchY_loc_depth2$Counozouls_1 <- asv_bchY_loc_depth2$Counozouls_1  *100 / sum (asv_bchY_loc_depth2$Counozouls_1)
asv_bchY_loc_depth2$Counozouls_2 <- asv_bchY_loc_depth2$Counozouls_2  *100 / sum (asv_bchY_loc_depth2$Counozouls_2)
asv_bchY_loc_depth2$Counozouls_3 <- asv_bchY_loc_depth2$Counozouls_3  *100 / sum (asv_bchY_loc_depth2$Counozouls_3)

asv_bchY_loc_depth2$Mann_1 <- asv_bchY_loc_depth2$Mann_1 *100 / sum (asv_bchY_loc_depth2$Mann_1)
asv_bchY_loc_depth2$Mann_2 <- asv_bchY_loc_depth2$Mann_2 *100 / sum (asv_bchY_loc_depth2$Mann_2)
asv_bchY_loc_depth2$Mann_3 <- asv_bchY_loc_depth2$Mann_3 *100 / sum (asv_bchY_loc_depth2$Mann_3)

asv_bchY_loc_depth2$Siikaneva_1 <- asv_bchY_loc_depth2$Siikaneva_1 *100 / sum (asv_bchY_loc_depth2$Siikaneva_1)
asv_bchY_loc_depth2$Siikaneva_2 <- asv_bchY_loc_depth2$Siikaneva_2 *100 / sum (asv_bchY_loc_depth2$Siikaneva_2)
asv_bchY_loc_depth2$Siikaneva_3 <- asv_bchY_loc_depth2$Siikaneva_3 *100 / sum (asv_bchY_loc_depth2$Siikaneva_3)

asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Abisko_1 <= 0.05] <- "Others" # if relative abundance is too low, go on other
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Abisko_2 <= 0.05]<- "Others" # if relative abundance is too low, go on other
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Abisko_3 <= 0.05] <- "Others" # if relative abundance is too low, go on other

asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Counozouls_1 <= 0.05] <- "Others"
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Counozouls_2 <= 0.05] <- "Others"
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Counozouls_3 <= 0.05] <- "Others"

asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Mann_1 <= 0.05] <- "Others"
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Mann_2 <= 0.05] <- "Others"
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Mann_3 <= 0.05] <- "Others"

asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Siikaneva_1 <= 0.05] <- "Others"
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Siikaneva_2 <= 0.05] <- "Others"
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Siikaneva_3 <= 0.05] <- "Others"

# reshape :
asv_bchY_loc_depth2_melt <- reshape2::melt(asv_bchY_loc_depth2)

# plot :
asv_bchY_loc_depth2_melt <- as.data.frame(asv_bchY_loc_depth2_melt)
asv_bchY_loc_depth2_melt$variable

asv_bchY_loc_depth2_melt$variable <- fct_relevel(asv_bchY_loc_depth2_melt$variable, 
                                                 c("Counozouls_1","Counozouls_2","Counozouls_3",
                                                   "Mann_1","Mann_2","Mann_3",
                                                   "Siikaneva_1","Siikaneva_2","Siikaneva_3",
                                                   "Abisko_1","Abisko_2","Abisko_3"))


bchY_heatmap <- ggplot(asv_bchY_loc_depth2_melt, 
                       aes(x = variable , y = Group.1, fill = value)) +
  geom_tile(color = "black",
            #lwd = 1.5,
            #linetype = 1
  ) +
  scale_fill_gradient2(low = "white", 
                       high = "#3E938BFF") +
  coord_fixed()
bchY_heatmap


# All graph together
heatmap_class <-  S23_heatmap / cbbL_heatmap / bchY_heatmap
heatmap_class


# 3. Heatmaps at family level ----
# | | 23S ---- 
asv_23S_loc_depth <- aggregate(ASV_23S, by = list(as.character(metadata_23S$Location_Depth)), FUN = sum)
rownames(asv_23S_loc_depth) <- asv_23S_loc_depth$Group.1
asv_23S_loc_depth <- asv_23S_loc_depth[,-1]

## 2 .aggregation by taxononic rank : 
asv_23S_loc_depth2 <- aggregate(as.data.frame(t(asv_23S_loc_depth)), by = list(as.character(Taxo_23S$Family)), FUN = sum)

# Relative abundance :
asv_23S_loc_depth2 <- subset(asv_23S_loc_depth2, 
                             asv_23S_loc_depth2$Group.1 != "Multi-affiliation")

asv_23S_loc_depth2$Abisko_1 <- asv_23S_loc_depth2$Abisko_1 *100 / sum (asv_23S_loc_depth2$Abisko_1)
asv_23S_loc_depth2$Abisko_2 <- asv_23S_loc_depth2$Abisko_2 *100 / sum (asv_23S_loc_depth2$Abisko_2)
asv_23S_loc_depth2$Abisko_3 <- asv_23S_loc_depth2$Abisko_3 *100 / sum (asv_23S_loc_depth2$Abisko_3)

asv_23S_loc_depth2$Counozouls_1 <- asv_23S_loc_depth2$Counozouls_1  *100 / sum (asv_23S_loc_depth2$Counozouls_1)
asv_23S_loc_depth2$Counozouls_2 <- asv_23S_loc_depth2$Counozouls_2  *100 / sum (asv_23S_loc_depth2$Counozouls_2)
asv_23S_loc_depth2$Counozouls_3 <- asv_23S_loc_depth2$Counozouls_3  *100 / sum (asv_23S_loc_depth2$Counozouls_3)

asv_23S_loc_depth2$Mann_1 <- asv_23S_loc_depth2$Mann_1 *100 / sum (asv_23S_loc_depth2$Mann_1)
asv_23S_loc_depth2$Mann_2 <- asv_23S_loc_depth2$Mann_2 *100 / sum (asv_23S_loc_depth2$Mann_2)
asv_23S_loc_depth2$Mann_3 <- asv_23S_loc_depth2$Mann_3 *100 / sum (asv_23S_loc_depth2$Mann_3)

asv_23S_loc_depth2$Siikaneva_1 <- asv_23S_loc_depth2$Siikaneva_1 *100 / sum (asv_23S_loc_depth2$Siikaneva_1)
asv_23S_loc_depth2$Siikaneva_2 <- asv_23S_loc_depth2$Siikaneva_2 *100 / sum (asv_23S_loc_depth2$Siikaneva_2)
asv_23S_loc_depth2$Siikaneva_3 <- asv_23S_loc_depth2$Siikaneva_3 *100 / sum (asv_23S_loc_depth2$Siikaneva_3)

asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Abisko_1 <= 0.05] <- "Others" # if relative abundance is too low, go on other
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Abisko_2 <= 0.05]<- "Others" # if relative abundance is too low, go on other
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Abisko_3 <= 0.05] <- "Others" # if relative abundance is too low, go on other

asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Counozouls_1 <= 0.05] <- "Others"
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Counozouls_2 <= 0.05] <- "Others"
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Counozouls_3 <= 0.05] <- "Others"

asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Mann_1 <= 0.05] <- "Others"
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Mann_2 <= 0.05] <- "Others"
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Mann_3 <= 0.05] <- "Others"

asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Siikaneva_1 <= 0.05] <- "Others"
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Siikaneva_2 <= 0.05] <- "Others"
asv_23S_loc_depth2$Group.1[asv_23S_loc_depth2$Siikaneva_3 <= 0.05] <- "Others"

#asv_23S_loc_depth2$Total_1 <-  rowSums(asv_23S_agg_loc2[,2:ncol(asv_23S_agg_loc2)]) *100 / sum(asv_23S_agg_loc2[,2:ncol(asv_23S_agg_loc2)]) 

# reshape :
asv_23S_loc_depth2_melt <- reshape2::melt(asv_23S_loc_depth2)

# plot :
asv_23S_loc_depth2_melt <- as.data.frame(asv_23S_loc_depth2_melt)
asv_23S_loc_depth2_melt$variable

asv_23S_loc_depth2_melt$variable <- fct_relevel(asv_23S_loc_depth2_melt$variable, 
                                                c("Counozouls_1","Counozouls_2","Counozouls_3",
                                                  "Mann_1","Mann_2","Mann_3",
                                                  "Siikaneva_1","Siikaneva_2","Siikaneva_3",
                                                  "Abisko_1","Abisko_2","Abisko_3"))
S23_heatmap_fam <- ggplot(asv_23S_loc_depth2_melt, 
                      aes(x = variable , y = Group.1, fill = value)) +
  geom_tile(color = "black",
            #lwd = 1.5,
            #linetype = 1
  ) +
  scale_fill_gradient2(low = "white", 
                       high = "#3E938BFF") +
  coord_fixed()
S23_heatmap_fam


# | |  cbbL ----
asv_cbbL_loc_depth <- aggregate(ASV_cbbL, by = list(as.character(metadata_cbbL$Location_Depth)), FUN = sum)
rownames(asv_cbbL_loc_depth) <- asv_cbbL_loc_depth$Group.1
asv_cbbL_loc_depth <- asv_cbbL_loc_depth[,-1]

## 2 .aggregation by taxononic rank : Class
asv_cbbL_loc_depth2 <- aggregate(as.data.frame(t(asv_cbbL_loc_depth)), by = list(as.character(Taxo_cbbL$Family)), FUN = sum)

# Relative abundance :
asv_cbbL_loc_depth2 <- subset(asv_cbbL_loc_depth2, 
                              asv_cbbL_loc_depth2$Group.1 != "Multi-affiliation")

asv_cbbL_loc_depth2$Abisko_1 <- asv_cbbL_loc_depth2$Abisko_1 *100 / sum (asv_cbbL_loc_depth2$Abisko_1)
asv_cbbL_loc_depth2$Abisko_2 <- asv_cbbL_loc_depth2$Abisko_2 *100 / sum (asv_cbbL_loc_depth2$Abisko_2)
asv_cbbL_loc_depth2$Abisko_3 <- asv_cbbL_loc_depth2$Abisko_3 *100 / sum (asv_cbbL_loc_depth2$Abisko_3)

asv_cbbL_loc_depth2$Counozouls_1 <- asv_cbbL_loc_depth2$Counozouls_1  *100 / sum (asv_cbbL_loc_depth2$Counozouls_1)
asv_cbbL_loc_depth2$Counozouls_2 <- asv_cbbL_loc_depth2$Counozouls_2  *100 / sum (asv_cbbL_loc_depth2$Counozouls_2)
asv_cbbL_loc_depth2$Counozouls_3 <- asv_cbbL_loc_depth2$Counozouls_3  *100 / sum (asv_cbbL_loc_depth2$Counozouls_3)

asv_cbbL_loc_depth2$Mann_1 <- asv_cbbL_loc_depth2$Mann_1 *100 / sum (asv_cbbL_loc_depth2$Mann_1)
asv_cbbL_loc_depth2$Mann_2 <- asv_cbbL_loc_depth2$Mann_2 *100 / sum (asv_cbbL_loc_depth2$Mann_2)
asv_cbbL_loc_depth2$Mann_3 <- asv_cbbL_loc_depth2$Mann_3 *100 / sum (asv_cbbL_loc_depth2$Mann_3)

asv_cbbL_loc_depth2$Siikaneva_1 <- asv_cbbL_loc_depth2$Siikaneva_1 *100 / sum (asv_cbbL_loc_depth2$Siikaneva_1)
asv_cbbL_loc_depth2$Siikaneva_2 <- asv_cbbL_loc_depth2$Siikaneva_2 *100 / sum (asv_cbbL_loc_depth2$Siikaneva_2)
asv_cbbL_loc_depth2$Siikaneva_3 <- asv_cbbL_loc_depth2$Siikaneva_3 *100 / sum (asv_cbbL_loc_depth2$Siikaneva_3)

asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Abisko_1 <= 0.05] <- "Others" # if relative abundance is too low, go on other
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Abisko_2 <= 0.05]<- "Others" # if relative abundance is too low, go on other
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Abisko_3 <= 0.05] <- "Others" # if relative abundance is too low, go on other

asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Counozouls_1 <= 0.05] <- "Others"
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Counozouls_2 <= 0.05] <- "Others"
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Counozouls_3 <= 0.05] <- "Others"

asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Mann_1 <= 0.05] <- "Others"
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Mann_2 <= 0.05] <- "Others"
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Mann_3 <= 0.05] <- "Others"

asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Siikaneva_1 <= 0.05] <- "Others"
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Siikaneva_2 <= 0.05] <- "Others"
asv_cbbL_loc_depth2$Group.1[asv_cbbL_loc_depth2$Siikaneva_3 <= 0.05] <- "Others"

# reshape :
asv_cbbL_loc_depth2_melt <- reshape2::melt(asv_cbbL_loc_depth2)

# plot :
asv_cbbL_loc_depth2_melt <- as.data.frame(asv_cbbL_loc_depth2_melt)
asv_cbbL_loc_depth2_melt$variable
asv_cbbL_loc_depth2_melt$variable <- fct_relevel(asv_cbbL_loc_depth2_melt$variable, 
                                                 c("Counozouls_1","Counozouls_2","Counozouls_3",
                                                   "Mann_1","Mann_2","Mann_3",
                                                   "Siikaneva_1","Siikaneva_2","Siikaneva_3",
                                                   "Abisko_1","Abisko_2","Abisko_3"))

cbbL_heatmap_fam <- ggplot(asv_cbbL_loc_depth2_melt, 
                       aes(x = variable , y = Group.1, fill = value)) +
  geom_tile(color = "black",
            #lwd = 1.5,
            #linetype = 1
  ) +
  scale_fill_gradient2(low = "white", 
                       high = "#3E938BFF") +
  coord_fixed()
cbbL_heatmap_fam 


#  | | bchY ----
asv_bchY_loc_depth <- aggregate(ASV_bchY, by = list(as.character(metadata_bchY$Location_Depth)), FUN = sum)
rownames(asv_bchY_loc_depth) <- asv_bchY_loc_depth$Group.1
asv_bchY_loc_depth <- asv_bchY_loc_depth[,-1]

## 2 .aggregation by taxononic rank : Class
asv_bchY_loc_depth2 <- aggregate(as.data.frame(t(asv_bchY_loc_depth)), by = list(as.character(Taxo_bchY$Family)), FUN = sum)

# Relative abundance :
asv_bchY_loc_depth2 <- subset(asv_bchY_loc_depth2, 
                              asv_bchY_loc_depth2$Group.1 != "Multi-affiliation")

asv_bchY_loc_depth2$Abisko_1 <- asv_bchY_loc_depth2$Abisko_1 *100 / sum (asv_bchY_loc_depth2$Abisko_1)
asv_bchY_loc_depth2$Abisko_2 <- asv_bchY_loc_depth2$Abisko_2 *100 / sum (asv_bchY_loc_depth2$Abisko_2)
asv_bchY_loc_depth2$Abisko_3 <- asv_bchY_loc_depth2$Abisko_3 *100 / sum (asv_bchY_loc_depth2$Abisko_3)

asv_bchY_loc_depth2$Counozouls_1 <- asv_bchY_loc_depth2$Counozouls_1  *100 / sum (asv_bchY_loc_depth2$Counozouls_1)
asv_bchY_loc_depth2$Counozouls_2 <- asv_bchY_loc_depth2$Counozouls_2  *100 / sum (asv_bchY_loc_depth2$Counozouls_2)
asv_bchY_loc_depth2$Counozouls_3 <- asv_bchY_loc_depth2$Counozouls_3  *100 / sum (asv_bchY_loc_depth2$Counozouls_3)

asv_bchY_loc_depth2$Mann_1 <- asv_bchY_loc_depth2$Mann_1 *100 / sum (asv_bchY_loc_depth2$Mann_1)
asv_bchY_loc_depth2$Mann_2 <- asv_bchY_loc_depth2$Mann_2 *100 / sum (asv_bchY_loc_depth2$Mann_2)
asv_bchY_loc_depth2$Mann_3 <- asv_bchY_loc_depth2$Mann_3 *100 / sum (asv_bchY_loc_depth2$Mann_3)

asv_bchY_loc_depth2$Siikaneva_1 <- asv_bchY_loc_depth2$Siikaneva_1 *100 / sum (asv_bchY_loc_depth2$Siikaneva_1)
asv_bchY_loc_depth2$Siikaneva_2 <- asv_bchY_loc_depth2$Siikaneva_2 *100 / sum (asv_bchY_loc_depth2$Siikaneva_2)
asv_bchY_loc_depth2$Siikaneva_3 <- asv_bchY_loc_depth2$Siikaneva_3 *100 / sum (asv_bchY_loc_depth2$Siikaneva_3)

asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Abisko_1 <= 0.05] <- "Others" # if relative abundance is too low, go on other
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Abisko_2 <= 0.05]<- "Others" # if relative abundance is too low, go on other
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Abisko_3 <= 0.05] <- "Others" # if relative abundance is too low, go on other

asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Counozouls_1 <= 0.05] <- "Others"
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Counozouls_2 <= 0.05] <- "Others"
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Counozouls_3 <= 0.05] <- "Others"

asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Mann_1 <= 0.05] <- "Others"
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Mann_2 <= 0.05] <- "Others"
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Mann_3 <= 0.05] <- "Others"

asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Siikaneva_1 <= 0.05] <- "Others"
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Siikaneva_2 <= 0.05] <- "Others"
asv_bchY_loc_depth2$Group.1[asv_bchY_loc_depth2$Siikaneva_3 <= 0.05] <- "Others"

# reshape :
asv_bchY_loc_depth2_melt <- reshape2::melt(asv_bchY_loc_depth2)

# plot :
asv_bchY_loc_depth2_melt <- as.data.frame(asv_bchY_loc_depth2_melt)
asv_bchY_loc_depth2_melt$variable

asv_bchY_loc_depth2_melt$variable <- fct_relevel(asv_bchY_loc_depth2_melt$variable, 
                                                 c("Counozouls_1","Counozouls_2","Counozouls_3",
                                                   "Mann_1","Mann_2","Mann_3",
                                                   "Siikaneva_1","Siikaneva_2","Siikaneva_3",
                                                   "Abisko_1","Abisko_2","Abisko_3"))


bchY_heatmap_fam <- ggplot(asv_bchY_loc_depth2_melt, 
                       aes(x = variable , y = Group.1, fill = value)) +
  geom_tile(color = "black",
            #lwd = 1.5,
            #linetype = 1
  ) +
  scale_fill_gradient2(low = "white", 
                       high = "#3E938BFF") +
  coord_fixed()
bchY_heatmap_fam


# All graph together ----
heatmap_family <-  S23_heatmap_fam / cbbL_heatmap_fam / bchY_heatmap_fam
heatmap_family

# 4. Add statistics for classes ----
# | | 23S ----
EF_heatmap_23S <- aggregate(ASV_23S, by = list(as.character(metadata_23S$ID_23S)), FUN = sum)
rownames(EF_heatmap_23S) <- EF_heatmap_23S$Group.1
EF_heatmap_23S <- EF_heatmap_23S[,-1]

# Data preparation:
# Aggregation by taxononic rank : 
EF_heatmap_23S_2 <- aggregate(as.data.frame(t(EF_heatmap_23S)), by = list(as.character(Taxo_23S$Class)), FUN = sum)

# Relative abundance :
EF_heatmap_23S_2   <- subset(EF_heatmap_23S_2  , 
                             EF_heatmap_23S_2$Group.1 != "Multi-affiliation")


# Creation of a phyloseq object :
EF_heatmap_23S_2_mat <- as.data.frame(EF_heatmap_23S_2)
rownames(EF_heatmap_23S_2_mat) <- EF_heatmap_23S_2_mat$Group.1
EF_heatmap_23S_2_mat<- EF_heatmap_23S_2_mat[,-1]

EF_heatmap_23S_OTU = otu_table(EF_heatmap_23S_2_mat, taxa_are_rows = TRUE)
EF_heatmap_23S_samples = sample_data(metadata_23S)
rownames(EF_heatmap_23S_samples) <- EF_heatmap_23S_samples$ID_23S
Phyloseq_23S <- phyloseq(EF_heatmap_23S_OTU, EF_heatmap_23S_samples)
Phyloseq_23S

Phyloseq_23S_relab <- transform_sample_counts(Phyloseq_23S, function(x) x/sum(x))

# merge samples
melted <- psmelt(Phyloseq_23S_relab)
hist(melted$Abundance)

# check test
melted %>%
  group_by(Sample) %>%
  summarise_at(c("Abundance"), sum)

# save results
groupData_23Sphoto_rar <- as.data.frame(melted)

# replace Nan by 0
groupData_23Sphoto_rar <- groupData_23Sphoto_rar %>% mutate_at(vars(Abundance), ~ replace(., is.nan(.), 0))

#create Class data set
class23 <- as.data.frame(groupData_23Sphoto_rar)
class23$Class <- as.factor(class23$OTU)

# Model:
mod <-lme(Abundance ~ Class*Location + Class*Depth, random=~1|Replicate, data = na.omit(class23), na.action=na.omit) 
anova(mod); #summary(mod)
hist(mod$residuals)
hist(as.numeric(class23$Abundance))

# per site - loop
list_class23 <- split(class23, f = class23$Class)
class23 <- na.omit(class23)

library(effectsize)
class <- 1:41
output <- lapply(X = class, FUN = function(x){
  mod <- lme(Abundance ~ Location * Depth, random=~1|Replicate, data = list_class23[[x]], na.action=na.omit) 
  fvalueL <- anova(mod)$`F-value`[2]
  fvalueD <- anova(mod)$`F-value`[3]
  fvalueLD <- anova(mod)$`F-value`[4]
  pvalueL <- anova(mod)$`p-value`[2]
  pvalueD <- anova(mod)$`p-value`[3]
  pvalueLD <- anova(mod)$`p-value`[4]
  eta2 <- omega_squared(mod)
  esL <- eta2$Omega2_partial[1]
  esD <- eta2$Omega2_partial[2]
  esLD <- eta2$Omega2_partial[3]
  eslowL <- eta2$CI_low[1]
  eslowD <- eta2$CI_low[2]
  eslowLD <- eta2$CI_low[3]
  eshighL <- eta2$CI_high[1]
  eshighD <- eta2$CI_high[2]
  eshighLD <- eta2$CI_high[3]
  ind.names <- c("fvalueL", "fvalueD", 'fvalueLD', "pvalueL", "pvalueD", "pvalueLD",
                 "esL", "esD", "esLD","eslowL", "eslowD", "eslowLD",
                 "eshighL", "eshighD", "eshighLD")
  ind.values <- c(fvalueL, fvalueD, fvalueLD, pvalueL, pvalueD, pvalueLD,
                  esL, esD, esLD,eslowL, eslowD, eslowLD,
                  eshighL, eshighD, eshighLD)
  es_mat <- as.data.frame(cbind(ind.names, ind.values))
  es_mat$name <- list_class23[[x]]$Class[1]
  return(es_mat)
})
output

es_class23 <- do.call(rbind, output) #to merge all the tables
es_class23 <- spread(es_class23, ind.names, ind.values) # reorganise

# plot:
plot_location <- ggplot(es_class23, aes(x = name,  y = as.numeric(pvalueL))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_location

plot_depth<- ggplot(es_class23, aes(x = name,  y = as.numeric(pvalueD))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_depth


plot_loc_depth <- ggplot(es_class23, aes(x = name,  y = as.numeric(pvalueLD))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_loc_depth

# |  | cbbL ----
EF_heatmap_cbbL <- aggregate(ASV_cbbL, by = list(as.character(metadata_cbbL$ID_cbbL)), FUN = sum)
rownames(EF_heatmap_cbbL) <- EF_heatmap_cbbL$Group.1
EF_heatmap_cbbL <- EF_heatmap_cbbL[,-1]

# Data prep:
# Aggregation by taxononic rank
EF_heatmap_cbbL_2 <- aggregate(as.data.frame(t(EF_heatmap_cbbL)), by = list(as.character(Taxo_cbbL$Class)), FUN = sum)

EF_heatmap_cbbL_2 <- subset(EF_heatmap_cbbL_2, 
                            EF_heatmap_cbbL_2$Group.1 != "Multi-affiliation")

# Creation of a phyloseq object :
EF_heatmap_cbbL_2_mat <- as.data.frame(EF_heatmap_cbbL_2 )
rownames(EF_heatmap_cbbL_2_mat) <- EF_heatmap_cbbL_2_mat$Group.1
EF_heatmap_cbbL_2_mat<- EF_heatmap_cbbL_2_mat[,-1]

EF_heatmap_cbbL_OTU = otu_table(EF_heatmap_cbbL_2_mat, taxa_are_rows = TRUE)
EF_heatmap_cbbL_samples = sample_data(metadata_cbbL)
rownames(EF_heatmap_cbbL_samples) <- EF_heatmap_cbbL_samples$ID_cbbL
Phyloseq_cbbL<- phyloseq(EF_heatmap_cbbL_OTU, EF_heatmap_cbbL_samples)
Phyloseq_cbbL

Phyloseq_cbbL_relab <- transform_sample_counts(Phyloseq_cbbL, function(x) x/sum(x))
Phyloseq_cbbL_relab@otu_table

# merge samples
melted_cbbL <- psmelt(Phyloseq_cbbL_relab)
hist(melted_cbbL$Abundance)

# checktest
melted_cbbL %>%
  group_by(Sample) %>%
  summarise_at(c("Abundance"), sum)

# save results
groupData_cbbL <- as.data.frame(melted_cbbL)

# replace Nan by 0
groupData_cbbL <- groupData_cbbL %>% mutate_at(vars(Abundance), ~ replace(., is.nan(.), 0))

# create Class data set
classcbbL <- as.data.frame(groupData_cbbL)
classcbbL$Class <- as.factor(classcbbL$OTU)

# Model:
mod <-lme(Abundance ~ Class*Location + Class*Depth, random=~1|Replicate, data = na.omit(classcbbL), na.action=na.omit) 
anova(mod); #summary(mod)
hist(mod$residuals)

# per site - loop
list_classcbbL <- split(classcbbL, f = classcbbL$Class)
classcbbL <- na.omit(classcbbL)

library(effectsize)
class <- 1:4
output_cbbL <- lapply(X = class, FUN = function(x){
  mod <- lme(Abundance ~ Location * Depth, random=~1|Replicate, data = list_classcbbL[[x]], na.action=na.omit) 
  fvalueL <- anova(mod)$`F-value`[2]
  fvalueD <- anova(mod)$`F-value`[3]
  fvalueLD <- anova(mod)$`F-value`[4]
  pvalueL <- anova(mod)$`p-value`[2]
  pvalueD <- anova(mod)$`p-value`[3]
  pvalueLD <- anova(mod)$`p-value`[4]
  eta2 <- omega_squared(mod)
  esL <- eta2$Omega2_partial[1]
  esD <- eta2$Omega2_partial[2]
  esLD <- eta2$Omega2_partial[3]
  eslowL <- eta2$CI_low[1]
  eslowD <- eta2$CI_low[2]
  eslowLD <- eta2$CI_low[3]
  eshighL <- eta2$CI_high[1]
  eshighD <- eta2$CI_high[2]
  eshighLD <- eta2$CI_high[3]
  ind.names <- c("fvalueL", "fvalueD", 'fvalueLD', "pvalueL", "pvalueD", "pvalueLD",
                 "esL", "esD", "esLD","eslowL", "eslowD", "eslowLD",
                 "eshighL", "eshighD", "eshighLD")
  ind.values <- c(fvalueL, fvalueD, fvalueLD, pvalueL, pvalueD, pvalueLD,
                  esL, esD, esLD,eslowL, eslowD, eslowLD,
                  eshighL, eshighD, eshighLD)
  es_mat <- as.data.frame(cbind(ind.names, ind.values))
  es_mat$name <- list_classcbbL[[x]]$Class[1]
  return(es_mat)
})
output_cbbL

es_classcbbL <- do.call(rbind, output_cbbL) #to merge all the tables
es_classcbbL <- spread(es_classcbbL, ind.names, ind.values) # reorganise

# plot:
plot_location <- ggplot(es_classcbbL, aes(x = name,  y = as.numeric(pvalueL))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_location

plot_depth<- ggplot(es_classcbbL, aes(x = name,  y = as.numeric(pvalueD))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_depth

plot_loc_depth <- ggplot(es_classcbbL, aes(x = name,  y = as.numeric(pvalueLD))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_loc_depth

# |  | bchY ----
EF_heatmap_bchY <- aggregate(ASV_bchY, by = list(as.character(metadata_bchY$ID_bchY)), FUN = sum)
rownames(EF_heatmap_bchY) <- EF_heatmap_bchY$Group.1
EF_heatmap_bchY <- EF_heatmap_bchY[,-1]

# Data preparation:
# Aggregation by taxononic rank
EF_heatmap_bchY_2 <- aggregate(as.data.frame(t(EF_heatmap_bchY)), by = list(as.character(Taxo_bchY$Class)), FUN = sum)

EF_heatmap_bchY_2 <- subset(EF_heatmap_bchY_2, 
                            EF_heatmap_bchY_2$Group.1 != "Multi-affiliation")

# Creation of a phyloseq object :
EF_heatmap_bchY_2_mat <- as.data.frame(EF_heatmap_bchY_2 )
rownames(EF_heatmap_bchY_2_mat) <- EF_heatmap_bchY_2_mat$Group.1
EF_heatmap_bchY_2_mat<- EF_heatmap_bchY_2_mat[,-1]

EF_heatmap_bchY_OTU = otu_table(EF_heatmap_bchY_2_mat, taxa_are_rows = TRUE)
EF_heatmap_bchY_samples = sample_data(metadata_bchY)
rownames(EF_heatmap_bchY_samples) <- EF_heatmap_bchY_samples$ID_bchY
Phyloseq_bchY<- phyloseq(EF_heatmap_bchY_OTU, EF_heatmap_bchY_samples)
Phyloseq_bchY

Phyloseq_bchY_relab <- transform_sample_counts(Phyloseq_bchY, function(x) x/sum(x))
Phyloseq_bchY_relab@otu_table

# merge samples
melted_bchY <- psmelt(Phyloseq_bchY_relab)
hist(melted_bchY$Abundance)

# checktest
melted_bchY %>%
  group_by(Sample) %>%
  summarise_at(c("Abundance"), sum)

# save results
groupData_bchY <- as.data.frame(melted_bchY)

# replace Nan by 0
groupData_bchY <- groupData_bchY %>% mutate_at(vars(Abundance), ~ replace(., is.nan(.), 0))

# create Class data set
classbchY <- as.data.frame(groupData_bchY)
classbchY$Class <- as.factor(classbchY$OTU)

# Model:
mod <-lme(Abundance ~ Class*Location + Class*Depth, random=~1|Replicate, data = na.omit(classbchY), na.action=na.omit) 
anova(mod); #summary(mod)
hist(mod$residuals)

# per site - loop
list_classbchY <- split(classbchY, f = classbchY$Class)
classbchY <- na.omit(classbchY)

library(effectsize)
class <- 1:22
output_bchY <- lapply(X = class, FUN = function(x){
  mod <- lme(Abundance ~ Location * Depth, random=~1|Replicate, data = list_classbchY[[x]], na.action=na.omit) 
  fvalueL <- anova(mod)$`F-value`[2]
  fvalueD <- anova(mod)$`F-value`[3]
  fvalueLD <- anova(mod)$`F-value`[4]
  pvalueL <- anova(mod)$`p-value`[2]
  pvalueD <- anova(mod)$`p-value`[3]
  pvalueLD <- anova(mod)$`p-value`[4]
  eta2 <- omega_squared(mod)
  esL <- eta2$Omega2_partial[1]
  esD <- eta2$Omega2_partial[2]
  esLD <- eta2$Omega2_partial[3]
  eslowL <- eta2$CI_low[1]
  eslowD <- eta2$CI_low[2]
  eslowLD <- eta2$CI_low[3]
  eshighL <- eta2$CI_high[1]
  eshighD <- eta2$CI_high[2]
  eshighLD <- eta2$CI_high[3]
  ind.names <- c("fvalueL", "fvalueD", 'fvalueLD', "pvalueL", "pvalueD", "pvalueLD",
                 "esL", "esD", "esLD","eslowL", "eslowD", "eslowLD",
                 "eshighL", "eshighD", "eshighLD")
  ind.values <- c(fvalueL, fvalueD, fvalueLD, pvalueL, pvalueD, pvalueLD,
                  esL, esD, esLD,eslowL, eslowD, eslowLD,
                  eshighL, eshighD, eshighLD)
  es_mat <- as.data.frame(cbind(ind.names, ind.values))
  es_mat$name <- list_classbchY[[x]]$Class[1]
  return(es_mat)
})
output_bchY

es_classbchY <- do.call(rbind, output_bchY) #to merge all the tables
es_classbchY <- spread(es_classbchY, ind.names, ind.values) # reorganise

# plot:
plot_location <- ggplot(es_classbchY, aes(x = name,  y = as.numeric(pvalueL))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_location

plot_depth <- ggplot(es_classbchY, aes(x = name,  y = as.numeric(pvalueD))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_depth


plot_loc_depth <- ggplot(es_classbchY, aes(x = name,  y = as.numeric(pvalueLD))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_loc_depth

# Add statistics for family ----
# | | 23S ----
EF_heatmap_23S <- aggregate(ASV_23S, by = list(as.character(metadata_23S$ID_23S)), FUN = sum)
rownames(EF_heatmap_23S) <- EF_heatmap_23S$Group.1
EF_heatmap_23S <- EF_heatmap_23S[,-1]

# Data preparation:
# Aggregation by taxononic rank : 
EF_heatmap_23S_2 <- aggregate(as.data.frame(t(EF_heatmap_23S)), by = list(as.character(Taxo_23S$Family)), FUN = sum)

# Relative abundance :
EF_heatmap_23S_2   <- subset(EF_heatmap_23S_2  , 
                             EF_heatmap_23S_2$Group.1 != "Multi-affiliation")


# Creation of a phyloseq object :
EF_heatmap_23S_2_mat <- as.data.frame(EF_heatmap_23S_2)
rownames(EF_heatmap_23S_2_mat) <- EF_heatmap_23S_2_mat$Group.1
EF_heatmap_23S_2_mat<- EF_heatmap_23S_2_mat[,-1]

EF_heatmap_23S_OTU = otu_table(EF_heatmap_23S_2_mat, taxa_are_rows = TRUE)
EF_heatmap_23S_samples = sample_data(metadata_23S)
rownames(EF_heatmap_23S_samples) <- EF_heatmap_23S_samples$ID_23S
Phyloseq_23S <- phyloseq(EF_heatmap_23S_OTU, EF_heatmap_23S_samples)
Phyloseq_23S

Phyloseq_23S_relab <- transform_sample_counts(Phyloseq_23S, function(x) x/sum(x))

# merge samples
melted <- psmelt(Phyloseq_23S_relab)
hist(melted$Abundance)

# check test
melted %>%
  group_by(Sample) %>%
  summarise_at(c("Abundance"), sum)

# save results
groupData_23Sphoto_rar <- as.data.frame(melted)

# replace Nan by 0
groupData_23Sphoto_rar <- groupData_23Sphoto_rar %>% mutate_at(vars(Abundance), ~ replace(., is.nan(.), 0))

#create Class data set
class23 <- as.data.frame(groupData_23Sphoto_rar)
class23$Class <- as.factor(class23$OTU)

# Model:
mod <-lme(Abundance ~ Class*Location + Class*Depth, random=~1|Replicate, data = na.omit(class23), na.action=na.omit) 
anova(mod); #summary(mod)
hist(mod$residuals)
hist(as.numeric(class23$Abundance))

# per site - loop
list_class23 <- split(class23, f = class23$Class)
class23 <- na.omit(class23)

library(effectsize)
class <- 1:81
output <- lapply(X = class, FUN = function(x){
  mod <- lme(Abundance ~ Location * Depth, random=~1|Replicate, data = list_class23[[x]], na.action=na.omit) 
  fvalueL <- anova(mod)$`F-value`[2]
  fvalueD <- anova(mod)$`F-value`[3]
  fvalueLD <- anova(mod)$`F-value`[4]
  pvalueL <- anova(mod)$`p-value`[2]
  pvalueD <- anova(mod)$`p-value`[3]
  pvalueLD <- anova(mod)$`p-value`[4]
  eta2 <- omega_squared(mod)
  esL <- eta2$Omega2_partial[1]
  esD <- eta2$Omega2_partial[2]
  esLD <- eta2$Omega2_partial[3]
  eslowL <- eta2$CI_low[1]
  eslowD <- eta2$CI_low[2]
  eslowLD <- eta2$CI_low[3]
  eshighL <- eta2$CI_high[1]
  eshighD <- eta2$CI_high[2]
  eshighLD <- eta2$CI_high[3]
  ind.names <- c("fvalueL", "fvalueD", 'fvalueLD', "pvalueL", "pvalueD", "pvalueLD",
                 "esL", "esD", "esLD","eslowL", "eslowD", "eslowLD",
                 "eshighL", "eshighD", "eshighLD")
  ind.values <- c(fvalueL, fvalueD, fvalueLD, pvalueL, pvalueD, pvalueLD,
                  esL, esD, esLD,eslowL, eslowD, eslowLD,
                  eshighL, eshighD, eshighLD)
  es_mat <- as.data.frame(cbind(ind.names, ind.values))
  es_mat$name <- list_class23[[x]]$Class[1]
  return(es_mat)
})
output

es_class23 <- do.call(rbind, output) #to merge all the tables
es_class23 <- spread(es_class23, ind.names, ind.values) # reorganise

# plot:
plot_location <- ggplot(es_class23, aes(x = name,  y = as.numeric(pvalueL))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_location

plot_depth<- ggplot(es_class23, aes(x = name,  y = as.numeric(pvalueD))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_depth


plot_loc_depth <- ggplot(es_class23, aes(x = name,  y = as.numeric(pvalueLD))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_loc_depth

# |  | cbbL ----
EF_heatmap_cbbL <- aggregate(ASV_cbbL, by = list(as.character(metadata_cbbL$ID_cbbL)), FUN = sum)
rownames(EF_heatmap_cbbL) <- EF_heatmap_cbbL$Group.1
EF_heatmap_cbbL <- EF_heatmap_cbbL[,-1]

# Data prep:
# Aggregation by taxononic rank
EF_heatmap_cbbL_2 <- aggregate(as.data.frame(t(EF_heatmap_cbbL)), by = list(as.character(Taxo_cbbL$Family)), FUN = sum)

EF_heatmap_cbbL_2 <- subset(EF_heatmap_cbbL_2, 
                            EF_heatmap_cbbL_2$Group.1 != "Multi-affiliation")

# Creation of a phyloseq object :
EF_heatmap_cbbL_2_mat <- as.data.frame(EF_heatmap_cbbL_2 )
rownames(EF_heatmap_cbbL_2_mat) <- EF_heatmap_cbbL_2_mat$Group.1
EF_heatmap_cbbL_2_mat<- EF_heatmap_cbbL_2_mat[,-1]

EF_heatmap_cbbL_OTU = otu_table(EF_heatmap_cbbL_2_mat, taxa_are_rows = TRUE)
EF_heatmap_cbbL_samples = sample_data(metadata_cbbL)
rownames(EF_heatmap_cbbL_samples) <- EF_heatmap_cbbL_samples$ID_cbbL
Phyloseq_cbbL<- phyloseq(EF_heatmap_cbbL_OTU, EF_heatmap_cbbL_samples)
Phyloseq_cbbL

Phyloseq_cbbL_relab <- transform_sample_counts(Phyloseq_cbbL, function(x) x/sum(x))
Phyloseq_cbbL_relab@otu_table

# merge samples
melted_cbbL <- psmelt(Phyloseq_cbbL_relab)
hist(melted_cbbL$Abundance)

# checktest
melted_cbbL %>%
  group_by(Sample) %>%
  summarise_at(c("Abundance"), sum)

# save results
groupData_cbbL <- as.data.frame(melted_cbbL)

# replace Nan by 0
groupData_cbbL <- groupData_cbbL %>% mutate_at(vars(Abundance), ~ replace(., is.nan(.), 0))

# create Class data set
classcbbL <- as.data.frame(groupData_cbbL)
classcbbL$Class <- as.factor(classcbbL$OTU)

# Model:
mod <-lme(Abundance ~ Class*Location + Class*Depth, random=~1|Replicate, data = na.omit(classcbbL), na.action=na.omit) 
anova(mod); #summary(mod)
hist(mod$residuals)

# per site - loop
list_classcbbL <- split(classcbbL, f = classcbbL$Class)
classcbbL <- na.omit(classcbbL)

library(effectsize)
class <- 1:15
output_cbbL <- lapply(X = class, FUN = function(x){
  mod <- lme(Abundance ~ Location * Depth, random=~1|Replicate, data = list_classcbbL[[x]], na.action=na.omit) 
  fvalueL <- anova(mod)$`F-value`[2]
  fvalueD <- anova(mod)$`F-value`[3]
  fvalueLD <- anova(mod)$`F-value`[4]
  pvalueL <- anova(mod)$`p-value`[2]
  pvalueD <- anova(mod)$`p-value`[3]
  pvalueLD <- anova(mod)$`p-value`[4]
  eta2 <- omega_squared(mod)
  esL <- eta2$Omega2_partial[1]
  esD <- eta2$Omega2_partial[2]
  esLD <- eta2$Omega2_partial[3]
  eslowL <- eta2$CI_low[1]
  eslowD <- eta2$CI_low[2]
  eslowLD <- eta2$CI_low[3]
  eshighL <- eta2$CI_high[1]
  eshighD <- eta2$CI_high[2]
  eshighLD <- eta2$CI_high[3]
  ind.names <- c("fvalueL", "fvalueD", 'fvalueLD', "pvalueL", "pvalueD", "pvalueLD",
                 "esL", "esD", "esLD","eslowL", "eslowD", "eslowLD",
                 "eshighL", "eshighD", "eshighLD")
  ind.values <- c(fvalueL, fvalueD, fvalueLD, pvalueL, pvalueD, pvalueLD,
                  esL, esD, esLD,eslowL, eslowD, eslowLD,
                  eshighL, eshighD, eshighLD)
  es_mat <- as.data.frame(cbind(ind.names, ind.values))
  es_mat$name <- list_classcbbL[[x]]$Class[1]
  return(es_mat)
})
output_cbbL

es_classcbbL <- do.call(rbind, output_cbbL) #to merge all the tables
es_classcbbL <- spread(es_classcbbL, ind.names, ind.values) # reorganise

# plot:
plot_location <- ggplot(es_classcbbL, aes(x = name,  y = as.numeric(pvalueL))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_location

plot_depth<- ggplot(es_classcbbL, aes(x = name,  y = as.numeric(pvalueD))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_depth

plot_loc_depth <- ggplot(es_classcbbL, aes(x = name,  y = as.numeric(pvalueLD))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_loc_depth

# |  | bchY ----
EF_heatmap_bchY <- aggregate(ASV_bchY, by = list(as.character(metadata_bchY$ID_bchY)), FUN = sum)
rownames(EF_heatmap_bchY) <- EF_heatmap_bchY$Group.1
EF_heatmap_bchY <- EF_heatmap_bchY[,-1]

# Data preparation:
# Aggregation by taxononic rank
EF_heatmap_bchY_2 <- aggregate(as.data.frame(t(EF_heatmap_bchY)), by = list(as.character(Taxo_bchY$Family)), FUN = sum)

EF_heatmap_bchY_2 <- subset(EF_heatmap_bchY_2, 
                            EF_heatmap_bchY_2$Group.1 != "Multi-affiliation")

# Creation of a phyloseq object :
EF_heatmap_bchY_2_mat <- as.data.frame(EF_heatmap_bchY_2 )
rownames(EF_heatmap_bchY_2_mat) <- EF_heatmap_bchY_2_mat$Group.1
EF_heatmap_bchY_2_mat<- EF_heatmap_bchY_2_mat[,-1]

EF_heatmap_bchY_OTU = otu_table(EF_heatmap_bchY_2_mat, taxa_are_rows = TRUE)
EF_heatmap_bchY_samples = sample_data(metadata_bchY)
rownames(EF_heatmap_bchY_samples) <- EF_heatmap_bchY_samples$ID_bchY
Phyloseq_bchY<- phyloseq(EF_heatmap_bchY_OTU, EF_heatmap_bchY_samples)
Phyloseq_bchY

Phyloseq_bchY_relab <- transform_sample_counts(Phyloseq_bchY, function(x) x/sum(x))
Phyloseq_bchY_relab@otu_table

# merge samples
melted_bchY <- psmelt(Phyloseq_bchY_relab)
hist(melted_bchY$Abundance)

# checktest
melted_bchY %>%
  group_by(Sample) %>%
  summarise_at(c("Abundance"), sum)

# save results
groupData_bchY <- as.data.frame(melted_bchY)

# replace Nan by 0
groupData_bchY <- groupData_bchY %>% mutate_at(vars(Abundance), ~ replace(., is.nan(.), 0))

# create Class data set
classbchY <- as.data.frame(groupData_bchY)
classbchY$Class <- as.factor(classbchY$OTU)

# Model:
mod <-lme(Abundance ~ Class*Location + Class*Depth, random=~1|Replicate, data = na.omit(classbchY), na.action=na.omit) 
anova(mod); #summary(mod)
hist(mod$residuals)

# per site - loop
list_classbchY <- split(classbchY, f = classbchY$Class)
classbchY <- na.omit(classbchY)

library(effectsize)
class <- 1:37
output_bchY <- lapply(X = class, FUN = function(x){
  mod <- lme(Abundance ~ Location * Depth, random=~1|Replicate, data = list_classbchY[[x]], na.action=na.omit) 
  fvalueL <- anova(mod)$`F-value`[2]
  fvalueD <- anova(mod)$`F-value`[3]
  fvalueLD <- anova(mod)$`F-value`[4]
  pvalueL <- anova(mod)$`p-value`[2]
  pvalueD <- anova(mod)$`p-value`[3]
  pvalueLD <- anova(mod)$`p-value`[4]
  eta2 <- omega_squared(mod)
  esL <- eta2$Omega2_partial[1]
  esD <- eta2$Omega2_partial[2]
  esLD <- eta2$Omega2_partial[3]
  eslowL <- eta2$CI_low[1]
  eslowD <- eta2$CI_low[2]
  eslowLD <- eta2$CI_low[3]
  eshighL <- eta2$CI_high[1]
  eshighD <- eta2$CI_high[2]
  eshighLD <- eta2$CI_high[3]
  ind.names <- c("fvalueL", "fvalueD", 'fvalueLD', "pvalueL", "pvalueD", "pvalueLD",
                 "esL", "esD", "esLD","eslowL", "eslowD", "eslowLD",
                 "eshighL", "eshighD", "eshighLD")
  ind.values <- c(fvalueL, fvalueD, fvalueLD, pvalueL, pvalueD, pvalueLD,
                  esL, esD, esLD,eslowL, eslowD, eslowLD,
                  eshighL, eshighD, eshighLD)
  es_mat <- as.data.frame(cbind(ind.names, ind.values))
  es_mat$name <- list_classbchY[[x]]$Class[1]
  return(es_mat)
})
output_bchY

es_classbchY <- do.call(rbind, output_bchY) #to merge all the tables
es_classbchY <- spread(es_classbchY, ind.names, ind.values) # reorganise

# plot:
plot_location <- ggplot(es_classbchY, aes(x = name,  y = as.numeric(pvalueL))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_location

plot_depth <- ggplot(es_classbchY, aes(x = name,  y = as.numeric(pvalueD))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_depth


plot_loc_depth <- ggplot(es_classbchY, aes(x = name,  y = as.numeric(pvalueLD))) +
  geom_point(aes(color = name), position = position_dodge(0.3), size = 3)+
  geom_hline(yintercept = 0.05, lty = 2, col = "red") + 
  geom_hline(yintercept = 0.01, lty = 2, col = "blue") + 
  geom_hline(yintercept = 0.001, lty = 2, col = "green") + 
  theme_classic() +
  coord_flip()
plot_loc_depth

# 6. Co-occurence graphics ----
# | | 23S ----
asv_23S_loc_depth <- aggregate(ASV_23S, by = list(as.character(metadata_23S$Location_Depth)), FUN = sum)
rownames(asv_23S_loc_depth) <- asv_23S_loc_depth$Group.1
asv_23S_loc_depth <- asv_23S_loc_depth[,-1]
asv_23S_loc_depth2 <- aggregate(as.data.frame(t(asv_23S_loc_depth)), by = list(as.character(Taxo_23S$Class)), FUN = sum)
asv_23S_loc_depth2 <- subset(asv_23S_loc_depth2, 
                             asv_23S_loc_depth2$Group.1 != "Multi-affiliation")

test_23S <- asv_23S_loc_depth2

test_23S$Abisko_1[test_23S$Abisko_1 > 0] <- 1
test_23S$Abisko_2[test_23S$Abisko_2 > 0] <- 1
test_23S$Abisko_3[test_23S$Abisko_3 > 0] <- 1

test_23S$Counozouls_1[test_23S$Counozouls_1 > 0] <- 1
test_23S$Counozouls_2[test_23S$Counozouls_2 > 0] <- 1
test_23S$Counozouls_3[test_23S$Counozouls_3 > 0] <- 1

test_23S$Mann_1[test_23S$Mann_1 > 0] <- 1
test_23S$Mann_2[test_23S$Mann_2 > 0] <- 1
test_23S$Mann_3[test_23S$Mann_3 > 0] <- 1

test_23S$Siikaneva_1[test_23S$Siikaneva_1 > 0] <- 1
test_23S$Siikaneva_2[test_23S$Siikaneva_2 > 0] <- 1
test_23S$Siikaneva_3[test_23S$Siikaneva_3 > 0] <- 1

upset(test_23S,nsets = 15,
      order.by = c("freq", "degree"), decreasing = c(FALSE,TRUE))


# | | cbbL ----
asv_cbbL_loc_depth <- aggregate(ASV_cbbL, by = list(as.character(metadata_cbbL$Location_Depth)), FUN = sum)
rownames(asv_cbbL_loc_depth) <- asv_cbbL_loc_depth$Group.1
asv_cbbL_loc_depth <- asv_cbbL_loc_depth[,-1]
asv_cbbL_loc_depth2 <- aggregate(as.data.frame(t(asv_cbbL_loc_depth)), by = list(as.character(Taxo_cbbL$Class)), FUN = sum)
asv_cbbL_loc_depth2 <- subset(asv_cbbL_loc_depth2, 
                              asv_cbbL_loc_depth2$Group.1 != "Multi-affiliation")

test_cbbL <- asv_cbbL_loc_depth2
test_cbbL$Abisko_1[test_cbbL$Abisko_1 > 0] <- 1
test_cbbL$Abisko_2[test_cbbL$Abisko_2 > 0] <- 1
test_cbbL$Abisko_3[test_cbbL$Abisko_3 > 0] <- 1

test_cbbL$Counozouls_1[test_cbbL$Counozouls_1 > 0] <- 1
test_cbbL$Counozouls_2[test_cbbL$Counozouls_2 > 0] <- 1
test_cbbL$Counozouls_3[test_cbbL$Counozouls_3 > 0] <- 1

test_cbbL$Mann_1[test_cbbL$Mann_1 > 0] <- 1
test_cbbL$Mann_2[test_cbbL$Mann_2 > 0] <- 1
test_cbbL$Mann_3[test_cbbL$Mann_3 > 0] <- 1

test_cbbL$Siikaneva_1[test_cbbL$Siikaneva_1 > 0] <- 1
test_cbbL$Siikaneva_2[test_cbbL$Siikaneva_2 > 0] <- 1
test_cbbL$Siikaneva_3[test_cbbL$Siikaneva_3 > 0] <- 1
test_cbbL$Group.1 <- as.factor(test_cbbL$Group.1)

upset(test_cbbL,nsets = 15,
      order.by = c("freq", "degree"), decreasing = c(FALSE,TRUE))


# | | bchY ----
asv_bchY_loc_depth <- aggregate(ASV_bchY, by = list(as.character(metadata_bchY$Location_Depth)), FUN = sum)
rownames(asv_bchY_loc_depth) <- asv_bchY_loc_depth$Group.1
asv_bchY_loc_depth <- asv_bchY_loc_depth[,-1]
asv_bchY_loc_depth2 <- aggregate(as.data.frame(t(asv_bchY_loc_depth)), by = list(as.character(Taxo_bchY$Class)), FUN = sum)
asv_bchY_loc_depth2 <- subset(asv_bchY_loc_depth2, 
                              asv_bchY_loc_depth2$Group.1 != "Multi-affiliation")

test_bchY <- asv_bchY_loc_depth2
test_bchY$Abisko_1[test_bchY$Abisko_1 > 0] <- 1
test_bchY$Abisko_2[test_bchY$Abisko_2 > 0] <- 1
test_bchY$Abisko_3[test_bchY$Abisko_3 > 0] <- 1

test_bchY$Counozouls_1[test_bchY$Counozouls_1 > 0] <- 1
test_bchY$Counozouls_2[test_bchY$Counozouls_2 > 0] <- 1
test_bchY$Counozouls_3[test_bchY$Counozouls_3 > 0] <- 1

test_bchY$Mann_1[test_bchY$Mann_1 > 0] <- 1
test_bchY$Mann_2[test_bchY$Mann_2 > 0] <- 1
test_bchY$Mann_3[test_bchY$Mann_3 > 0] <- 1

test_bchY$Siikaneva_1[test_bchY$Siikaneva_1 > 0] <- 1
test_bchY$Siikaneva_2[test_bchY$Siikaneva_2 > 0] <- 1
test_bchY$Siikaneva_3[test_bchY$Siikaneva_3 > 0] <- 1
test_bchY$Group.1 <- as.factor(test_bchY$Group.1)

upset(test_bchY,nsets = 15,
      order.by = c("freq", "degree"), decreasing = c(FALSE,TRUE))
