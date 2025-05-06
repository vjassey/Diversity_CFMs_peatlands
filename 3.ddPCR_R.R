# Diversity, abundance, and biogeography of CO2-fixing microorganisms in peatlands
# Marie Le Geay, Kyle Mayers, Anna Sytiuk, Ellen Dorrepaal, Martin Küttim, 
# Mariusz Lamentowicz, Eeva-Stiina Tuittila, Béatrice Lauga, and Vincent E.J. Jassey

# Author: Marie Le Geay and Vincent Jassey

################################################################################
# 3. Absolute quantification of microbial genes 

# Main text:
## Figure 2

# Supplementary Information:
## Tables S8 and S9
## Figures S7 and S8

################################################################################
# |   |  Load packages ----
library(ggplot2); library(plyr); library(reshape2); library(grid); library(gridExtra); library(ape); library(picante)
library(seqinr); library(phytools); library(scales); library(nlme); library(vegan); library(RColorBrewer)
library("BiocManager"); library("Biobase"); library("rhdf5"); library("devtools"); library(phyloseq); library(msa)
library(ggtree); library(tidyjson); library(RColorBrewer);library(tidyverse); library(viridis); library(randomForest); library(caret); library(cluster)
library(tidyverse); library(cluster);  library(parallelDist); library(boral);  library(corrplot); library(igraph); library(car); library(fossil); library(ROCR); library(doMC)
library(effectsize); library(phyloseq); library(Biostrings); library(dplyr); library(tidyr); library(bioseq)
library(vegan); library(ade4); library(BAT); library(patchwork); library(hrbrthemes); library(gcookbook)
library(htmltools); library(HH); library(ggpubr); library(nlme)
library(nlme); library(emmeans); library(multcomp)


# |   |  Load data ----
getwd()
setwd("C:/Users/mlegeay/OneDrive/Documents/chapter_2/Redaction/Révision/RStudio-to-publish/Data")# set up a new path
data_ddPCR <- read.table("ddPCR.csv", header = T, sep = ";")
data_ddPCR <- as.data.frame(data_ddPCR)
data_ddPCR$Location <- fct_relevel(data_ddPCR$Location, 
                          c("COU","MAN",
                            "SIIK","ABI"))
data_ddPCR$Depth <- fct_relevel(data_ddPCR$Depth, c("Surface","Intermediate","Deep"))


metadata_16S <- read.table("meta_16S.csv", header = T, sep = ";")
metadata_16S$Location <- fct_relevel(metadata_16S$Location, 
                                     c("Counozouls","Mann",
                                       "Siikaneva","Abisko"))
metadata_23S <- read.table("meta_23S.csv", header = T, sep = ";")
metadata_23S$Location <- fct_relevel(metadata_23S$Location, 
                                     c("Counozouls","Mann",
                                       "Siikaneva","Abisko"))
metadata_cbbL <- read.table("meta_cbbL.csv", header = T, sep = ";")
metadata_cbbL$Location <- fct_relevel(metadata_cbbL$Location, 
                                     c("Counozouls","Mann",
                                       "Siikaneva","Abisko"))
metadata_bchY <- read.table("meta_bchY.csv", header = T, sep = ";")
metadata_bchY$Location <- fct_relevel(metadata_bchY$Location, 
                                     c("Counozouls","Mann",
                                       "Siikaneva","Abisko"))


# subset the cyanobacteria (gene targeting cynaobacteria)
cyano_ddPCR <- subset(data_ddPCR, Gene == "16S-cyano")
ddPCR_1 <- subset(data_ddPCR, Gene != "16S-cyano")
data_ddPCR_tot_abund <-subset(ddPCR_1, Gene != "16S")
ddPCR <- subset(ddPCR_1, Gene != "total")


# transform in log
ddPCR$log_copies <- log10(ddPCR$Copies.g.drypeat) 
cyano_ddPCR$log_copies <- log(cyano_ddPCR$Copies.g.drypeat)
data_ddPCR_tot_abund$log_copies <- log10(data_ddPCR_tot_abund$Copies.g.drypeat)


# 1. Visualisation for each gene ----
dodge2 <- position_dodge(width = 0.8)
g_16S <- ggplot(data = metadata_16S, aes(x = as.factor(Depth), y = log(as.numeric(AQ_16S)), fill = as.factor(Depth)))+
  geom_violin(position = dodge2, alpha = 0.2) +
  geom_boxplot(width =0.4,position = dodge2)+
  scale_fill_manual(values = c("#FED789FF","#3E938BFF", "#72874EFF","#70684DFF"))+
  labs(title="16S")+
  theme_classic()+
  facet_wrap(~Gene, nrow = 1)+
  facet_wrap(~Location, nrow = 1)
g_16S

g_23S <- ggplot(data = metadata_23S, aes(x = as.factor(Depth), y = log(as.numeric(AQ_23S)), fill = as.factor(Depth)))+
  geom_violin(position = dodge2, alpha = 0.2) +
  geom_boxplot(width =0.4,position = dodge2)+
  scale_fill_manual(values = c("#FED789FF","#3E938BFF", "#72874EFF","#70684DFF"))+
  labs(title="23S")+
  theme_classic()+
  facet_wrap(~Gene, nrow = 1)+
  facet_wrap(~Location, nrow = 1)
g_23S

g_cbbL <- ggplot(data = metadata_cbbL, aes(x = as.factor(Depth), y = log(as.numeric(AQ_cbbL)), fill = as.factor(Depth)))+
  geom_violin(position = dodge2, alpha = 0.2) +
  geom_boxplot(width =0.4,position = dodge2)+
  scale_fill_manual(values = c("#FED789FF","#3E938BFF", "#72874EFF","#70684DFF"))+
  labs(title="cbbL")+
  theme_classic()+
  facet_wrap(~Gene, nrow = 1)+
  facet_wrap(~Location, nrow = 1)
g_cbbL

g_pufM <- ggplot(data = metadata_bchY, aes(x = as.factor(Depth), y = log(as.numeric(AQ_pufM)), fill = as.factor(Depth)))+
  geom_violin(position = dodge2, alpha = 0.2) +
  geom_boxplot(width =0.4,position = dodge2)+
  scale_fill_manual(values = c("#FED789FF","#3E938BFF", "#72874EFF","#70684DFF"))+
  labs(title="pufM")+
  theme_classic()+
  facet_wrap(~Gene, nrow = 1)+
  facet_wrap(~Location, nrow = 1)
g_pufM

##  All graph together ----
p_all_AQ <- (g_16S | g_23S) / (g_cbbL | g_pufM)
p_all_AQ
p_all_AQ  + plot_annotation(
  title = 'ddPCR',
  tag_levels = 'a'
)

## Cyanobacteria ----
g_cyano <- ggplot(data = cyano_ddPCR, aes(x = as.factor(Depth), y = log_copies, fill = as.factor(Depth)))+
  geom_violin(position = dodge2, alpha = 0.2) +
  geom_boxplot(width =0.4,position = dodge2)+
  scale_fill_manual(values = c("#FED789FF","#3E938BFF", "#72874EFF","#70684DFF"))+
  labs(title="pufM")+
  theme_classic()+
  facet_wrap(~Gene, nrow = 1)+
  facet_wrap(~Location, nrow = 1)
g_cyano


# 2. Visualization per site and with total gene abundance ----
data_ddPCR_tot_abund$Gene <- fct_relevel(data_ddPCR_tot_abund$Gene, 
                                         c("total","23S",
                                           "cbbL","pufM"))
data_ddPCR_tot_abund$Location <- fct_relevel(data_ddPCR_tot_abund$Location , 
                                             c("COU","MAN",
                                               "SIIK","ABI"))
data_ddPCR_tot_abund$Depth <- fct_relevel(data_ddPCR_tot_abund$Depth , 
                                          c("Surface","Intermediate",
                                            "Deep"))

my_comparisons <- list( c("23S", "cbbL"), c("cbbL", "pufM"), c("23S", "pufM") )

g_tot_abund <- ggplot(data = data_ddPCR_tot_abund, aes(x = as.factor(Gene), y = log(as.numeric(Copies.g.drypeat)), fill = as.factor(Gene)))+
  geom_violin(position = dodge2, alpha = 0.2) +
  geom_boxplot(width =0.4,position = dodge2)+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  scale_fill_manual(values = c("#FED789FF","#3E938BFF", "#72874EFF","#70684DFF"))+
  labs(title="total")+
  theme_classic() +
  facet_wrap(~Location, nrow = 1)
g_tot_abund

# 3. Add statistical analyses ----
# Linear mixed models with lme function of nlme package
# Replicate by location as random term
ddPCR$Depth <- as.factor(ddPCR$Depth)
ddPCR_16S <- subset(ddPCR, ddPCR$Gene == "16S")
ddPCR_16S_Cou <- subset(ddPCR_16S, ddPCR_16S$Location == "COU")
ddPCR_16S_Man <- subset(ddPCR_16S, ddPCR_16S$Location == "MAN")
ddPCR_16S_Siik <- subset(ddPCR_16S, ddPCR_16S$Location == "SIIK")
ddPCR_16S_Abi <- subset(ddPCR_16S, ddPCR_16S$Location == "ABI")


ddPCR_23S <- subset(ddPCR, ddPCR$Gene == "23S")
ddPCR_23S_Cou <- subset(ddPCR_23S, ddPCR_23S$Location == "COU")
ddPCR_23S_Man <- subset(ddPCR_23S, ddPCR_23S$Location == "MAN")
ddPCR_23S_Siik <- subset(ddPCR_23S, ddPCR_23S$Location == "SIIK")
ddPCR_23S_Abi <- subset(ddPCR_23S, ddPCR_23S$Location == "ABI")

ddPCR_cbbL <- subset(ddPCR, ddPCR$Gene == "cbbL")
ddPCR_cbbL_Cou <- subset(ddPCR_cbbL, ddPCR_cbbL$Location == "COU")
ddPCR_cbbL_Man <- subset(ddPCR_cbbL, ddPCR_cbbL$Location == "MAN")
ddPCR_cbbL_Siik <- subset(ddPCR_cbbL, ddPCR_cbbL$Location == "SIIK")
ddPCR_cbbL_Abi <- subset(ddPCR_cbbL, ddPCR_cbbL$Location == "ABI")

ddPCR_pufM <- subset(ddPCR, ddPCR$Gene == "pufM")
ddPCR_pufM_Cou <- subset(ddPCR_pufM, ddPCR_pufM$Location == "COU")
ddPCR_pufM_Man <- subset(ddPCR_pufM, ddPCR_pufM$Location == "MAN")
ddPCR_pufM_Siik <- subset(ddPCR_pufM, ddPCR_pufM$Location == "SIIK")
ddPCR_pufM_Abi <- subset(ddPCR_pufM, ddPCR_pufM$Location == "ABI")

cyano_ddPCR_Cou <- subset(cyano_ddPCR, cyano_ddPCR$Location == "COU")
cyano_ddPCR_Man <- subset(cyano_ddPCR, cyano_ddPCR$Location == "MAN")
cyano_ddPCR_Siik <- subset(cyano_ddPCR, cyano_ddPCR$Location == "SIIK")
cyano_ddPCR_Abi <- subset(cyano_ddPCR, cyano_ddPCR$Location == "ABI")

# | | Effect of location ----
# Compute the model:
mod_loc <-lme(log_copies ~ Location, random=~1|Replicate/Location, 
              data = na.omit(ddPCR_16S), na.action=na.omit) 
# change 16S by other genes (23S, cbbL and bchY)
# Look at model results:
anova(mod_loc) 
# Check residuals:
hist(mod_loc$residuals)
# Post hoc analysis:
emmeans(mod_loc, pairwise ~ Location)

# | | Effect of depth ----
# Compute the model for each location:
mod_depth <-lme(log_copies ~ Depth, random=~1|Replicate/Location, 
                data = na.omit(ddPCR_16S_Cou), na.action=na.omit)
# change 23S_Cou by other gene x location (cbbL_Cou, cbbL_Abi...)
# Look at model results:
anova(mod_depth)
# Check residuals:
hist(mod_depth$residuals)
# Post hoc analysis:
emmeans(mod_depth, pairwise ~ Depth)


