# Diversity, abundance, and biogeography of CO2-fixing microorganisms in peatlands
# Marie Le Geay, Kyle Mayers, Anna Sytiuk, Ellen Dorrepaal, Martin Küttim, 
# Mariusz Lamentowicz, Eeva-Stiina Tuittila, Béatrice Lauga, and Vincent E.J. Jassey

# Author: Marie Le Geay and Vincent Jassey

################################################################################
# 2. Diversity analysis 

# Main text:
## Figure 2

# Supplementary Information:
## Tables S10 and S11
## Figures S9 and S10

################################################################################
# |   |  Load packages ----

library(ggplot2); library(plyr); library(reshape2); library(grid); library(gridExtra); library(ape); library(picante)
library(seqinr); library(phytools); library(scales); library(nlme); library(vegan); library(RColorBrewer)
library("BiocManager"); library("Biobase"); library("rhdf5"); library("devtools"); library(phyloseq); library(msa)
library(ggtree); library(tidyjson); library(RColorBrewer);library(tidyverse); library(viridis); library(randomForest); library(caret); library(cluster)
library(tidyverse); library(cluster);  library(parallelDist); library(boral);  library(corrplot); library(igraph); library(car); library(fossil); library(ROCR); library(doMC)
library(effectsize); library(phyloseq); library(Biostrings); library(dplyr); library(tidyr); library(bioseq)
library(vegan); library(ade4); library(BAT); library(patchwork); library(hrbrthemes); library(gcookbook)
library(htmltools); library(HH); library(ggpubr);library(nlme)
library(emmeans)
library(multcomp)

dodge2 <- position_dodge(width = 0.8)

# |   |  Load data ----
setwd("C:/Users/mlegeay/OneDrive/Documents/chapter_2/Redaction/Révision/RStudio-to-publish/Data")# set up a new path

metadata_23S <- read.table("meta_23S.csv", header = T, sep = ";")
metadata_cbbL <- read.table("meta_cbbL.csv", header = T, sep = ";")
metadata_bchY <- read.table("meta_bchY.csv", header = T, sep = ";")


ASV_23S <- read.table("t_ASV_23S.csv", header = T, sep = ",")
rownames(ASV_23S) <- ASV_23S$X
ASV_23S <- ASV_23S[,-1]

ASV_cbbL <- read.table("t_ASV_cbbL.csv", header = T, sep = ",")
rownames(ASV_cbbL) <- ASV_cbbL$X
ASV_cbbL <- ASV_cbbL[,-1]

ASV_bchY <- read.table("t_ASV_bchY.csv", header = T, sep = ",")
rownames(ASV_bchY) <- ASV_bchY$X
ASV_bchY <- ASV_bchY[,-1]


# 1. Alpha diversity ----
#### |  |  23S ----
# richness :
metadata_23S$S <- specnumber(ASV_23S) 

# composite :
metadata_23S$shannon <- diversity(ASV_23S,index="shannon")

# Visualisation :
metadata_23S$Location <- fct_relevel(metadata_23S$Location, 
                                 c("Counozouls","Mann",
                                   "Siikaneva","Abisko"))

g_23S_S <- ggplot(data = metadata_23S, aes(x = as.factor(Depth), y = S, fill = as.factor(Depth)))+
  geom_violin(position = dodge2, alpha = 0.2) +
  geom_boxplot(width =0.4,position = dodge2)+
  scale_fill_manual(values = c("#FED789FF","#3E938BFF", "#72874EFF"))+
  ylim(50,1000) +
  labs(title="23S")+
  theme_classic()+
  facet_wrap(~Location, nrow = 1)
g_23S_S

g_23S_shannon <- ggplot(data = metadata_23S, aes(x = as.factor(Depth), y = shannon, fill = as.factor(Depth)))+
  geom_violin(position = dodge2, alpha = 0.2) +
  geom_boxplot(width =0.4,position = dodge2)+
  scale_fill_manual(values = c("#FED789FF","#3E938BFF", "#72874EFF"))+
  labs(title="23S")+
  ylim(0,7) +
  theme_classic()+
  facet_wrap(~Location, nrow = 1)
g_23S_shannon


#### |  |  cbbL ----
# richness :
metadata_cbbL$S <- specnumber(ASV_cbbL) 

# composite :
metadata_cbbL$shannon <- diversity(ASV_cbbL,index="shannon")

# Visualisation :
metadata_cbbL$Location <- fct_relevel(metadata_cbbL$Location, 
                                     c("Counozouls","Mann",
                                       "Siikaneva","Abisko"))

g_cbbL_S <- ggplot(data = metadata_cbbL, aes(x = as.factor(Depth), y = S, fill = as.factor(Depth)))+
  geom_violin(position = dodge2, alpha = 0.2) +
  geom_boxplot(width =0.4,position = dodge2)+
  scale_fill_manual(values = c("#FED789FF","#3E938BFF", "#72874EFF"))+
  ylim(50,1000) +
  labs(title="cbbL")+
  theme_classic()+
  facet_wrap(~Location, nrow = 1)
g_cbbL_S 

g_cbbL_shannon <- ggplot(data = metadata_cbbL, aes(x = as.factor(Depth), y = shannon, fill = as.factor(Depth)))+
  geom_violin(position = dodge2, alpha = 0.2) +
  geom_boxplot(width =0.4,position = dodge2)+
  scale_fill_manual(values = c("#FED789FF","#3E938BFF", "#72874EFF"))+
  ylim(0,7) +
  labs(title="cbbL")+
  theme_classic()+
  facet_wrap(~Location, nrow = 1)
g_cbbL_shannon


#### |  |  bchY ----
# richness :
metadata_bchY$S <- specnumber(ASV_bchY) 

# composite :
metadata_bchY$shannon <- diversity(ASV_bchY,index="shannon")

# Visualisation :
metadata_bchY$Location <- fct_relevel(metadata_bchY$Location, 
                                      c("Counozouls","Mann",
                                        "Siikaneva","Abisko"))

g_bchY_S <- ggplot(data = metadata_bchY, aes(x = as.factor(Depth), y = S, fill = as.factor(Depth)))+
  geom_violin(position = dodge2, alpha = 0.2) +
  geom_boxplot(width =0.4,position = dodge2)+
  scale_fill_manual(values = c("#FED789FF","#3E938BFF", "#72874EFF"))+
  labs(title="bchY")+
  ylim(0,1000) +
  theme_classic()+
  facet_wrap(~Location, nrow = 1)
g_bchY_S

g_bchY_shannon <- ggplot(data = metadata_bchY, aes(x = as.factor(Depth), y = shannon, fill = as.factor(Depth)))+
  geom_violin(position = dodge2, alpha = 0.2) +
  geom_boxplot(width =0.4,position = dodge2)+
  scale_fill_manual(values = c("#FED789FF","#3E938BFF", "#72874EFF"))+
  labs(title="bchY")+
  ylim(0,7) +
  theme_classic()+
  facet_wrap(~Location, nrow = 1)
g_bchY_shannon


#### All graph together ----
p_all <- (g_23S_S | g_cbbL_S | g_bchY_S) / (g_23S_shannon  | g_cbbL_shannon | g_bchY_shannon) 
p_all
p_all + plot_annotation(
  title = 'Alpha div',
  tag_levels = 'a'
)


# 2. Richness total ----
data_rich_tot <- read.table("rich_all.csv", header = T, sep = ";")
data_rich_tot$Gene <- fct_relevel(data_rich_tot$Gene, 
                                  c("total","23S",
                                    "cbbL","bchY"))
data_rich_tot$Location <- fct_relevel(data_rich_tot$Location , 
                                      c("Counozouls","Mann",
                                        "Siikaneva","Abisko"))
data_rich_tot$Depth <- fct_relevel(data_rich_tot$Depth , 
                                   c("D1","D2",
                                     "D3"))

my_comparisons <- list( c("23S", "cbbL"), c("cbbL", "bchY"), c("23S", "bchY") )

g_rich <- ggplot(data = data_rich_tot, aes(x = as.factor(Gene), y = S, fill = as.factor(Gene)))+
  geom_violin(position = dodge2, alpha = 0.2) +
  geom_boxplot(width =0.4,position = dodge2)+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  scale_fill_manual(values = c("#FED789FF","#3E938BFF", "#72874EFF","#70684DFF"))+
  labs(title="total")+
  theme_classic() +
  facet_wrap(~Location, nrow = 1)
g_rich


# 3. Add statistical analyses ----
# Linear mixed models with lme function of nlme package
# Replicate by location as random term

# | | Effect of location ----
# Compute the model:
mod_loc <-lme(S ~ Location, random=~1|Replicate/Location, 
              data = na.omit(metadata_23S), na.action=na.omit) 
# change 23S by other genes (cbbL and bchY)
# change S by shannon
# Look at model results:
anova(mod_loc) 
# Check residuals:
hist(mod_loc$residuals)
# Post hoc analysis:
emmeans(mod_loc, pairwise ~ Location)

# | | Effect of depth ----
metadata_23S$Depth <- as.factor(metadata_23S$Depth)
metadata_cbbL$Depth <- as.factor(metadata_cbbL$Depth)
metadata_bchY$Depth <- as.factor(metadata_bchY$Depth)

metadata_23S_Cou <- subset(metadata_23S, metadata_23S$Location == "Counozouls")
metadata_23S_Man <- subset(metadata_23S, metadata_23S$Location == "Mann")
metadata_23S_Siik <- subset(metadata_23S, metadata_23S$Location == "Siikaneva")
metadata_23S_Abi <- subset(metadata_23S, metadata_23S$Location == "Abisko")

metadata_cbbL_Cou <- subset(metadata_cbbL, metadata_cbbL$Location == "Counozouls")
metadata_cbbL_Man <- subset(metadata_cbbL, metadata_cbbL$Location == "Mann")
metadata_cbbL_Siik <- subset(metadata_cbbL, metadata_cbbL$Location == "Siikaneva")
metadata_cbbL_Abi <- subset(metadata_cbbL, metadata_cbbL$Location == "Abisko")

metadata_bchY_Cou <- subset(metadata_bchY, metadata_bchY$Location == "Counozouls")
metadata_bchY_Man <- subset(metadata_bchY, metadata_bchY$Location == "Mann")
metadata_bchY_Siik <- subset(metadata_bchY, metadata_bchY$Location == "Siikaneva")
metadata_bchY_Abi <- subset(metadata_bchY, metadata_bchY$Location == "Abisko")

# Compute the model for each location:
mod_depth <-lme(shannon ~ Depth, random=~1|Replicate/Location, 
                data = na.omit(metadata_cbbL_Abi), na.action=na.omit)  
# change 23S_Cou by other gene x location (cbbL_Cou, cbbL_Abi...)
# change S by shannon
# Look at model results:
anova(mod_depth)
# Check residuals:
hist(mod_depth$residuals)
# Post hoc analysis:
emmeans(mod_depth, pairwise ~ Depth)

# 4. Beta diversity ----
#### | |  23S ----
## Presence/absence > Jaccard
jac_23S <- vegdist(ASV_23S,method = "jaccard",binary = T)
## Abundance > Bray-Curtis
bray_23S <- vegdist(ASV_23S,method="bray")
range(bray_23S)

## NMDS : 
mds_bray_23S <- metaMDS(bray_23S,k=2)
mds_bray_23S$stress
stressplot(mds_bray_23S)
ordihull(mds_bray_23S,groups=metadata_23S$Location,draw="polygon",col = rep(brewer.pal(4,"Set3")),label=F)
orditorp(mds_bray_23S,display="sites",col="black",air=0.01)

# NMDS better plotting
data_scores_23S <- vegan::scores(mds_bray_23S, display = "sites")
data_scores_23S <- as.data.frame(data_scores_23S)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data_scores_23S$site <- metadata_23S$Location  # create a column of site names, from the rownames of data.scores
data_scores_23S$grp <- metadata_23S$Depth  #  add the grp variable created earlier
data_scores_23S$loc_depth <- metadata_23S$Location_Depth

grp_23S_cou1 <- data_scores_23S[data_scores_23S$loc_depth == "Counozouls_1", ][chull(data_scores_23S[data_scores_23S$loc_depth == 
                                                                                                       "Counozouls_1", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_23S_cou2 <-  data_scores_23S[data_scores_23S$loc_depth == "Counozouls_2", ][chull(data_scores_23S[data_scores_23S$loc_depth == 
                                                                                                        "Counozouls_2", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_23S_cou3 <- data_scores_23S[data_scores_23S$loc_depth == "Counozouls_3", ][chull(data_scores_23S[data_scores_23S$loc_depth == 
                                                                                                       "Counozouls_3", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_23S_man1 <- data_scores_23S[data_scores_23S$loc_depth == "Mann_1", ][chull(data_scores_23S[data_scores_23S$loc_depth == 
                                                                                                 "Mann_1", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_23S_man2 <- data_scores_23S[data_scores_23S$loc_depth == "Mann_2", ][chull(data_scores_23S[data_scores_23S$loc_depth == 
                                                                                                 "Mann_2", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_23S_man3 <- data_scores_23S[data_scores_23S$loc_depth == "Mann_3", ][chull(data_scores_23S[data_scores_23S$loc_depth == 
                                                                                                 "Mann_3", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_23S_siik1 <- data_scores_23S[data_scores_23S$loc_depth == "Siikaneva_1", ][chull(data_scores_23S[data_scores_23S$loc_depth == 
                                                                                                       "Siikaneva_1", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_23S_siik2 <- data_scores_23S[data_scores_23S$loc_depth == "Siikaneva_2", ][chull(data_scores_23S[data_scores_23S$loc_depth == 
                                                                                                       "Siikaneva_2", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_23S_siik3 <- data_scores_23S[data_scores_23S$loc_depth == "Siikaneva_3", ][chull(data_scores_23S[data_scores_23S$loc_depth == 
                                                                                                       "Siikaneva_3", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_23S_abi1 <- data_scores_23S[data_scores_23S$loc_depth == "Abisko_1", ][chull(data_scores_23S[data_scores_23S$loc_depth == 
                                                                                                   "Abisko_1", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_23S_abi2 <- data_scores_23S[data_scores_23S$loc_depth == "Abisko_2", ][chull(data_scores_23S[data_scores_23S$loc_depth == 
                                                                                                   "Abisko_2", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_23S_abi3 <- data_scores_23S[data_scores_23S$loc_depth == "Abisko_3", ][chull(data_scores_23S[data_scores_23S$loc_depth == 
                                                                                                   "Abisko_3", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
hull_data_23s <- rbind(grp_23S_cou1, grp_23S_cou2, grp_23S_cou3,
                       grp_23S_man1,grp_23S_man2,grp_23S_man3,
                       grp_23S_siik1,grp_23S_siik2,grp_23S_siik3,
                       grp_23S_abi1, grp_23S_abi2,grp_23S_abi3)  #combine grp.a and grp.b
hull_data_23s

NMDS_23s <- ggplot() + 
  geom_polygon(data = hull_data_23s,aes(x=NMDS1,y=NMDS2,fill=site,group=loc_depth),alpha=0.30) + # add the convex hulls
  geom_point(data = data_scores_23S, aes(x = NMDS1, y = NMDS2, shape=as.factor(grp),colour=site),size=3) + # add the point markers
  coord_equal() +
  theme_bw() + 
  scale_color_manual(values = c("#FC4E07","#059748", "#00AFBB","#B27"))+
  scale_shape_manual(values = c(15,16,17))+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
NMDS_23s

# testing the NMDS with PERMANOVA

dist <- metadata_23S[,c(4,6)]
dist <- sapply(dist, as.numeric)

per_bc <- adonis2(bray_23S ~ dist, permutations = 9999)
per_bc


####| |  cbbL ----
## Presence/absence > Jaccard
jac_cbbL <- vegdist(ASV_cbbL,method = "jaccard",binary = T)
## Abundance > Bray Curtis
bray_cbbL <- vegdist(ASV_cbbL,method="bray")
range(bray_cbbL)
## NMDS : 
mds_bray_cbbL<- metaMDS(bray_cbbL,k=2)
mds_bray_cbbL$stress
stressplot(mds_bray_cbbL)
ordihull(mds_bray_cbbL,groups=metadata_cbbL$Location,draw="polygon",col = rep(brewer.pal(4,"Set3")),label=F)
orditorp(mds_bray_cbbL,display="sites",col="black",air=0.01)

# NMDS better plotting
data_scores_cbbL <- vegan::scores(mds_bray_cbbL, display = "sites")
data_scores_cbbL <- as.data.frame(data_scores_cbbL)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data_scores_cbbL$site <- metadata_cbbL$Location  # create a column of site names, from the rownames of data.scores
data_scores_cbbL$grp <- metadata_cbbL$Depth  #  add the grp variable created earlier
data_scores_cbbL$loc_depth <- metadata_cbbL$Location_Depth

grp_cbbL_cou1 <- data_scores_cbbL[data_scores_cbbL$loc_depth == "Counozouls_1", ][chull(data_scores_cbbL[data_scores_cbbL$loc_depth == 
                                                                                                           "Counozouls_1", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_cbbL_cou2 <-  data_scores_cbbL[data_scores_cbbL$loc_depth == "Counozouls_2", ][chull(data_scores_cbbL[data_scores_cbbL$loc_depth == 
                                                                                                            "Counozouls_2", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_cbbL_cou3 <- data_scores_cbbL[data_scores_cbbL$loc_depth == "Counozouls_3", ][chull(data_scores_cbbL[data_scores_cbbL$loc_depth == 
                                                                                                           "Counozouls_3", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_cbbL_man1 <- data_scores_cbbL[data_scores_cbbL$loc_depth == "Mann_1", ][chull(data_scores_cbbL[data_scores_cbbL$loc_depth == 
                                                                                                     "Mann_1", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_cbbL_man2 <- data_scores_cbbL[data_scores_cbbL$loc_depth == "Mann_2", ][chull(data_scores_cbbL[data_scores_cbbL$loc_depth == 
                                                                                                     "Mann_2", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_cbbL_man3 <- data_scores_cbbL[data_scores_cbbL$loc_depth == "Mann_3", ][chull(data_scores_cbbL[data_scores_cbbL$loc_depth == 
                                                                                                     "Mann_3", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_cbbL_siik1 <- data_scores_cbbL[data_scores_cbbL$loc_depth == "Siikaneva_1", ][chull(data_scores_cbbL[data_scores_cbbL$loc_depth == 
                                                                                                           "Siikaneva_1", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_cbbL_siik2 <- data_scores_cbbL[data_scores_cbbL$loc_depth == "Siikaneva_2", ][chull(data_scores_cbbL[data_scores_cbbL$loc_depth == 
                                                                                                           "Siikaneva_2", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_cbbL_siik3 <- data_scores_cbbL[data_scores_cbbL$loc_depth == "Siikaneva_3", ][chull(data_scores_cbbL[data_scores_cbbL$loc_depth == 
                                                                                                           "Siikaneva_3", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_cbbL_abi1 <- data_scores_cbbL[data_scores_cbbL$loc_depth == "Abisko_1", ][chull(data_scores_cbbL[data_scores_cbbL$loc_depth == 
                                                                                                       "Abisko_1", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_cbbL_abi2 <- data_scores_cbbL[data_scores_cbbL$loc_depth == "Abisko_2", ][chull(data_scores_cbbL[data_scores_cbbL$loc_depth == 
                                                                                                       "Abisko_2", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_cbbL_abi3 <- data_scores_cbbL[data_scores_cbbL$loc_depth == "Abisko_3", ][chull(data_scores_cbbL[data_scores_cbbL$loc_depth == 
                                                                                                       "Abisko_3", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
hull_data_cbbL <- rbind(grp_cbbL_cou1,grp_cbbL_cou2, grp_cbbL_cou3,
                        grp_cbbL_man1,grp_cbbL_man2,grp_cbbL_man3,
                        grp_cbbL_siik1,grp_cbbL_siik2,grp_cbbL_siik3,
                        grp_cbbL_abi1, grp_cbbL_abi2,grp_cbbL_abi3)  #combine grp.a and grp.b
hull_data_cbbL

NMDS_cbbL <- ggplot() + 
  geom_polygon(data = hull_data_cbbL,aes(x=NMDS1,y=NMDS2,fill=site,group=loc_depth),alpha=0.30) + # add the convex hulls
  geom_point(data = data_scores_cbbL, aes(x = NMDS1, y = NMDS2, shape=as.factor(grp),colour=site),size=3) + # add the point markers
  coord_equal() +
  theme_bw() + 
  scale_color_manual(values = c("#FC4E07","#059748", "#00AFBB","#B27"))+
  scale_shape_manual(values = c(15,16,17))+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
NMDS_cbbL

# testing NMDS with PERMANOVA
dist <- metadata_cbbL[,c(4,6)]
dist <- sapply(dist, as.numeric)

per_bc <- adonis2(bray_cbbL ~ dist, permutations = 9999)
per_bc


#### | | bchY ----
## Presence/absence > Jaccard
jac_bchY <- vegdist(ASV_bchY,method = "jaccard",binary = T)
## Abundance > Bray-Curtis
bray_bchY <- vegdist(ASV_bchY,method="bray")
range(bray_bchY)
## NMDS : 
mds_bray_bchY<- metaMDS(bray_bchY,k=2)
mds_bray_bchY$stress
stressplot(mds_bray_bchY)
ordihull(mds_bray_bchY,groups=metadata_bchY$Location,draw="polygon",col = rep(brewer.pal(4,"Set3")),label=F)
orditorp(mds_bray_bchY,display="sites",col="black",air=0.01)

# NMDS better plotting
data_scores_bchY <- vegan::scores(mds_bray_bchY, display = "sites")
data_scores_bchY <- as.data.frame(data_scores_bchY)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data_scores_bchY$site <- metadata_bchY$Location  # create a column of site names, from the rownames of data.scores
data_scores_bchY$grp <- metadata_bchY$Depth  #  add the grp variable created earlier
data_scores_bchY$loc_depth <- metadata_bchY$Location_Depth

grp_bchY_cou1 <- data_scores_bchY[data_scores_bchY$loc_depth == "Counozouls_1", ][chull(data_scores_bchY[data_scores_bchY$loc_depth == 
                                                                                                           "Counozouls_1", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_bchY_cou2 <-  data_scores_bchY[data_scores_bchY$loc_depth == "Counozouls_2", ][chull(data_scores_bchY[data_scores_bchY$loc_depth ==
                                                                                                            "Counozouls_2", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_bchY_cou3 <- data_scores_bchY[data_scores_bchY$loc_depth == "Counozouls_3", ][chull(data_scores_bchY[data_scores_bchY$loc_depth == 
                                                                                                           "Counozouls_3", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_bchY_man1 <- data_scores_bchY[data_scores_bchY$loc_depth == "Mann_1", ][chull(data_scores_bchY[data_scores_bchY$loc_depth ==
                                                                                                     "Mann_1", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_bchY_man2 <- data_scores_bchY[data_scores_bchY$loc_depth == "Mann_2", ][chull(data_scores_bchY[data_scores_bchY$loc_depth == 
                                                                                                     "Mann_2", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_bchY_man3 <- data_scores_bchY[data_scores_bchY$loc_depth == "Mann_3", ][chull(data_scores_bchY[data_scores_bchY$loc_depth == 
                                                                                                     "Mann_3", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_bchY_siik1 <- data_scores_bchY[data_scores_bchY$loc_depth == "Siikaneva_1", ][chull(data_scores_bchY[data_scores_bchY$loc_depth == 
                                                                                                           "Siikaneva_1", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_bchY_siik2 <- data_scores_bchY[data_scores_bchY$loc_depth == "Siikaneva_2", ][chull(data_scores_bchY[data_scores_bchY$loc_depth == 
                                                                                                           "Siikaneva_2", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_bchY_siik3 <- data_scores_bchY[data_scores_bchY$loc_depth == "Siikaneva_3", ][chull(data_scores_bchY[data_scores_bchY$loc_depth == 
                                                                                                           "Siikaneva_3", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_bchY_abi1 <- data_scores_bchY[data_scores_bchY$loc_depth == "Abisko_1", ][chull(data_scores_bchY[data_scores_bchY$loc_depth ==
                                                                                                       "Abisko_1", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_bchY_abi2 <- data_scores_bchY[data_scores_bchY$loc_depth == "Abisko_2", ][chull(data_scores_bchY[data_scores_bchY$loc_depth ==
                                                                                                       "Abisko_2", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp_bchY_abi3 <- data_scores_bchY[data_scores_bchY$loc_depth =="Abisko_3", ][chull(data_scores_bchY[data_scores_bchY$loc_depth ==
                                                                                                      "Abisko_3", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
hull_data_bchY <- rbind(grp_bchY_cou1,grp_bchY_cou2, grp_bchY_cou3,
                        grp_bchY_man1,grp_bchY_man2,grp_bchY_man3,
                        grp_bchY_siik1,grp_bchY_siik2,grp_bchY_siik3,
                        grp_bchY_abi1, grp_bchY_abi2,grp_bchY_abi3)  #combine grp.a and grp.b
hull_data_bchY

NMDS_bchY <- ggplot() + 
  geom_polygon(data = hull_data_bchY,aes(x=NMDS1,y=NMDS2,fill=site,group=loc_depth),alpha=0.30) + # add the convex hulls
  geom_point(data = data_scores_bchY, aes(x = NMDS1, y = NMDS2, shape=as.factor(grp),colour=site),size=3) + # add the point markers
  coord_equal() +
  theme_bw() + 
  scale_color_manual(values = c("#FC4E07","#059748", "#00AFBB","#B27"))+
  scale_shape_manual(values = c(15,16,17))+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
NMDS_bchY

# testing NMDS with PERMANOVA
dist <- metadata_bchY[,c(4,6)]
dist <- sapply(dist, as.numeric)

per_bc <- adonis2(bray_bchY ~ dist, permutations = 9999)
per_bc

##  All graph together ----
p_all_div <- (NMDS_23s | NMDS_cbbL) / (NMDS_bchY)
p_all_div 
p_all_div  + plot_annotation(
  title = 'NMDS according to gene- depth',
  tag_levels = 'a'
)


# 5. Correlation between env data and abundance + index of diversity ----
# creation of a table with everything :

meta_abund_div <- as.data.frame(metadata_23S)

meta_abund_div$shannon_23S <- metadata_23S$shannon
meta_abund_div$shannon_cbbL <- metadata_cbbL$shannon
meta_abund_div$shannon_bchY <- metadata_bchY$shannon

meta_abund_div$NMDS1_23S <- data_scores_23S$NMDS1
meta_abund_div$NMDS1_cbbL <- data_scores_cbbL$NMDS1
meta_abund_div$NMDS1_bchY <- data_scores_bchY$NMDS1

meta_abund_div$NMDS2_23S <- data_scores_23S$NMDS2
meta_abund_div$NMDS2_cbbL <- data_scores_cbbL$NMDS2
meta_abund_div$NMDS2_bchY <- data_scores_bchY$NMDS2

meta_abund_div$S_23S <- metadata_23S$S
meta_abund_div$S_cbbL <- metadata_cbbL$S
meta_abund_div$S_bchY <- metadata_bchY$S

library(corrplot)

meta_abund_div_num <- as.data.frame(sapply(meta_abund_div[,9:67], as.numeric))
meta_abund_div_num <-meta_abund_div_num[,-5]
corrplot(cor(meta_abund_div_num))


## ggplot2 with selected variables: pH, DOC, TN, K, Br, PO4, water cont, phenols, RFE and SST----
library(ggcorrplot)
# Compute a correlation matrix
meta_abund_div_num2 <- meta_abund_div_num[,c(2,3,4,7,13,16,17,20,29,36,39,41:44,47:58)]

corr <- round(cor(meta_abund_div_num2), 1)
head(corr[, 1:6])
tail(corr[, 1:6])

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(meta_abund_div_num2)
head(p.mat[, 1:4])
tail(p.mat[, 1:4])
# Visualize the correlation matrix
# method = "square" (default)
ggcorrplot(corr)

# method = "circle"
ggcorrplot(corr, method = "circle")

# pvalue and Leave blank on no significant coefficient
ggcorrplot(corr, method = "circle",   p.mat = p.mat, hc.order = TRUE,
           insig = "blank", type = "full", show.diag = T)


# simpler graph :
ggcorrplot(corr, method = "circle", hc.order = TRUE,
           type = "full", show.diag = T) +
  scale_fill_gradientn(colors = viridis(256, option = 'D'))


ggcorrplot(corr, method = "circle", hc.order = TRUE,
           type = "full", show.diag = T, colors = c("#3E938BFF","white","darkolivegreen4"))

