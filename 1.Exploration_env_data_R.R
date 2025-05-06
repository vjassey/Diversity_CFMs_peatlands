# Diversity, abundance, and biogeography of CO2-fixing microorganisms in peatlands
# Marie Le Geay, Kyle Mayers, Anna Sytiuk, Ellen Dorrepaal, Martin Küttim, 
# Mariusz Lamentowicz, Eeva-Stiina Tuittila, Béatrice Lauga, and Vincent E.J. Jassey

# Author: Marie Le Geay and Vincent Jassey

################################################################################
# 1. Exploration of environmental data 

# Main text:
## Figure 1

# Supplementary Information:
## Tables S5, S6 and S7 
## Figures S2, S3, S4, S5 and S6

################################################################################
# |   |  Load packages ----

library(ggplot2); library(plyr); library(reshape2); library(grid); library(gridExtra); library(ape); library(picante)
library(seqinr); library(phytools); library(scales); library(nlme); library(vegan); library(RColorBrewer)
library("BiocManager"); library("Biobase"); library("rhdf5"); library("devtools"); library(phyloseq); library(msa)
library(ggtree); library(tidyjson); library(RColorBrewer);library(tidyverse); library(viridis); library(randomForest); library(caret); library(cluster)
library(tidyverse); library(cluster);  library(parallelDist); library(boral);  library(corrplot); library(igraph); library(car); library(fossil); library(ROCR); library(doMC)
library(effectsize); library(phyloseq); library(Biostrings); library(dplyr); library(tidyr); library(bioseq)
library(vegan); library(ade4); library(BAT); library(patchwork); library(hrbrthemes); library(gcookbook)
library(htmltools); library(HH);library(FactoMineR);library(factoextra)
library(MASS);library(reshape2);library(cowplot);library(nlme);library(emmeans);library(multcomp)

# |   |  Load data ----
getwd() # See in which folder you are located
setwd("C:/Users/mlegeay/OneDrive/Documents/chapter_2/Redaction/Révision/RStudio-to-publish/Data")# set up a new path

metadata <- read.table("Metadata.csv", header = T, sep = ";")
metadata$Location <- fct_relevel(metadata$Location, 
                                 c("Counozouls","Mann",
                                   "Siikaneva","Abisko"))


# 1. Visualisation of metadata ----
# | | Global characteristics ----
g_c1 <- ggplot(data = metadata, aes(x = Location, y = pH, fill = as.factor(Location))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="pH")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_c1 

g_c2 <- ggplot(data = metadata, aes(x = Location, y = DOC, fill = as.factor(Location))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="DOC")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_c2

g_c3 <- ggplot(data = metadata, aes(x = Location, y = TN, fill = as.factor(Location))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="TN")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_c3

g_c4 <- ggplot(data = metadata, aes(x = Location, y = WTD, fill = as.factor(Location))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="WTD")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_c4

p_all_climate <- (g_c1 |g_c2)  / (g_c3 |g_c4)
p_all_climate
p_all_climate  + plot_annotation(
  title = 'Global characteristics',
  tag_levels = 'a'
)

# | | Nutrients ----
g_1 <- ggplot(data = metadata, aes(x = Location, y = Na, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="Na")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_1

g_2 <- ggplot(data = metadata, aes(x = Location, y = NH4, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="NH4")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_2

g_3 <- ggplot(data = metadata, aes(x = Location, y = K, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="K")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_3

g_4 <- ggplot(data = metadata, aes(x = Location, y = Mg, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="Mg")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_4

g_5 <- ggplot(data = metadata, aes(x = Location, y = Ca, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="Ca")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_5

g_6 <- ggplot(data = metadata, aes(x = Location, y = F, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="F")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_6

g_7 <- ggplot(data = metadata, aes(x = Location, y = Cl, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="Cl")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_7

g_8 <- ggplot(data = metadata, aes(x = Location, y = NO2, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="NO2")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_8

g_9 <- ggplot(data = metadata, aes(x = Location, y = Br, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="Br")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_9

g_10 <- ggplot(data = metadata, aes(x = Location, y = NO3, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="NO3")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_10

g_11 <- ggplot(data = metadata, aes(x = Location, y = SO4, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="SO4")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_11

g_12 <- ggplot(data = metadata, aes(x = Location, y = PO4, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="PO4")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
g_12

p_all_nutrient2 <- (g_1 |g_2 |g_3)  / (g_4 |g_5 | g_6)  / (g_7| g_8| g_9) / (g_10| g_11| g_12)
p_all_nutrient2
p_all_nutrient2  + plot_annotation(
  title = 'Nutrients',
  tag_levels = 'a'
)


# | |  Sphagnum metabolites ----
gm_1 <- ggplot(data = metadata, aes( x = Location, y = carbohydrates, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="carbohydrates")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
gm_1

gm_2 <- ggplot(data = metadata, aes( x = Location, y = flavonoids, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="flavonoids")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
gm_2

gm_3 <- ggplot(data = metadata, aes( x = Location, y = phenols, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="phenols")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
gm_3

gm_4 <- ggplot(data = metadata, aes( x = Location, y = tanins, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="tannins")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
gm_4

gm_5 <- ggplot(data = metadata, aes( x = Location, y = water_phenols, fill = as.factor(Depth))) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="water_phenols")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
gm_5

p_all_metabolites  <- (gm_1 |gm_2) / (gm_3|gm_4) / (gm_5/plot_spacer()) 
p_all_metabolites  + plot_annotation(
  title = 'Metabolites',
  tag_levels = 'a'
)


# | |  OM quality ----
gq_1 <- ggplot(data = metadata, aes( x = Location, y = DOCq, fill = Location)) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="DOCq")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
gq_1

gq_2 <- ggplot(data = metadata, aes( x = Location, y = Peak_A, fill = Location)) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="Peak_A")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
gq_2

gq_3 <- ggplot(data = metadata, aes( x = Location, y = Peak_C, fill = Location)) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="Peak_C")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
gq_3

gq_4 <- ggplot(data = metadata, aes( x = Location, y = Peak_M, fill = Location)) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="Peak_M")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
gq_4

gq_5 <- ggplot(data = metadata, aes( x = Location, y = RFE, fill = Location)) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="RFE")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
gq_5

gq_6 <- ggplot(data = metadata, aes( x = Location, y = Freshness, fill = Location)) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="Freshness")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
gq_6

gq_7 <- ggplot(data = metadata, aes( x = Location, y = BIX, fill = Location)) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="BIX")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
gq_7

gq_8 <- ggplot(data = metadata, aes( x = Location, y = FI, fill = Location)) + 
  geom_boxplot(alpha = 0.4)+
  labs(title="FI")+
  scale_fill_brewer(palette = 10)+
  theme_classic()
gq_8

p_all_docq  <- (gq_1 |gq_2 |gq_3 | gq_4) / (gq_5 |gq_6 |gq_7 | gq_8) 
p_all_docq  + plot_annotation(
  title = 'docq',
  tag_levels = 'a'
)


################################################################################
# 2. Add statistical analyses ----
# Linear mixed models with lme funciton of nlme package
# Replicate by location as random term

# | | Effect of location ----

# Compute the model:
mod_loc <-lme(DOCq ~ Location, random=~1|Replicate/Location, 
              data = na.omit(metadata), na.action=na.omit) 
# change Na by other va (NH4,F, tanins...)
# Look at model results:
anova(mod_loc) 
# Check residuals:
hist(mod_loc$residuals)
# Post hoc analysis:
emmeans(mod_loc, pairwise ~ Location)


# | | Effect of depth ----
metadata$Depth <- as.factor(metadata$Depth)
metadata_Cou <- subset(metadata, metadata$Location == "Counozouls")
metadata_Man <- subset(metadata, metadata$Location == "Mann")
metadata_Siik <- subset(metadata, metadata$Location == "Siikaneva")
metadata_Abi <- subset(metadata, metadata$Location == "Abisko")
# Compute the model for each location:
mod_depth <-lme(tanins ~ Depth, random=~1|Replicate/Location, 
                data = na.omit(metadata_Cou), na.action=na.omit) 

# change Na by other va (NH4,F, tanins...)
# change _Cou by _Man, _Siik or _Abi
# Look at model results:
anova(mod_depth)
# Check residuals:
hist(mod_depth$residuals)
# Post hoc analysis:
emmeans(mod_depth, pairwise ~ Depth)


################################################################################
# 3. PCA analysis----
# | | Normalisation of the data ----
metadata_sel <- metadata[,-c(48,51:55)]
metadata_SD <- decostand(metadata_sel[,12:49], "standardize",na.rm = TRUE)

# | | Compute PCA ----
# PCA with the function PCA of the FactoMineR package

pca1 <- PCA(metadata_SD, graph = F)
pca1$ind
pca1$eig # to get the score of each axes (second column, third column = sum of axes)

# extract pc scores for first two component and add to dataframe
metadata$pc1 <- pca1$ind$coord[, 1] # indexing the first column
metadata$pc2 <- pca1$ind$coord[, 2]  # indexing the second column
metadata$pc3 <- pca1$ind$coord[, 3]  # indexing the third column

# extract the data for the variable contributions to each of the pc axes
pca.vars <- pca1$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars.m <- melt(pca.vars, id.vars = "vars")

# | | Plot the PCA ----
# By convention, the variable contribution plot has a circle around the variables that has a radius of 1
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

circ <- circleFun(c(0,0),2,npoints = 500)


p <- ggplot(data = metadata, 
            aes(x = pc1, y = pc2, color = Location, shape = as.factor(Depth))) +
  scale_shape_manual(values = c(15,16,17))+
  scale_color_manual(values = c("#FC4E07","#059748", "#00AFBB","#B27"))+
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(alpha = 0.8) 
p

p2 <- p + stat_ellipse(geom="polygon", aes(fill = Location), 
                       alpha = 0.2, 
                       show.legend = FALSE, 
                       level = 0.95) +
  xlab("PC 1 (28.07%)") + 
  ylab("PC 2 (20.13%)") +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))
p2


vars.p <-  ggplot() +
  geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_segment(data = pca.vars, aes(x = 0, xend = Dim.1, y = 0, yend = Dim.2),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
               lwd = 1) + 
  geom_text(data = pca.vars, 
            aes(x = Dim.1, y =  Dim.2, 
                label = c("WTD", "pH", "DOC", "	TN", "	Na", "NH4", "K",
                          "Mg", "Ca", "F","Cl","NO2", "Br", "NO3",
                          "SO4","PO4","water_content", "carbohydrates", "flavonoids",
                          "phenols","tanins","water_phenols", "DOCq","Peak_A","Peak_C",
                          "Peak_M","FI","Freshness","RFE","BIX","winter_air_temp",
                          "winter_prec","winter_WTD","winter_soil_temp","spring_air_temp",
                          "spring_prec","spring_WTD","spring_soil_temp")), 
            check_overlap = F, size = 3) +
  xlab("PC 1") + 
  ylab("PC2") +
  coord_equal() +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))
vars.p


# 4. Selection of variables for further analysis ----
ASV_23S <- read.table("t_ASV_23S.csv", header = T, sep = ",")
rownames(ASV_23S) <- ASV_23S$X
ASV_23S <- ASV_23S[,-1]
bray_23S <- vegdist(ASV_23S,method="bray")

# | | Pre-selection of parameters ----
# will select only one or two
library(corrplot)
corrplot(cor(metadata[,13:51]), type = "lower", diag = F)

par(mfrow=c(2,3))
for (i in 13:51)hist(metadata[,i],main=colnames(metadata)[i])

# | | Correction of data distribution ----
metadata$K <- log(metadata$K)
metadata$Mg <- log(metadata$Mg)
metadata$Ca <- log(metadata$Ca)
metadata$Br <- sqrt(metadata$Br)
metadata$NO3 <- sqrt(metadata$NO3)
metadata$PO4 <- sqrt(metadata$PO4)
metadata$carbohydrates <- sqrt(metadata$carbohydrates)
metadata$phenols <- log(metadata$phenols)
metadata$water_phenols <- log(metadata$water_phenols)
metadata$F <- sqrt(metadata$F)
metadata$Cl <- sqrt(metadata$Cl)


# New table to keep variables not too colinear
metadata_2 <- metadata[,c(1,2,3,4,6,8,9,10,11,12,13,14,15,
                                    16,17,18,19,20,21,22,24,27,28,29,31,39,40,41,42,47)]
corrplot(cor(metadata_2[,10:30]), type = "lower", diag = F)

# | | Selection and corrplot ----
# Selected variables based on all previous results:
# pH, DOC, TN, K, Br, PO4, Water content, Phenols, RFE and SST
metadata_3 <- metadata_2[,c(11,12,13,16,21,22,23,25,27,30)]
corrplot(cor(metadata_3), type = "lower", diag = T)
