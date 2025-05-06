# Diversity, abundance, and biogeography of CO2-fixing microorganisms in peatlands
# Marie Le Geay, Kyle Mayers, Anna Sytiuk, Ellen Dorrepaal, Martin Küttim, 
# Mariusz Lamentowicz, Eeva-Stiina Tuittila, Béatrice Lauga, and Vincent E.J. Jassey

# Author: Marie Le Geay and Vincent Jassey

################################################################################
# 6. MFA

# Main text:
## Figure 3

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
library("FactoMineR"); library("factoextra"); library(lmom); library(vegan)
library(ggforce); library(scatterpie)

# |   |  Load data ----
setwd("C:/Users/mlegeay/OneDrive/Documents/chapter_2/Redaction/Révision/RStudio-to-publish/Data")
ASV_23S <- read.table("t_ASV_23S.csv", header = T, sep = ",")
rownames(ASV_23S) <- ASV_23S$X
ASV_23S <- ASV_23S[,-1]

ASV_cbbL <- read.table("t_ASV_cbbL.csv", header = T, sep = ",")
rownames(ASV_cbbL) <- ASV_cbbL$X
ASV_cbbL <- ASV_cbbL[,-1]

ASV_bchY <- read.table("t_ASV_bchY.csv", header = T, sep = ",")
rownames(ASV_bchY) <- ASV_bchY$X
ASV_bchY <- ASV_bchY[,-1]

cluster_data <- read.table("JSDMxMFA.csv", head = T, sep = ",")


# 1. MFA on all AsvS ----
# | | Data preparation ----
matrix_23S <- decostand(ASV_23S, method = "hellinger")
matrix_cbbL <- decostand(ASV_cbbL, method = "hellinger")
matrix_bchY <- decostand(ASV_bchY, method = "hellinger")

tab.win = cbind(matrix_23S,matrix_cbbL,matrix_bchY)

gr = c(ncol(matrix_23S), ncol(matrix_cbbL),  ncol(matrix_bchY))
gr

# | | Run the MFA and plot ----
allmfa = MFA(tab.win, group=gr, type=c("c","c","c"), ncp=2, 
             name.group=c("ASV_23S","ASV_cbbL","ASV_bchY")) # all plots

# Eigenvalues and % of variation
allmfa$eig
ev = allmfa$eig[,1]
names(ev) = 1:nrow(allmfa$eig)
evplot(ev)

### RV coefficients (correlations between groups of variables) with tests
rvp = allmfa$group$RV
rvp[1,2] = coeffRV(matrix_23S, matrix_cbbL)$p.value
rvp[1,3] = coeffRV(matrix_23S,matrix_bchY)$p.value

rvp[2,3] = coeffRV(matrix_cbbL, matrix_bchY)$p.value

rvp

# Site scores 
sit.sco = allmfa$global.pca$ind$coord

#  Hierarchical clustering from site scores 
mfa.de = dist(sit.sco)
mfa.ew = hclust(mfa.de, "ward.D2")
plot(mfa.ew)

# Plot the dendrogram
k = 6
gr = cutree(mfa.ew, k)

# | | Plot ----
pl = ordiplot(sit.sco[,1:2],  main="MFA new analysis")
abline(h=0, lty=3)
abline(v=0, lty=3)
points(sit.sco[,1:2], cex=2, col=1+c(1:k)[gr], pch=14+c(1:k)[gr])
text(sit.sco[,1:2], rownames(sit.sco[,1:2]), pos=4, cex=0.7)
# Add clustering dendrogram
ordicluster(pl, mfa.ew, col="dark grey")
# Add a legend for groups
legend("bottomleft", paste("Group",c(1:k)), pch=14+c(1:k), col=1+c(1:k), pt.cex=2)


# 2. Pie chart plot based on JSDMs within the MFA ----
# Merge PCA and cluster data
MFA_data <- data.frame(
  site = rownames(sit.sco),         # Site names
  PC1 = sit.sco[,1],                # PCA coordinates for PC1
  PC2 = sit.sco[,2]                # PCA coordinates for PC2
)
rownames(MFA_data) <- rownames(cluster_data)
MFA_cluster_data <- cbind(MFA_data, cluster_data)
MFA_cluster_data <- MFA_cluster_data[,-4]

# Reshape data for plotting pies
pca_cluster_long <- MFA_cluster_data %>%
  pivot_longer(cols = starts_with("PresC1"),
               names_to = "PresC",
               values_to = "Percentage")

# Plot PCA with pie charts
ggplot(MFA_cluster_data, aes(x = PC1, y = PC2)) +
  geom_point(color = "white") + 
  geom_hline(yintercept=0, linetype="dashed", color = "grey")+
  geom_vline(xintercept=0, linetype="dashed", color = "grey")+
  geom_scatterpie(aes(x = PC1, y = PC2, r = 0.1), 
                  data = MFA_cluster_data, 
                  cols = c("PresC1", "PresC2", "PresC3", "PresC4", "PresC5", "PresC6")) +
  coord_fixed() +
  labs(title = "PCA with Cluster Percentages as Pie Charts") +
  theme_classic()

