# Diversity, abundance, and biogeography of CO2-fixing microorganisms in peatlands
# Marie Le Geay, Kyle Mayers, Anna Sytiuk, Ellen Dorrepaal, Martin Küttim, 
# Mariusz Lamentowicz, Eeva-Stiina Tuittila, Béatrice Lauga, and Vincent E.J. Jassey

# Author: Marie Le Geay and Vincent Jassey

################################################################################
# 5. JSDM

# Main text:
## Figure 5

# Supplementary Information:
## Figures S15 and S16

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
library(future.apply); library(jSDM); library(dplyr); library(igraph); library(tidyverse); library(netplot); library(vegan)
library(magrittr)

# |   |  Load data ----
setwd("C:/Users/mlegeay/OneDrive/Documents/chapter_2/Redaction/Révision/RStudio-to-publish/Data")

# load R environment as it took ages for this computer to run the JSDMs
# >>> JSD.RData
data23s <- read.delim("ASV_23S_rar.csv", sep = ";")
data23s <- t(data23s)
data23s <- as.data.frame(data23s)
colnames(data23s) <- data23s[1,]
data23s <- data23s[-1,]

databchY <- read.delim("ASV_bchY_rar.csv", sep = ";")
databchY <- t(databchY)
databchY <- as.data.frame(databchY)
colnames(databchY) <- databchY[1,]
databchY <- databchY[-1,]

datacbbL <- read.delim("ASV_cbbL_rar.csv", sep = ";")
datacbbL <- t(datacbbL)
datacbbL <- as.data.frame(datacbbL)
colnames(datacbbL) <- datacbbL[1,]
datacbbL <- datacbbL[-1,]

env <- read.delim("Metadata.csv", sep = ';')
env <- env[-61,]
env.sel <- env %>% dplyr::select(pH, DOC, K, spring_soil_temp, RFE, TN, phenols, PO4, tanins, water_content, spring_WTD)
env.sel <- sapply(env.sel, as.numeric)
env.sel.sd <- decostand(env.sel, 'standardize')

# 1. Data preparation ----
# Merge ASV files
allasv <- cbind(data23s, databchY, datacbbL)
allasv_num <- sapply(allasv, as.numeric)

# Calculate the relative abundance of each species
rela_ab <- prop.table(allasv_num, 1)*100

## Filter species with a relative abundance > 1%
allasv_sorted <- allasv_num[,apply(rela_ab, 2, function(x) any(x>0.1))]
str(allasv_sorted)


# 2. jSDM modelling ----
X = env.sel.sd
Y = allasv_sorted

#= Number of sites
nsite <- 60

#= Set seed for repeatability
seed <- 1234
set.seed(seed)

#= Number of species
nsp <- 1309

#= Number of latent variables
n_latent <- 2
np <- ncol(allasv_sorted)

#= Latent variables W
W <- matrix(rnorm(nsite*n_latent,0,1), nsite, n_latent)

#= Fixed species effect beta 
beta.target <- t(matrix(runif(nsp*np,-1,1),
                        byrow=TRUE, nrow=nsp))

#= Factor loading lambda  
lambda.target <- matrix(0, n_latent, nsp)
mat <- t(matrix(runif(nsp*n_latent, -1, 1), byrow=TRUE, nrow=nsp))
lambda.target[upper.tri(mat, diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
diag(lambda.target) <- runif(n_latent, 0, 2)

#= Variance of random site effect 
V_alpha.target <- 0.2

#= Random site effect alpha
alpha.target <- rnorm(nsite,0 , sqrt(V_alpha.target))

# |   | run jSDM (DO NOT RUN IF NOT ENOUGHT MEMORY ON THE COMPUTER)
mod <- jSDM_gaussian(# Iteration
  burnin=2000,
  mcmc=20000,
  thin=100,
  # Response variable
  response_data=Y,
  # Explanatory variables
  site_formula=~.,
  site_data = X,
  n_latent=2,
  site_effect="random",
  # Starting values
  alpha_start=0,
  beta_start=0,
  lambda_start=0,
  W_start=0,
  V_alpha=1,
  V_start=1 ,
  # Priors
  shape_Valpha=0.5, rate_Valpha=0.0005,
  shape_V=0.5, rate_V=0.0005,
  mu_beta=0, V_beta=1,
  mu_lambda=0, V_lambda=1,
  seed=1234, verbose=1)


mod<-jSDM_binomial_probit(# Iteration
  burnin=500,
  mcmc=1000,
  thin=10,
  # Response variable
  presence_data=decostand(Y, 'pa'),
  # Explanatory variables
  site_formula=~.,
  site_data = X,
  n_latent=2,
  site_effect="random",
  # Starting values
  alpha_start=0,
  beta_start=0,
  lambda_start=0,
  W_start=0,
  V_alpha=1,
  # Priors
  shape_Valpha=0.5,
  rate_Valpha=0.0005,
  mu_beta=0, V_beta=1,
  mu_lambda=0, V_lambda=1,
  seed=1234, verbose=1)

# |   | Outputs ----

#= Parameter estimates
oldpar <- par(no.readonly = TRUE)
## beta_j
mean_beta <- matrix(0,nsp,ncol(env.sel.sd))
pdf(file=file.path(tempdir(), "Posteriors_beta_jSDM_probit.pdf"))
par(mfrow=c(ncol(env.sel.sd),2))
for (j in 1:nsp) {
  mean_beta[j,] <- apply(mod$mcmc.sp[[j]]
                         [,1:ncol(env.sel.sd)], 2, mean)  
  for (p in 1:ncol(env.sel.sd)){
    coda::traceplot(mod$mcmc.sp[[j]][,p])
    coda::densplot(mod$mcmc.sp[[j]][,p],
                   main = paste(colnames(mod$mcmc.sp[[j]])[p],", species : ",j))
    abline(v=beta.target[p,j],col='red')
  }
}
dev.off()

## lambda_j
mean_lambda <- matrix(0,nsp,n_latent)
pdf(file=file.path(tempdir(), "Posteriors_lambda_jSDM_probit.pdf"))
par(mfrow=c(n_latent*2,2))
for (j in 1:nsp) {
  mean_lambda[j,] <- apply(mod$mcmc.sp[[j]]
                           [,(ncol(X)+1):(ncol(X)+n_latent)], 2, mean)  
  for (l in 1:n_latent) {
    coda::traceplot(mod$mcmc.sp[[j]][,ncol(X)+l])
    coda::densplot(mod$mcmc.sp[[j]][,ncol(X)+l],
                   main=paste(colnames(mod$mcmc.sp[[j]])
                              [ncol(X)+l],", species : ",j))
    abline(v=lambda.target[l,j],col='red')
  }
}
dev.off()

## W latent variables
par(mfrow=c(1,2))
for (l in 1:n_latent) {
  plot(W[,l],
       summary(mod$mcmc.latent[[paste0("lv_",l)]])[[1]][,"Mean"],
       main = paste0("Latent variable W_", l),
       xlab ="obs", ylab ="fitted")
  abline(a=0,b=1,col='red')
}

## alpha
par(mfrow=c(1,3))
plot(alpha.target, summary(mod$mcmc.alpha)[[1]][,"Mean"],
     xlab ="obs", ylab ="fitted", main="site effect alpha")
abline(a=0,b=1,col='red')

## Valpha
coda::traceplot(mod$mcmc.V_alpha)
coda::densplot(mod$mcmc.V_alpha)
abline(v=V_alpha.target,col='red')

## Variance of residuals
par(mfrow=c(1,2))
coda::traceplot(mod$mcmc.V)
coda::densplot(mod$mcmc.V,
               main="Variance of residuals")
abline(v=V.target, col='red')

## Deviance
summary(mod$mcmc.Deviance)
plot(mod$mcmc.Deviance)

#= Predictions
par(mfrow=c(1,1))
plot(allasv_sorted, mod$Y_pred,
     main="Response variable",xlab="obs",ylab="fitted")
abline(a=0,b=1,col='red')

par(oldpar)


# |   | Beta coefficients from the model and make a mean through iterations ----
betacoefs <- lapply(mod$mcmc.sp, function(x) colMeans(x))
bcoefs <- do.call(rbind, betacoefs)
rownames(bcoefs) <- colnames(allasv_sorted)
head(bcoefs)
bcoefs <- as.data.frame(bcoefs)
bcoefs$CI <- rownames(bcoefs)
hist(bcoefs$beta_spring_soil_temp)


# |   | Environmental correlations ----
# Also very long
envirocor <- get_enviro_cor(mod)


# |   | Residual correlations
residcor <- get_residual_cor(mod)


## |   | Build co-occurrence network ----
Matrix_cor_full <- as.data.frame(envirocor$cor.sig)
diag(Matrix_cor_full) = 0

# 2.bis Run the code from here after having loading the environment if computer not powerful enought ----
Matrix_cor <- data.frame(row=rownames(Matrix_cor_full)[row(Matrix_cor_full)[upper.tri(Matrix_cor_full)]], 
                         col=colnames(Matrix_cor_full)[col(Matrix_cor_full)[upper.tri(Matrix_cor_full)]], 
                         corr=Matrix_cor_full [upper.tri(Matrix_cor_full)])

colnames(Matrix_cor)[1:3] <- c("from","to","weight")

Matrix_cor_cut <- Matrix_cor[Matrix_cor$weight>0.6,]
hist(Matrix_cor_cut$weight) #define the threshold

g <- graph_from_data_frame(Matrix_cor_cut, directed = F)
mst <- igraph::mst(g)
layout  <- layout.kamada.kawai(mst)

#network plot ----

wc <- cluster_fast_greedy(g)

#assign names to vertex
test <- as.data.frame(g[1])
test$names <- rownames(test)
V(g)$names <- test$names

#create a vector with averaged relative abundance of each ASV in the network
rela_ab2 <- prop.table(allasv_sorted, 1) ##Calculate the relative abundance of each species
relab <- as.data.frame(colMeans(rela_ab2))
relab$names <- row.names(relab)
colnames(relab) <- c("relab", "names")

relab_net <- inner_join(relab, test)
range(relab_net$relab)

# plotting the network
plot(wc, g, vertex.size = 4)


# |   |Extract data on clusters ----
clusters_network <- as.data.frame(cbind(wc$names, wc$membership)) #create the matrix with cluster and ASVs
colnames(clusters_network) <- c("CI", "cluster")
mem <- clusters_network %>% group_by(as.factor(cluster)) %>% tally() #overview of the cluster size

cluster_comp <- bcoefs %>%
  inner_join(clusters_network, key=CI)
cluster_comp$cluster <- as.numeric(cluster_comp$cluster)
str(cluster_comp)

cluster_comp2 <- subset(cluster_comp, cluster < 7) # six major groups!

car::Anova(aov(cluster_comp2$beta_DOC~as.factor(cluster_comp2$cluster))) #stats with unbalanced size
TukeyHSD(aov(cluster_comp2$beta_DOC~as.factor(cluster_comp2$cluster)))


# |   | Beta-coefs better graphs per cluster ----
bcoefs_mean<- cluster_comp2 %>% dplyr::select(-CI) %>%
  group_by(cluster) %>%
  na.omit() %>% 
  summarise(across(everything(), c(mean)), .groups = 'drop') %>% as.data.frame() %>% gather(group)

std <- function(x) sd(x)/sqrt(length(x))
bcoefs_sd<- cluster_comp2 %>% dplyr::select(-CI) %>%
  group_by(cluster) %>%
  na.omit() %>% 
  summarise(across(everything(), std), .groups = 'drop') %>% as.data.frame() %>% gather(group)

bcoefs_stats <- bcoefs_mean
bcoefs_stats$sd <- bcoefs_sd$value
colnames(bcoefs_stats) <- c("driver", "mean", "sd")
bcoefs_stats$cluster <- rep(1:6, 15)
bcoefs_stats <- subset(bcoefs_stats, driver != 'cluster' & driver != 'lambda_1_1' & driver != 'lambda_2_1' & driver != 'beta_(Intercept)_1')
levels(droplevels(as.factor(bcoefs_stats$driver)))

cols <- c("darkgreen","cadetblue","#C42503FF","#E1DD37FF", "purple","#FEA632FF")

ggplot(bcoefs_stats, aes(x=driver, y=mean, group=cluster)) + 
  geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd, colour = factor(cluster), group = cluster)) + scale_color_manual(values = cols) + theme_classic() + geom_hline(yintercept = -0,
                                                                                                                                                                 linetype = 3,
                                                                                                                                                                 lwd = 0.5) + facet_wrap(~cluster, scales = "free") + coord_flip()
# 3. Part of network at each site ----
# Size of each cluster:
length(wc[[1]])
length(wc[[2]])
length(wc[[3]])
length(wc[[4]])
length(wc[[5]])
length(wc[[6]])
length(wc[[7]])
length(wc[[1]])+length(wc[[2]])+length(wc[[3]])+length(wc[[4]])+length(wc[[5]])+length(wc[[6]])+length(wc[[7]])
length(wc[[1]])+length(wc[[2]])+length(wc[[3]])+length(wc[[4]])

# for each site how many ASVs by cluster ?
spe.pa <- decostand(allasv_sorted, "pa") #matrix of PA

Clusters = data.frame(
  Treatment = env$ID_23S,
  PresC1=rowSums(spe.pa[,wc[[1]]])/rowSums(spe.pa),
  PresC2=rowSums(spe.pa[,wc[[2]]])/rowSums(spe.pa),
  PresC3=rowSums(spe.pa[,wc[[3]]])/rowSums(spe.pa),
  PresC4=rowSums(spe.pa[,wc[[4]]])/rowSums(spe.pa),
  PresC5=rowSums(spe.pa[,wc[[5]]])/rowSums(spe.pa),
  PresC6=rowSums(spe.pa[,wc[[6]]])/rowSums(spe.pa)
)

Clusters_count = data.frame(
  Treatment = env$ID_23S,
  PresC1=rowSums(spe.pa[,wc[[1]]]),
  PresC2=rowSums(spe.pa[,wc[[2]]]),
  PresC3=rowSums(spe.pa[,wc[[3]]]),
  PresC4=rowSums(spe.pa[,wc[[4]]]),
  PresC5=rowSums(spe.pa[,wc[[5]]]),
  PresC6=rowSums(spe.pa[,wc[[6]]])
)

ClustersMH=Clusters %>%
  group_by(Treatment) %>%
  summarise(
    PC1=mean(PresC1),
    PC2=mean(PresC2),
    PC3=mean(PresC3),
    PC4=mean(PresC4),
    PC5=mean(PresC5),
    PC6=mean(PresC6),
    PC7=1 - (PC1 + PC2 + PC3 + PC4 + PC4 + PC5 + PC6)
  )

# Graphs
ClustP = ClustersMH[,c(1:8)] %>%
  pivot_longer(c(2:8), names_to = "Cluster", values_to = "Presence")

g=ggplot(ClustP, aes(x=Treatment, y=Presence, fill=Cluster))
g + geom_bar(stat="identity") + ylab('')+
  xlab('')+
  coord_cartesian(ylim=c(0,1))+
  theme(
    axis.line = element_line(colour = "black",size = 0.5, linetype = "solid"),
    panel.background=element_blank(), #rect(fill='lavenderblush3', colour='black', size=1, linetype='solid'),
    panel.grid.major=element_blank(), #element_line(size=0.4, colour='grey', linetype="solid"),
    panel.grid.minor=element_blank(), #element_line(size=0.4, colour='white', linetype="solid"),
    panel.border=element_rect(fill=NA, colour='grey'),
    strip.text=element_text(size=13, face='italic'),
    #strip.background = element_blank(),
    axis.title.x=element_text(size=14),
    axis.title.y=element_text(size=14),
    axis.text.x=element_text(size=12, color='black'),
    axis.text.y=element_text(size=12, color='black'))

# group count ----
getwd()
count <- read.table("group_count.csv", head = T, sep = ",")
g2 = ggplot(count, aes(x=Cluster, y=relab, fill=Gene))
g2 + geom_bar(stat="identity") + ylab('')+
  xlab('')+
  coord_cartesian(ylim=c(0,1))

#  4. Testing the residues ----
Matrix_cor_res <- as.data.frame(residcor$cor.sig)
diag(Matrix_cor_res) = 0

Matrix_res <- data.frame(row=rownames(Matrix_cor_res)[row(Matrix_cor_res)[upper.tri(Matrix_cor_res)]], 
                         col=colnames(Matrix_cor_res)[col(Matrix_cor_res)[upper.tri(Matrix_cor_res)]], 
                         corr=Matrix_cor_res [upper.tri(Matrix_cor_res)])

colnames(Matrix_res)[1:3] <- c("from","to","weight")

Matrix_res_cut <- Matrix_res[Matrix_res$weight>0.2,]
hist(Matrix_res_cut$weight) #define the threshold
hist(Matrix_res$weight)

g_res <- graph_from_data_frame(Matrix_res_cut, directed = F)
mst_res <- igraph::mst(g)
layout_res  <- layout.kamada.kawai(mst)

#network plot

wc_res <- cluster_fast_greedy(g_res)
plot(wc_res, g_res, vertex.size = 4, vertex.label=NA)


# 5. Relative abundance by cluster ----
relab_net
clusters_network

names(relab_net)[2] <- "CI"

relativ_abund_cluster <- inner_join(relab_net,clusters_network)

Taxo_23S <- read.table("Taxo_23S.csv", header = T, sep = ",")
rownames(Taxo_23S) <- Taxo_23S$X
Taxo_23S <- Taxo_23S[,-1]

Taxo_cbbL <- read.table("Taxo_cbbL.csv", header = T, sep = ",")
rownames(Taxo_cbbL) <- Taxo_cbbL$X
Taxo_cbbL <- Taxo_cbbL[,-1]

Taxo_bchY <- read.table("Taxo_bchY.csv", header = T, sep = ",")
rownames(Taxo_bchY) <- Taxo_bchY$X
Taxo_bchY <- Taxo_bchY[,-1]

Taxo_23S$ID <- rownames(Taxo_23S)

abund_comp_clust <- read.table("relab_comm_comp_clusters.csv", head = T , sep = ",")

g=ggplot(abund_comp_clust, aes(x=cluster, y=relab, fill=class))
g + geom_bar(stat="identity") + ylab('')+
  xlab('')+
  coord_cartesian(ylim=c(0,1))+
  theme(
    axis.line = element_line(colour = "black",size = 0.5, linetype = "solid"),
    panel.background=element_blank(), #rect(fill='lavenderblush3', colour='black', size=1, linetype='solid'),
    panel.grid.major=element_blank(), #element_line(size=0.4, colour='grey', linetype="solid"),
    panel.grid.minor=element_blank(), #element_line(size=0.4, colour='white', linetype="solid"),
    panel.border=element_rect(fill=NA, colour='grey'),
    strip.text=element_text(size=13, face='italic'),
    #strip.background = element_blank(),
    axis.title.x=element_text(size=14),
    axis.title.y=element_text(size=14),
    axis.text.x=element_text(size=12, color='black'),
    axis.text.y=element_text(size=12, color='black')) +
  facet_wrap(~Gene, nrow = 3)

