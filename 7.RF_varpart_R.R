# Diversity, abundance, and biogeography of CO2-fixing microorganisms in peatlands
# Marie Le Geay, Kyle Mayers, Anna Sytiuk, Ellen Dorrepaal, Martin Küttim, 
# Mariusz Lamentowicz, Eeva-Stiina Tuittila, Béatrice Lauga, and Vincent E.J. Jassey

# Author: Marie Le Geay and Vincent Jassey

################################################################################
# 7. Random Forest and Variance partitioning

# Main text:
## Figure 6

################################################################################
# |   |  Load packages ----
library(ggplot2); library(plyr); library(reshape2); library(grid); library(gridExtra); library(ape); library(picante)
library(ggtree); library(tidyjson); library(RColorBrewer);library(tidyverse); library(viridis); library(randomForest); library(caret); library(cluster)
library(tidyverse); library(cluster);  library(parallelDist); library(boral);  library(corrplot); library(igraph); library(car); library(fossil); library(ROCR); library(doMC)
library(effectsize); library(phyloseq); library(Biostrings); library(dplyr); library(tidyr); library(bioseq)
library(vegan); library(ade4); library(BAT); library(patchwork); library(hrbrthemes); library(gcookbook)
library(data.table); library (kableExtra); library(randomForestExplainer); library(spatialRF)

# |   |  Load data ----
setwd("C:/Users/mlegeay/OneDrive/Documents/chapter_2/Redaction/Révision/RStudio-to-publish/Data")
data_RF <- read.table("RF.csv", head = T, sep = ";")
data_varpart <- read.table("varpart.csv", head = T, sep = ";")

# 1. RF for total CFM abundance ----
# |   | Data
data_RF2 <- data_RF[,c(5:15,25,33:37,40,41,45)]
data_RF2$tot <- log(data_RF2$AQ_total)
data_RF_std <- decostand(data_RF2, method = "standardize")

# |   |  Variable names
dependent.variable.name <- "tot"
predictor.variable.names <- colnames(data_RF_std[,-c(13,21)])

# |   | Fitting an initial model to define the preference order of the predictors in the multicollinearity analysis
m <- spatialRF::rf(
  data = data_RF_std,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  ranger.arguments = list(
    num.trees = 999
  )
)

# |   | Repeat
m.repeat <- rf_repeat(
  model = m,
  repetitions = 999,
  #seed = random.seed,
  verbose = FALSE
)

spatialRF::plot_importance(
  m.repeat,
  verbose = FALSE
)

# |   | p.value
importance.df <- randomForestExplainer::measure_importance(
  m,
  measures = c("mean_min_depth", "no_of_nodes", "times_a_root", "p_value")
)

kableExtra::kbl(
  importance.df %>%
    dplyr::arrange(mean_min_depth) %>%
    dplyr::mutate(p_value = round(p_value, 4)),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# |   | variance partitioning

res = varpart(data_RF_std$tot, 
              ~ MFA_1+MFA_2+ S + shannon,
              ~ C1+C2+C3+C4,
              ~ DOC + K + pH + TN + PO4+ Br+ phenols + tannins,
              ~ WTD + spring_soil_temp + spring_prec, 
              data = data_RF_std)
res
plot(res)

ggplot(data_varpart)+
  geom_col(aes(x = Name, y = tot, fill = Name))


# 2. RF for oxygenic phototrophs abundance ----
# |   | Data
data_RF3 <- data_RF[,c(5:15,22,30,34:37,40,41,42)]
data_RF3$ox <- log(data_RF3$AQ_23S)
data_RF_ox_std <- decostand(data_RF3, method = "standardize")

# |   |  Variable names
dependent.variable.name <- "ox"
predictor.variable.names <- colnames(data_RF_ox_std[,-c(13,21)])

# |   | Fitting an initial model to define the preference order of the predictors in the multicollinearity analysis
m <- spatialRF::rf(
  data = data_RF_ox_std,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  ranger.arguments = list(
    num.trees = 999
  )
)

# |   | Repeat
m.repeat <- rf_repeat(
  model = m,
  repetitions = 999,
  #seed = random.seed,
  verbose = FALSE
)

spatialRF::plot_importance(
  m.repeat,
  verbose = FALSE
)

# |   | p.value
importance.df <- randomForestExplainer::measure_importance(
  m,
  measures = c("mean_min_depth", "no_of_nodes", "times_a_root", "p_value")
)

kableExtra::kbl(
  importance.df %>%
    dplyr::arrange(mean_min_depth) %>%
    dplyr::mutate(p_value = round(p_value, 4)),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# |   | variance partitioning
res_23S = varpart(data_RF_ox_std$ox, 
              ~ MFA_1+MFA_2+ S_23S + shannon_23S,
              ~ C1+C2+C3+C4,
              ~ DOC + K + pH + TN + PO4+ Br+ phenols + tannins,
              ~ WTD + spring_soil_temp + spring_prec, 
              data = data_RF_ox_std)
res_23S
plot(res)

ggplot(data_varpart)+
  geom_col(aes(x = Name, y = ox, fill = Name))


# 3. RF for chemoautotrophs abundance ----
# |   | Data
data_RF4 <- data_RF[,c(5:15,23,31,34:37,40,41,43)]
data_RF4$chem <- log(data_RF4$AQ_cbbL)
data_RF_chem_std <- decostand(data_RF4, method = "standardize")

# |   |  Variable names
dependent.variable.name <- "chem"
predictor.variable.names <- colnames(data_RF_chem_std[,-c(13,21)])

# |   | Fitting an initial model to define the preference order of the predictors in the multicollinearity analysis
m <- spatialRF::rf(
  data = data_RF_chem_std,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  ranger.arguments = list(
    num.trees = 999
  )
)

# |   | Repeat
m.repeat <- rf_repeat(
  model = m,
  repetitions = 999,
  #seed = random.seed,
  verbose = FALSE
)

spatialRF::plot_importance(
  m.repeat,
  verbose = FALSE
)

# |   | p.value
importance.df <- randomForestExplainer::measure_importance(
  m,
  measures = c("mean_min_depth", "no_of_nodes", "times_a_root", "p_value")
)

kableExtra::kbl(
  importance.df %>%
    dplyr::arrange(mean_min_depth) %>%
    dplyr::mutate(p_value = round(p_value, 4)),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# |   | variance partitioning
res_cbbL = varpart(data_RF_chem_std$chem, 
                  ~ MFA_1+MFA_2+ S_cbbL + shannon_cbbL,
                  ~ C1+C2+C3+C4,
                  ~ DOC + K + pH + TN + PO4+ Br+ phenols + tannins,
                  ~ WTD + spring_soil_temp + spring_prec, 
                  data = data_RF_chem_std)
res_cbbL
plot(res_cbbL)

ggplot(data_varpart)+
  geom_col(aes(x = Name, y = chem, fill = Name))


# 4. RF for AAnPBs abundance ----
# |   | Data
data_RF5 <- data_RF[,c(5:15,24,32,34:37,40,41,44)]
data_RF5$an <- log(data_RF5$AQ_bchY)
data_RF_an_std <- decostand(data_RF5, method = "standardize")

# |   |  Variable names
dependent.variable.name <- "an"
predictor.variable.names <- colnames(data_RF_an_std[,-c(13,21)])

# |   | Fitting an initial model to define the preference order of the predictors in the multicollinearity analysis
m <- spatialRF::rf(
  data = data_RF_an_std,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  ranger.arguments = list(
    num.trees = 999
  )
)

# |   | Repeat
m.repeat <- rf_repeat(
  model = m,
  repetitions = 999,
  #seed = random.seed,
  verbose = FALSE
)

spatialRF::plot_importance(
  m.repeat,
  verbose = FALSE
)

# |   | p.value
importance.df <- randomForestExplainer::measure_importance(
  m,
  measures = c("mean_min_depth", "no_of_nodes", "times_a_root", "p_value")
)

kableExtra::kbl(
  importance.df %>%
    dplyr::arrange(mean_min_depth) %>%
    dplyr::mutate(p_value = round(p_value, 4)),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# |   | variance partitioning
res_bchY = varpart(data_RF_an_std$an, 
                   ~ MFA_1+MFA_2+ S_pufM + shannon_bchY,
                   ~ C1+C2+C3+C4,
                   ~ DOC + K + pH + TN + PO4+ Br+ phenols + tannins,
                   ~ WTD + spring_soil_temp + spring_prec, 
                   data = data_RF_an_std)
res_bchY
plot(res_cbbL)

ggplot(data_varpart)+
  geom_col(aes(x = Name, y = an, fill = Name))
