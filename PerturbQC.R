# S Master, Children's Hospital of Philadelphia
# masters@chop.edu
# last edit: 2/6/26

library(tidyverse)
library(tidymodels)
library(Biobase)
library(patchwork)
require(xgboost)
require(ranger)
require(C50)

##########
# Load and process gene expression data
##########

# 70-gene breast cancer classifier data obtained from Marchionni et al, PMID: 23682826
# https://pmc.ncbi.nlm.nih.gov/articles/instance/3662649/bin/1471-2164-14-336-S2.gz

load("AdditionalData/objs/glasEset.rda")
load("AdditionalData/objs/buyseEset.rda")

#
# code to select...adapted from Marchionni et al, 2013: supplemental data
select70geneGlas <- featureData(glasEset)$genes70
glasDataAll.es <- glasEset[select70geneGlas, ]
glasCasesAll.ttm <- which(!is.na(pData(glasDataAll.es)$FiveYearMetastasis) &
                            !(pData(glasDataAll.es)$OS < 5 & pData(glasDataAll.es)$TTMevent == 0))
glasDataAll.es <- glasDataAll.es[, glasCasesAll.ttm]
glasDataAll <- as_tibble(t(exprs(glasDataAll.es)))

glasRisk <- pData(glasDataAll.es)$FiveYearMetastasis
glasDataAll <- glasDataAll |>
  add_column(Metastasis = factor(glasRisk), .before=1)

# process Buyse data
select70geneBuyse <- featureData(buyseEset)$genes70
buyseDataAll.es <- buyseEset[select70geneBuyse, ]
buyseCasesAll.ttm <- which(!is.na(pData(buyseDataAll.es)$FiveYearRecurrence))
buyseDataAll.es <- buyseDataAll.es[, buyseCasesAll.ttm]
buyseDataAll <- as_tibble(t(exprs(buyseDataAll.es)))

buyseRisk <- pData(buyseDataAll.es)$FiveYearRecurrence # Risk category

######
# The following code implements three different classification methods:
#.   1) Gradient boosted trees with XGBoost
#.   2) Random Forest
#    3) Gradient boosted trees with C5.0
#
# Tidymodels frameworks are used to train and test
######

##### XGB
mp_xgb <-
  boost_tree(
    trees = 100
  ) |>
  set_engine("xgboost") |>
  set_mode("classification")

mp_xgb_recipe <- recipe(Metastasis ~ ., data=glasDataAll) |> 
  step_zv(all_predictors()) 

mp_xgb_wf <- workflow() |>
  add_model(mp_xgb) |>
  add_recipe(mp_xgb_recipe)

set.seed(1228)
mp_xgb_fit <-
  mp_xgb_wf |>
  fit(glasDataAll)

mp_xgb_pred <-
  mp_xgb_fit |> 
  predict(buyseDataAll, type="prob")

xgb_out <- mp_xgb_pred |> 
  add_column(Outcome=factor(buyseRisk), .before=1) |>
  add_column(PredClass=factor(ifelse(mp_xgb_pred$`.pred_TRUE` >= 0.5, "TRUE", "FALSE")))

xgbroc <- xgb_out |>
  roc_curve(Outcome, .pred_FALSE)

xgbthresh <- xgbroc |>
  filter(sensitivity >= 0.75) |>
  summarise(xgbthres=max(.threshold)) |>
  unlist()

autoplot(xgbroc)

roc_auc(xgb_out, Outcome, .pred_FALSE)

##### RF
mp_rf <-
  rand_forest(
    trees = 1000
  ) |>
  set_engine("ranger") |>
  set_mode("classification")

mp_rf_recipe <- recipe(Metastasis ~ ., data=glasDataAll) |> 
  step_zv(all_predictors()) 

mp_rf_wf <- workflow() |>
  add_model(mp_rf) |>
  add_recipe(mp_rf_recipe)

set.seed(1228)
mp_rf_fit <-
  mp_rf_wf |>
  fit(glasDataAll)

mp_rf_pred <-
  mp_rf_fit |> 
  predict(buyseDataAll, type="prob")

rf_out <- mp_rf_pred |> 
  add_column(Outcome=factor(buyseRisk), .before=1) |>
  add_column(PredClass=factor(ifelse(mp_rf_pred$`.pred_TRUE` >= 0.5, "TRUE", "FALSE")))

rfroc <- rf_out |>
  roc_curve(Outcome, .pred_TRUE, event_level="second")

rfthresh <- rfroc |>
  filter(sensitivity >= 0.75) |>
  summarise(rfthres=max(.threshold)) |>
  unlist()

autoplot(rfroc)

roc_auc(rf_out, Outcome, .pred_FALSE)

## GB with C5.0
mp_gb5 <-
  boost_tree(
    trees = 100
  ) |>
  set_engine("C5.0") |>
  set_mode("classification")

mp_gb5_recipe <- recipe(Metastasis ~ ., data=glasDataAll) |> 
  step_zv(all_predictors()) 

mp_gb5_wf <- workflow() |>
  add_model(mp_gb5) |>
  add_recipe(mp_gb5_recipe)

set.seed(1228)
mp_gb5_fit <-
  mp_gb5_wf |>
  fit(glasDataAll)

mp_gb5_pred <-
  mp_gb5_fit |> 
  predict(buyseDataAll, type="prob")

gb5_out <- mp_gb5_pred |> 
  add_column(Outcome=factor(buyseRisk), .before=1) |>
  add_column(PredClass=factor(ifelse(mp_gb5_pred$`.pred_TRUE` >= 0.5, "TRUE", "FALSE")))

gb5roc <- gb5_out |>
  roc_curve(Outcome, .pred_FALSE)

gb5thresh <- gb5roc |>
  filter(sensitivity >= 0.75) |>
  summarise(gbthres=max(.threshold)) |>
  unlist()

autoplot(gb5roc)

roc_auc(gb5_out, Outcome, .pred_FALSE)

# combine output
xgb_out_t <- xgb_out |> 
  mutate(Model="XGBoost")

rf_out_t <- rf_out |> 
  mutate(Model="Random Forest")

gb5_out_t <- gb5_out |> 
  mutate(Model="Grad Boost (C5.0)")

all_out <- bind_rows(xgb_out_t, rf_out_t, gb5_out_t)

# Supplemental figure 7: combined ROC curves
pdf("Figures/Supp_Figure7.pdf")
all_out |> 
  group_by(Model) |> 
  roc_curve(Outcome, .pred_FALSE) |> 
  autoplot()
dev.off()

##########
# PERTURB CLASSIFIER (for these purposes, we focus on RF and XGB as examples)
#. 1.  Global perturbation
#. 2.  Local perturbation
##########

samps <- length(buyseRisk)
genes <- 70

####
# Global perturbation
####

iters <- 100000

# place to store results
gperturb_xgb_delta <- matrix(rep(0, iters*genes),nrow=iters, ncol=genes)
gperturb_rf_delta <- matrix(rep(0, iters*genes),nrow=iters, ncol=genes)
gperturb_gb5_delta <- matrix(rep(0, iters*genes),nrow=iters, ncol=genes)

set.seed(123)
for (i in 1:genes) { # from 1 to max genes perturbed
  print(str_c("Genes ", i))
  sel_samp <- sample.int(samps, iters, replace=T) # pick samples to perturb
  orig_score_xgb <- mp_xgb_pred$.pred_TRUE[sel_samp] # original predictions, xgb
  orig_score_rf <- mp_rf_pred$.pred_TRUE[sel_samp] # original predictions, rf
  orig_score_gb5 <- mp_gb5_pred$.pred_TRUE[sel_samp] # original predictions, C5.0
  
  perturb_data <- as.matrix(buyseDataAll[sel_samp,]) # grab all of the randomly selected samples
  
  # for each sample, select random genes (set number based on i) and replace with a values from other random samples
  for (j in 1:iters) {
    sel_genes <- sample.int(genes, i)
    sel_wrongsamp <- sample.int(samps-1, i, replace=T)
    drawfrom <- as.matrix(buyseDataAll[-sel_samp[j], ])
    perturb_data[j, sel_genes] <- drawfrom[matrix(c(sel_wrongsamp, sel_genes), ncol=2)]
  }
  
  perturb_data <- as_tibble(perturb_data) # |>
  #    mutate(across(contains("_"), as.numeric))
  
  ppred_xgb <- mp_xgb_fit |> 
    predict(as_tibble(perturb_data), type="prob")
  ppred_rf <- mp_rf_fit |> 
    predict(as_tibble(perturb_data), type="prob")
  ppred_gb5 <- mp_gb5_fit |> 
    predict(as_tibble(perturb_data), type="prob")
  
  gperturb_xgb_delta[,i] <- ppred_xgb$.pred_TRUE - orig_score_xgb
  gperturb_rf_delta[,i] <- ppred_rf$.pred_TRUE - orig_score_rf
  gperturb_gb5_delta[,i] <- ppred_gb5$.pred_TRUE - orig_score_gb5
}

emptyzero <- tibble(`1%`=0.0, `5%`=0.0, `50%`=0.0, `95%`=0.0, `99%`=0.0, Gene=0, `% Perturbed`=0.0)
dxdelta_xgb <- as_tibble(t(apply(gperturb_xgb_delta, 2, \(x) quantile(x, probs=c(0.01, 0.05, 0.5, 0.95, 0.99),
                                                                      type=2)))) |>
  add_column(Gene=1:genes, `% Perturbed` = (100*(1:genes)/genes)) |>
  bind_rows(emptyzero) |>
  add_column(Classifier="XGBoost")

dxdelta_rf <- as_tibble(t(apply(gperturb_rf_delta, 2, \(x) quantile(x, probs=c(0.01, 0.05, 0.5, 0.95, 0.99),
                                                                    type=2)))) |>
  add_column(Gene=1:genes, `% Perturbed` = (100*(1:genes)/genes)) |>
  bind_rows(emptyzero) |>
  add_column(Classifier="Random Forest") 

dxdelta_gb5 <- as_tibble(t(apply(gperturb_gb5_delta, 2, \(x) quantile(x, probs=c(0.01, 0.05, 0.5, 0.95, 0.99),
                                                                      type=2)))) |>
  add_column(Gene=1:genes, `% Perturbed` = (100*(1:genes)/genes)) |>
  bind_rows(emptyzero) |>
  add_column(Classifier="Grad Boost (C5.0)") 

dxdelta <- bind_rows(dxdelta_xgb, dxdelta_rf, dxdelta_gb5)

# ...and visualize; supplemental figure 5
fig8_s <- dxdelta |>
  rename(Median=`50%`) |>
  pivot_longer(`1%`:`99%`, names_to="Quantile", values_to="Change in Probability Score") |>
  mutate(Quantile=factor(Quantile,levels=c("Median","1%","5%","95%","99%"))) |>
  ggplot(aes(x=`% Perturbed`, y=`Change in Probability Score`)) +
  geom_line(aes(lty=Quantile, col=Classifier)) +
  geom_vline(aes(xintercept=(100*2/70)), lty=6) +
  annotate("text", x=3.5, y=0.8, hjust="left", label="2/70 (2.9%) perturbed")

pdf("Figures/Supp_Figure8.pdf", width=10, height=5)
fig8_s
dev.off()

# show underlying population distributions, supplemental figure 6

fig9a_s <- xgb_out |>
  mutate(`Prediction score`=.pred_TRUE) |>
  ggplot(aes(`Prediction score`)) +
  geom_histogram(bins=30) +
  labs(title="XGBoost") +
  #  xlim(-0.05,1.05) +
  geom_vline(aes(xintercept=xgbthresh), lwd=0.75) +
  geom_vline(aes(xintercept=xgbthresh+unlist(dxdelta_xgb[2,"1%"])), color="green", lty=2, lwd=0.75) +
  geom_vline(aes(xintercept=xgbthresh+unlist(dxdelta_xgb[2,"99%"])), color="green", lty=2, lwd=0.75) +
  geom_vline(aes(xintercept=xgbthresh+unlist(dxdelta_xgb[2,"5%"])), color="blue", lty=2, lwd=0.75) +
  geom_vline(aes(xintercept=xgbthresh+unlist(dxdelta_xgb[2,"95%"])), color="blue", lty=2, lwd=0.75)

fig9b_s <- rf_out |>
  mutate(`Prediction score`=.pred_TRUE) |>
  ggplot(aes(`Prediction score`)) +
  geom_histogram(bins=30) +
  labs(title="Random Forest") +
  #  xlim(-0.1,1.1) +
  geom_vline(aes(xintercept=rfthresh), lwd=0.75) +
  geom_vline(aes(xintercept=rfthresh+unlist(dxdelta_rf[2,"1%"])), color="green", lty=2, lwd=0.75) +
  geom_vline(aes(xintercept=rfthresh+unlist(dxdelta_rf[2,"99%"])), color="green", lty=2, lwd=0.75) +
  geom_vline(aes(xintercept=rfthresh+unlist(dxdelta_rf[2,"5%"])), color="blue", lty=2, lwd=0.75) +
  geom_vline(aes(xintercept=rfthresh+unlist(dxdelta_rf[2,"95%"])), color="blue", lty=2, lwd=0.75)


fig9c_s <- gb5_out |>
  mutate(`Prediction score`=.pred_TRUE) |>
  ggplot(aes(`Prediction score`)) +
  geom_histogram(bins=30) +
  labs(title="Grad Boost (C5.0)") +
  #  xlim(-0.05,1.05) +
  geom_vline(aes(xintercept=gb5thresh), lwd=0.75) +
  geom_vline(aes(xintercept=gb5thresh+unlist(dxdelta_gb5[2,"1%"])), color="green", lty=2, lwd=0.75) +
  geom_vline(aes(xintercept=gb5thresh+unlist(dxdelta_gb5[2,"99%"])), color="green", lty=2, lwd=0.75) +
  geom_vline(aes(xintercept=gb5thresh+unlist(dxdelta_gb5[2,"5%"])), color="blue", lty=2, lwd=0.75) +
  geom_vline(aes(xintercept=gb5thresh+unlist(dxdelta_gb5[2,"95%"])), color="blue", lty=2, lwd=0.75)

pdf("Figures/Supp_Figure9.pdf", width=8, height=3)
fig9a_s + fig9b_s + fig9c_s
dev.off()

# find % of cases that are still diagnosable at 98% confidence with ~3% analyte error rate
dx_98_xgb <- xgb_out |>
  filter((`.pred_TRUE` < xgbthresh+unlist(dxdelta_xgb[2,"1%"])) |
           (`.pred_TRUE` > xgbthresh+unlist(dxdelta_xgb[2,"99%"]))) |>
  dim()
dx_98_rf <- rf_out |>
  filter((`.pred_TRUE` < rfthresh+unlist(dxdelta_rf[2,"1%"])) |
           (`.pred_TRUE` > rfthresh+unlist(dxdelta_rf[2,"99%"]))) |>
  dim()
dx_98_gb5 <- gb5_out |>
  filter((`.pred_TRUE` < gb5thresh+unlist(dxdelta_gb5[2,"1%"])) |
           (`.pred_TRUE` > gb5thresh+unlist(dxdelta_gb5[2,"99%"]))) |>
  dim()

dx_98_xgb[1] / 307
dx_98_rf[1] / 307
dx_98_gb5[1] / 307

# ...and 90% confidence (with ~3% analyte error rate)
dx_90_xgb <- xgb_out |>
  filter((`.pred_TRUE` < xgbthresh+unlist(dxdelta_xgb[2,"5%"])) |
           (`.pred_TRUE` > xgbthresh+unlist(dxdelta_xgb[2,"95%"]))) |>
  dim()
dx_90_rf <- rf_out |>
  filter((`.pred_TRUE` < rfthresh+unlist(dxdelta_rf[2,"5%"])) |
           (`.pred_TRUE` > rfthresh+unlist(dxdelta_rf[2,"95%"]))) |>
  dim()
dx_90_gb5 <- gb5_out |>
  filter((`.pred_TRUE` < gb5thresh+unlist(dxdelta_gb5[2,"5%"])) |
           (`.pred_TRUE` > gb5thresh+unlist(dxdelta_gb5[2,"95%"]))) |>
  dim()

dx_90_xgb[1] / 307
dx_90_rf[1] / 307
dx_90_gb5[1] / 307

####
# Local perturbation
####

# we can do this with random draws, but in our current case (307 samples) it's tractable to
# exhaustively look at all possible replacements for each gene

# record class probability changes #####$**********
output_xgb_d <- matrix(vector(mode="numeric", length=(genes*samps*(samps-1))), ncol=genes)
output_gb5_d <- matrix(vector(mode="numeric", length=(genes*samps*(samps-1))), ncol=genes)
output_rf_d <- matrix(vector(mode="numeric", length=(genes*samps*(samps-1))), ncol=genes)

for (i in 1:samps) {
  print(str_c("Iter ", i))
  orig <- buyseDataAll[rep(i, samps-1),]
  for (j in 1:genes) {
    repvalues <- unlist(buyseDataAll[,j])[-i]
    perturb <- orig
    perturb[,j] <- repvalues
    idx <- (i-1)*(samps-1)
    output_xgb_d[idx+(1:(samps-1)), j] <- predict(mp_xgb_fit, perturb, type="prob")$`.pred_FALSE` - mp_xgb_pred$`.pred_FALSE`[i]
    output_gb5_d[idx+(1:(samps-1)), j] <- predict(mp_gb5_fit, perturb, type="prob")$`.pred_FALSE` - mp_gb5_pred$`.pred_FALSE`[i]
    output_rf_d[idx+(1:(samps-1)), j] <- predict(mp_rf_fit, perturb, type="prob")$`.pred_FALSE` - mp_rf_pred$`.pred_FALSE`[i]
  }
}

######### Process for display
df_rf_d <- as_tibble(output_rf_d)
names(df_rf_d) <- names(buyseDataAll)
df_rf_d <- pivot_longer(df_rf_d, everything(), names_to="Gene", values_to="Delta")

df_gb5_d <- as_tibble(output_gb5_d)
names(df_gb5_d) <- names(buyseDataAll)
df_gb5_d <- pivot_longer(df_gb5_d, everything(), names_to="Gene", values_to="Delta")

df_xgb_d <- as_tibble(output_xgb_d)
names(df_xgb_d) <- names(buyseDataAll)
df_xgb_d <- pivot_longer(df_xgb_d, everything(), names_to="Gene", values_to="Delta")


df_xgb_d_t <- df_xgb_d |> 
  mutate(Classifier = "XGBoost")

df_gb5_d_t <- df_gb5_d |> 
  mutate(Classifier = "Grad Boost (C5.0)")

df_rf_d_t <- df_rf_d |> 
  mutate(Classifier = "Random Forest")

df_d_t <- bind_rows(df_rf_d_t, df_gb5_d_t, df_xgb_d_t)

# With outliers
#
# ggplot(df_d_t, aes(x=Gene, y=Delta)) +
#   geom_boxplot(aes(color=Classifier), outlier.alpha = 0.03, outlier.size=1) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Without outliers: Figure 5
fig5 <- df_d_t |>
  mutate(Gene=str_remove(Gene, "_.*")) |>
  ggplot(aes(x=Gene, y=Delta)) +
  geom_boxplot(aes(color=Classifier), outliers=F) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("Figures/Figure5.pdf",width=10, height=5)
fig5
dev.off()

save.image(file="Results/PerturbQCresults.RData")
