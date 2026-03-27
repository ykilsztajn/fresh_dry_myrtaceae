############################################################
# Libraries
############################################################

library(ape)
library(dplyr)
library(tidyr)
library(randomForest)
library(gbm)
library(caret)

############################################################
# Load data
############################################################

df <- read.csv("spp_mean_data.csv")
trees <- read.tree("trees_with_inserted_species.tre")

############################################################
# Species names
############################################################

df$tip_name <- ifelse(
  is.na(df$superfilo_name),
  df$spp,
  df$superfilo_name
)

############################################################
# Remove weight traits
############################################################

df <- df %>%
  dplyr::select(-matches("weigth"))

############################################################
# Climate variables
############################################################

vars_climate <- c(
  "temp_annual",
  "temp_sazon",
  "precip_annual",
  "precip_sazon"
)

############################################################
# Trait base names
############################################################

trait_cols <- grep("(_fresh$|_dry$)", names(df), value=TRUE)

vars_base <- unique(gsub("_fresh|_dry","",trait_cols))

############################################################
# Function: phylogenetic eigenvectors
############################################################

get_phylo_vectors <- function(tree, species){
  
  tree <- drop.tip(tree, setdiff(tree$tip.label, species))
  
  dist_phy <- cophenetic(tree)
  
  pcoa <- ape::pcoa(dist_phy)
  
  var <- cumsum(pcoa$values$Relative_eig)
  
  k <- which(var >= 0.8)[1]
  k <- min(k,10)
  
  vec <- as.data.frame(pcoa$vectors[,1:k])
  
  colnames(vec) <- paste0("Axis",1:k)
  
  vec$tip_name <- rownames(vec)
  
  return(vec)
  
}

############################################################
# Function: dataset builder
############################################################

build_dataset <- function(df, trait, model){
  
  var_fresh <- paste0(trait,"_fresh")
  var_dry <- paste0(trait,"_dry")
  
  if(!all(c(var_fresh,var_dry) %in% names(df))) return(NULL)
  
  df2 <- df %>%
    rename(
      fresh = !!var_fresh,
      dry = !!var_dry
    )
  
  if(model=="M1"){
    
    dados <- df2 %>%
      dplyr::select(fresh,dry)
    
  }
  
  if(model=="M2"){
    
    dados <- df2 %>%
      dplyr::select(fresh,dry,starts_with("Axis"))
    
  }
  
  if(model=="M3"){
    
    dados <- df2 %>%
      dplyr::select(fresh,dry,all_of(vars_climate),starts_with("Axis"))
    
  }
  
  dados <- na.omit(dados)
  
  if(nrow(dados) < 30) return(NULL)
  
  return(dados)
  
}

############################################################
# Function: metrics
############################################################

compute_metrics <- function(obs,pred){
  
  rmse <- sqrt(mean((pred-obs)^2))
  
  r2 <- cor(pred,obs)^2
  
  range_obs <- max(obs) - min(obs)
  
  if(range_obs == 0) range_obs <- mean(obs)
  
  nrmse <- rmse / range_obs
  nrmse_mean <- rmse / mean(obs)
  
  percent_error <- nrmse * 100
  
  list(
    rmse = rmse,
    nrmse = nrmse,
    range = range_obs,
    nrmse_mean = nrmse_mean,
    mean = mean(obs),
    perc_error = percent_error,
    r2 = r2
  )
  
}

############################################################
# Random Forest
############################################################

run_rf <- function(df, trait, model){
  
  dados <- build_dataset(df,trait,model)
  
  y <- dados$fresh
  X <- dados %>% dplyr::select(-fresh)
  
  set.seed(123)
  train_id <- createDataPartition(y,p=0.8,list=FALSE)
  
  rf <- randomForest(
    x=X[train_id,,drop=FALSE],
    y=y[train_id],
    ntree=1000,
    mtry=max(1,floor(sqrt(ncol(X))))
  )
  
  pred <- predict(rf,X[-train_id,,drop=FALSE])
  obs <- y[-train_id]
  
  m <- compute_metrics(obs,pred)
  
  metrics <- data.frame(
    trait=trait,
    model=model,
    method="RF",
    rmse=m$rmse,
    nrmse=m$nrmse,
    range = m$range,
    nrmse_mean = m$nrmse_mean,
    mean = m$mean,
    percent_error=m$perc_error,
    r2=m$r2
  )
  
  return(metrics)
  
}

############################################################
# Gradient Boosting
############################################################

run_gbm <- function(df,trait,model){
  
  dados <- build_dataset(df,trait,model)
  if(is.null(dados)) return(NULL)
  
  y <- dados$fresh
  
  set.seed(123)
  train_id <- createDataPartition(y,p=0.8,list=FALSE)
  
  gb <- gbm(
    formula=fresh~.,
    data=dados[train_id,],
    distribution="gaussian",
    n.trees=3000,
    interaction.depth=3,
    shrinkage=0.01,
    bag.fraction=0.7,
    verbose=FALSE
  )
  
  pred <- predict(gb,dados[-train_id,],n.trees=3000)
  obs <- y[-train_id]
  
  m <- compute_metrics(obs,pred)
  
  metrics <- data.frame(
    trait=trait,
    model=model,
    method="GBM",
    rmse=m$rmse,
    nrmse=m$nrmse,
    range = m$range,
    nrmse_mean = m$nrmse_mean,
    mean = m$mean,
    percent_error=m$perc_error,
    r2=m$r2
  )
  
  return(metrics)
  
}

############################################################
# Run models across trees
############################################################

models <- c("M1","M2","M3")

results_metrics <- list()

for(i in seq_along(trees)){
  
  cat("Tree",i,"\n")
  
  tree <- trees[[i]]
  
  phylo_vectors <- get_phylo_vectors(tree,df$tip_name)
  
  df_model <- df %>%
    left_join(phylo_vectors,by="tip_name")
  
  for(trait in vars_base){
    
    for(model in models){
      
      res_rf <- tryCatch(run_rf(df_model,trait,model),error=function(e) NULL)
      res_gbm <- tryCatch(run_gbm(df_model,trait,model),error=function(e) NULL)
      
      for(res in list(res_rf,res_gbm)){
        
        if(!is.null(res)){
          
          res$tree <- i
          
          results_metrics[[length(results_metrics)+1]] <- res
          
        }
        
      }
      
    }
    
  }
  
}

results_metrics <- bind_rows(results_metrics)

############################################################
# Summaries
############################################################

summary_metrics <- results_metrics %>%
  group_by(trait,model,method) %>%
  summarise(
    mean_rmse=mean(rmse,na.rm=TRUE),
    sd_rmse=sd(rmse,na.rm=TRUE),
    mean_nrmse_range=mean(nrmse,na.rm=TRUE),
    sd_nrmse_range=sd(nrmse,na.rm=TRUE),
    mean_range=mean(range,na.rm=TRUE),
    mean_nrmse_mean=mean(nrmse_mean,na.rm=TRUE),
    sd_nrmse_mean=sd(nrmse_mean,na.rm=TRUE),
    mean_mean=mean(mean,na.rm=TRUE),
    mean_percent_error=mean(percent_error,na.rm=TRUE),
    sd_percent_error=sd(percent_error,na.rm=TRUE),
    mean_r2=mean(r2,na.rm=TRUE),
    sd_r2=sd(r2,na.rm=TRUE),
    r2_low=quantile(r2,0.025,na.rm=TRUE),
    r2_high=quantile(r2,0.975,na.rm=TRUE),
    .groups="drop"
  )

############################################################
# Export
############################################################

write.csv(summary_metrics,"trait_model_metrics.csv",row.names=FALSE)


#############################################################
#Making plots
##############################################################

library(ggplot2)
library(dplyr)

summary_metrics <- read.csv("trait_model_metrics.csv")
# Ordem dos traits
ordem_traits <- rev(c(
  "leaf_area", "leaf_perimeter", "leaf_length", "leaf_width",
  "leaf_length_per_width", "leaf_circularity_index", "leaf_lma_area",
  "flower_petal_size", "flower_length", "flower_diameter_petal",
  "flower_diameter_stamen", "flower_hipant_length", "flower_hipant_diameter",
  "fruit_length", "fruit_diameter", "seed_length", "seed_diameter"
))

# Preparar tabela
summary_metrics <- summary_metrics %>%
  filter(trait %in% ordem_traits) %>%
  mutate(
    trait = factor(trait, levels = ordem_traits),
    model_algo = paste(method, model, sep = " - ")
  )

# ---------- NRMSE pelo mean ----------
ggplot(summary_metrics, aes(x = trait, y = mean_nrmse_mean, color = model_algo)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = mean_nrmse_mean - sd_nrmse_mean,
                    ymax = mean_nrmse_mean + sd_nrmse_mean),
                width = 0.2, position = position_dodge(width = 0.5), alpha = 0.7) +
  coord_flip() +
  theme_classic() +
  labs(x = NULL, y = "NRMSE (RMSE / trait mean)", color = "Algorithm - Model") +
  scale_color_manual(values = c(
    "RF - M1" = "orange", "RF - M2" = "red", "RF - M3" = "darkred",
    "GBM - M1" = "lightblue", "GBM - M2" = "blue", "GBM - M3" = "darkblue"
  )) +
  geom_hline(yintercept = c(0.1,0.2,0.3), linetype = "dashed", alpha = 0.4)

# ---------- NRMSE pelo range ----------
ggplot(summary_metrics, aes(x = trait, y = mean_nrmse_range, color = model_algo)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = mean_nrmse_range - sd_nrmse_range,
                    ymax = mean_nrmse_range + sd_nrmse_range),
                width = 0.2, position = position_dodge(width = 0.5), alpha = 0.7) +
  coord_flip() +
  theme_classic() +
  labs(x = NULL, y = "NRMSE (RMSE / trait range)", color = "Algorithm - Model") +
  scale_color_manual(values = c(
    "RF - M1" = "orange", "RF - M2" = "red", "RF - M3" = "darkred",
    "GBM - M1" = "lightblue", "GBM - M2" = "blue", "GBM - M3" = "darkblue"
  )) +
  geom_hline(yintercept = c(0.1,0.2,0.3), linetype = "dashed", alpha = 0.4)

# ---------- R² ----------
ggplot(summary_metrics, aes(x = trait, y = mean_r2, fill = model_algo)) +
  geom_col(position = position_dodge(width = 0.7), width = 1, alpha = 0.9) +
  geom_errorbar(aes(ymin = mean_r2 - sd_r2, ymax = mean_r2 + sd_r2),
                width = 0.7, position = position_dodge(width = 0.7), alpha = 0.7) +
  coord_flip() +
  theme_classic() +
  labs(x = NULL, y = expression(R^2), fill = "Algorithm - Model") +
  scale_fill_manual(values = c(
    "RF - M1" = "orange", "RF - M2" = "red", "RF - M3" = "darkred",
    "GBM - M1" = "lightblue", "GBM - M2" = "blue", "GBM - M3" = "darkblue"
  )) +
  geom_hline(yintercept = c(0.25,0.5,0.75), linetype = "dashed", alpha = 0.4)

