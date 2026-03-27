library(dplyr)
library(stringr)
library(ggplot2)
library(broom)
library(tidyr)
library(MuMIn)
library(tidyverse)
library(ape)
library(caper)
library(rr2)
library(phylolm)

#######################################################
# Testing PGLS models to explain herbarium effect
#######################################################

# load data
df_media_spp <- read.csv("spp_mean_data.csv")

# load trees
trees <- read.tree("trees_with_inserted_species.tre")

# create column matching phylogeny tip labels
df_media_spp$phylo_name <- ifelse(
  is.na(df_media_spp$superfilo_name),
  df_media_spp$spp,
  df_media_spp$superfilo_name
)

# Structures to be analyzed
estruturas <- c("flower", "fruit", "seed", "leaf")

efeitos_por_estrutura <- list(
  flower = "_herb_effect",
  fruit  = "_herb_effect",
  seed   = "_herb_effect",
  leaf   = "_herb_effect"
)

#remove weigth_herb_effect
df_media_spp <- df_media_spp %>%
  dplyr::select(-flower_weigth_herb_effect,
         -fruit_weigth_herb_effect,
         -seed_weigth_herb_effect,
         -leaf_weigth_herb_effect)

# storage
resultados_modelos_herb <- list()
tabela_modelos_herb <- data.frame()

# loop over structures
for (estrutura in estruturas) {
  
  var_resposta <- grep(
    paste0("^", estrutura, ".*", efeitos_por_estrutura[[estrutura]], "$"),
    names(df_media_spp),
    value = TRUE
  )
  
  var_wp <- paste0(estrutura, "_water_percentage")
  
  if (!(var_wp %in% names(df_media_spp))) next
  
  for (var in var_resposta) {
    
    dados_brutos <- df_media_spp[, c(var, var_wp, "genus", "spp", "phylo_name")]
    
    colnames(dados_brutos)[1:2] <- c("herb_effect", "water_percentage")
    
    dados <- na.omit(dados_brutos)
    
    if (nrow(dados) > 10) {
      
      # vectors to store results across trees
      slopes <- c()
      pvals <- c()
      r2s <- c()
      r2s_lik <- c()
      
      for (tree in trees) {
        
        # prune tree to species present
        tree_pruned <- drop.tip(
          tree,
          setdiff(tree$tip.label, dados$phylo_name)
        )
        
        dados_tree <- dados[match(tree_pruned$tip.label, dados$phylo_name), ]
        
        rownames(dados_tree) <- dados_tree$phylo_name
        
        mod <- try(
          phylolm(
            herb_effect ~ water_percentage,
            data = dados_tree,
            phy = tree_pruned,
            model = "lambda"
          ),
          silent = TRUE
        )
        
        if (inherits(mod, "try-error")) next
        
        if(var(dados_tree$herb_effect)==0) next
        if(var(dados_tree$water_percentage)==0) next
        
        coef_table <- summary(mod)$coefficients
        
        slopes <- c(slopes, coef_table["water_percentage", "Estimate"])
        pvals  <- c(pvals, coef_table["water_percentage", "p.value"])
        
        r2_val <- R2_resid(mod, phy = tree_pruned)
        r2_val_lik <- R2_lik(mod)
        
        r2s <- c(r2s, r2_val)
        r2s_lik <- c(r2s_lik, r2_val_lik)
      }
      
      if (length(slopes) > 0) {
        
        resumo_modelos <- data.frame(
          estrutura = estrutura,
          efeito = "herb",
          var_resposta = var,
          n_trees = length(slopes),
          slope_mean = mean(slopes),
          slope_sd = sd(slopes),
          pvalue_mean = mean(pvals),
          r2res_mean = mean(r2s),
          r2res_sd = sd(r2s),
          r2lik_mean = mean(r2s_lik),
          r2lik_sd = sd(r2s_lik),
          stringsAsFactors = FALSE
        )
        
        tabela_modelos_herb <- bind_rows(
          tabela_modelos_herb,
          resumo_modelos
        )
      }
    }
  }
}

# export results
write.csv(
  tabela_modelos_herb,
  "result_herb_effect_pgls_models.csv",
  row.names = FALSE
)


#########################################################
###herb effect ~water percentage plots
#######################################################

library(dplyr)
library(ggplot2)

### adicionar coluna de significância

tabela_modelos_herb <- tabela_modelos_herb %>%
  mutate(
    signif = ifelse(pvalue_mean < 0.05, "significant", "ns")
  )


### construir dados para plot

dados_plot <- list()

for (estrutura in estruturas) {
  
  var_resposta <- grep(
    paste0("^", estrutura, ".*", efeitos_por_estrutura[[estrutura]], "$"),
    names(df_media_spp),
    value = TRUE
  )
  
  var_wp <- paste0(estrutura, "_water_percentage")
  
  if (!(var_wp %in% names(df_media_spp))) next
  
  for (var in var_resposta) {
    
    dados_temp <- df_media_spp[, c(var, var_wp)]
    
    colnames(dados_temp) <- c("herb_effect", "water_percentage")
    
    dados_temp$trait <- var
    
    dados_plot[[var]] <- dados_temp
  }
}

dados_plot <- bind_rows(dados_plot)

dados_plot <- dados_plot %>%
  filter(!is.na(herb_effect), !is.na(water_percentage))


### criar linhas de regressão médias

linhas_modelo <- list()

for(i in 1:nrow(tabela_modelos_herb)){
  
  var <- tabela_modelos_herb$var_resposta[i]
  slope <- tabela_modelos_herb$slope_mean[i]
  slope_sd <- tabela_modelos_herb$slope_sd[i]
  sig <- ifelse(tabela_modelos_herb$pvalue_mean[i] < 0.05,"significant","ns")
  
  dados_trait <- dados_plot %>%
    filter(trait == var)
  
  if(nrow(dados_trait) == 0) next
  
  x_mean <- mean(dados_trait$water_percentage)
  y_mean <- mean(dados_trait$herb_effect)
  
  intercept <- y_mean - slope * x_mean
  
  wp_seq <- seq(
    min(dados_trait$water_percentage),
    max(dados_trait$water_percentage),
    length.out = 100
  )
  
  pred_mean <- intercept + slope * wp_seq
  pred_upper <- intercept + (slope + slope_sd) * wp_seq
  pred_lower <- intercept + (slope - slope_sd) * wp_seq
  
  linhas_modelo[[var]] <- data.frame(
    trait = var,
    water_percentage = wp_seq,
    herb_effect = pred_mean,
    upper = pred_upper,
    lower = pred_lower,
    signif = sig
  )
}

linhas_modelo <- bind_rows(linhas_modelo)

### identificar estrutura

dados_plot <- dados_plot %>%
  mutate(
    estrutura = case_when(
      grepl("^leaf_", trait) ~ "Leaf",
      grepl("^flower_", trait) ~ "Flower",
      grepl("^fruit_", trait) ~ "Fruit",
      grepl("^seed_", trait) ~ "Seed"
    )
  )

linhas_modelo <- linhas_modelo %>%
  mutate(
    estrutura = case_when(
      grepl("^leaf_", trait) ~ "Leaf",
      grepl("^flower_", trait) ~ "Flower",
      grepl("^fruit_", trait) ~ "Fruit",
      grepl("^seed_", trait) ~ "Seed"
    )
  )



########################################
### FIGURA POR ESTRUTURA
########################################

fig_paineis <- ggplot() +
  
  geom_point(
    data = dados_plot,
    aes(x = water_percentage, y = herb_effect),
    alpha = 1,
    cex = 1.5,
    color = "black"
  ) +
  
  geom_line(
    data = linhas_modelo,
    aes(
      x = water_percentage,
      y = herb_effect,
      color = estrutura,
      linetype = signif,
      group = trait
    ),
    linewidth = 1
  ) +
  
  facet_wrap(~estrutura, scales = "free") +
  
  scale_color_manual(values = c(
    "Leaf" = "#45cf09",
    "Flower" = "#fbc0f3",
    "Fruit" = "#ff7979",
    "Seed" = "#ffb160"
  )) +
  
  scale_linetype_manual(values = c(
    "significant" = "solid",
    "ns" = "dashed"
  )) +
  
  theme_classic() +
  
  labs(
    x = "Water percentage",
    y = "Herbarium effect",
    color = "Structure",
    linetype = "Significance"
  )

fig_paineis



#######################################################
###herb effect plot
#######################################################


library(dplyr)
library(tidyr)
library(ggplot2)

# Identifica colunas com "_herb_effect" e remove weigth/petiole
cols_herb_effect <- grep("_herb_effect$", names(df_media_spp), value = TRUE)
cols_herb_effect <- cols_herb_effect[!grepl("weigth|petiole", cols_herb_effect)]

# Transforma os dados para formato longo
df_long <- df_media_spp %>%
  dplyr::select(all_of(cols_herb_effect)) %>%
  pivot_longer(cols = everything(), names_to = "variavel", values_to = "valor") %>%
  mutate(
    # Remove "_herb_effect" do nome das variáveis
    variavel = gsub("_herb_effect$", "", variavel)
  )

# Ordem manual das variáveis
ordem <- rev(c(
  "leaf_area",
  "leaf_perimeter",
  "leaf_length",
  "leaf_width",
  "leaf_length_per_width",
  "leaf_circularity_index",
  "leaf_lma_area",
  "flower_petal_size",
  "flower_length",
  "flower_diameter_petal",
  "flower_diameter_stamen",
  "flower_hipant_length",
  "flower_hipant_diameter",
  "fruit_length",
  "fruit_diameter",
  "seed_length",
  "seed_diameter"
))
df_long$variavel <- factor(df_long$variavel, levels = ordem)

# Cria o gráfico
grafico_herb_effect <- ggplot(df_long, aes(x = variavel, y = valor)) +
  geom_boxplot(fill = "gray",  color = "black", linewidth = 0.3, outlier.color = "black", outlier.shape = 16, outlier.size = 0.8) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.8) +
  labs(x = "Trait", y = "Herbarium effect") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(linewidth = 0.7),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(size = 13, color = "black")
  ) +
  coord_flip()

# Exibe o gráfico
grafico_herb_effect


#looking in to the fruit special case
library(dplyr)
library(tidyr)
library(ggplot2)

# Select fruit traits for different effects
cols_fruit_effects <- c(
  "fruit_length_press_effect",
  "fruit_diameter_press_effect",
  "fruit_length_drying_effect",
  "fruit_diameter_drying_effect",
  "fruit_length_herb_effect",
  "fruit_diameter_herb_effect"
)

# Keep only columns present in the dataset
cols_fruit_effects <- cols_fruit_effects[cols_fruit_effects %in% names(df_media_voucher)]

# Transform to long format
df_long_fruit <- df_media_voucher %>%
  dplyr::select(all_of(cols_fruit_effects)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "trait",
    values_to = "value"
  )

# Optional: define order of traits for plotting
ordem_traits <- rev(c(
  "fruit_length_press_effect",
  "fruit_diameter_press_effect",
  "fruit_length_drying_effect",
  "fruit_diameter_drying_effect",
  "fruit_length_herb_effect",
  "fruit_diameter_herb_effect"
))
df_long_fruit$trait <- factor(df_long_fruit$trait, levels = ordem_traits)

# Create boxplot
grafico_fruit_effect <- ggplot(df_long_fruit, aes(x = trait, y = value)) +
  geom_boxplot(
    fill = "gray",      # box fill color
    color = "black",      # box border color
    linewidth = 0.8,         # box border thickness
    outlier.shape = 16,      # solid circle for outliers
    outlier.size = 2,        # size of outliers
    outlier.color = "black"    # color of outliers
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  coord_flip() +
  theme_minimal() +
  labs(
    x = "Trait",
    y = "Effect size",
    title = "Fruit traits: Press, Drying, and Herbarium effects"
  ) +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

# Display plot
grafico_fruit_effect

