library(tidyverse)

#######################################################
# Preparing data
#######################################################

# Loading and organizing dataset
trait_data<-read.table(file="data_trait.txt", head=T, sep="\t",dec = ("."))
df<- trait_data
df[, 8:70] <- lapply(df[, 8:70], as.numeric)
df<-df[,-c(1,8)]
names(df) <- gsub("^petal_", "flower_petal_", names(df))
names(df) <- gsub("^hipant_", "flower_hipant_", names(df))
names(df) <- gsub("^circularity_", "leaf_circularity_", names(df))
names(df) <- gsub("^lma_", "leaf_lma_", names(df))
#df <- df[, !names(df) %in% c(
# "leaf_lma_area_fresh",
# "leaf_lma_area_dry"
#)]
str(df)

# Identify columns with suffix _fresh and _dry
col_fresh <- names(df)[str_detect(names(df), "_fresh$")]
col_dry   <- names(df)[str_detect(names(df), "_dry$")]

# Extract base name (without suffix) from fresh and dry columns
base_fresh <- str_replace(col_fresh, "_fresh$", "")
base_dry   <- str_replace(col_dry, "_dry$", "")

# Find matching pairs present in both
paired_vars <- intersect(base_fresh, base_dry)

# For each pair, create a new column with relative difference
for (var in paired_vars) {
  col_fresh_name <- paste0(var, "_fresh")
  col_dry_name   <- paste0(var, "_dry")
  new_col_name   <- paste0(var, "_herb_effect")
  
  df[[new_col_name]] <- (df[[col_dry_name]] / df[[col_fresh_name]]) - 1}

# saving checkpoint
write.csv(df,file = "full_data.csv",row.names = F)

# mean per voucher (individual)
df_media_voucher <- df %>%
  group_by(voucher) %>%
  summarise(
    spp = first(spp),
    section = first(section),
    genus = first(genus),
    superfilo_name = first(superfilo_name),
    # Add more categorical columns you want to keep here
    across(
      where(is.numeric), 
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop")

# merging climate data with voucher mean data
clim_data <- read.csv(file="climate_data.csv")
clim_data <- data.frame(
  voucher = clim_data$voucher,
  temp_annual=clim_data$bio1,
  temp_sazon=clim_data$bio4,
  precip_annual=clim_data$bio12,
  precip_sazon=clim_data$bio15)
df_media_voucher <- merge.data.frame(df_media_voucher,clim_data,by="voucher")

# mean per species
df_media_spp <- df_media_voucher %>%
  group_by(spp) %>%
  summarise(
    genus = first(genus),
    section = first(section),
    superfilo_name = first(superfilo_name),
    # Add more categorical columns you want to keep here
    across(
      where(is.numeric), 
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop")

# saving checkpoint
write.csv(df_media_voucher,file = "voucher_mean_data.csv",row.names = F)
write.csv(df_media_spp,file = "spp_mean_data.csv",row.names = F)