library(raster)
library(geodata)
library(dplyr)
library(sp)
library(progress)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# 2. Read and convert coordinates
coord <- read.table(
  file   = "data_coord.txt",
  header = TRUE,
  sep    = "\t",
  dec    = ".")

# Latitude in decimal degrees
coord$lat_dec <- with(
  coord,
  (lat.g + lat.m / 60 + lat.s / 3600) *
    ifelse(n.s == "S", -1, 1))

# Longitude in decimal degrees
coord$lon_dec <- with(
  coord,
  (long.g + long.m / 60 + long.s / 3600) *
    ifelse(w.e == "W", -1, 1))

# Final data frame of points
df <- data.frame(
  voucher = coord$voucher,
  lat     = coord$lat_dec,
  lon     = coord$lon_dec)

# 3. Download and prepare environmental layers
# WorldClim v2.1 (19 bioclimatic variables)
bio <- worldclim_global(
  var     = "bio",
  res     = 2.5,
  version = "2.1",
  path    = getwd())

bio_stack <- stack(bio)
names(bio_stack) <- paste0("bio", seq_len(nlayers(bio_stack)))
bio_layers <- unstack(bio_stack)


# Final list of layers
layers <- c(bio_layers)
names(layers) <- c(paste0("bio", 1:19))

# 4. Extract environmental values
points_sp <- df[, c("lon", "lat")]
coordinates(points_sp) <- ~ lon + lat
proj4string(points_sp) <- CRS("+proj=longlat +datum=WGS84")

env_data <- as.data.frame(points_sp)

pb <- progress_bar$new(
  format = "Extracting layers [:bar] :percent elapsed: :elapsed",
  total  = length(layers),
  width  = 60)

for (i in seq_along(layers)) {
  env_data[[names(layers)[i]]] <-
    raster::extract(layers[[i]], points_sp)
  pb$tick()
}

climate_data <- bind_cols(
  voucher = df$voucher,
  env_data)

write.csv(
  climate_data,
  "climate_data.csv",
  row.names = FALSE)

# 5. Map of sampling points
# 5. Aggregate vouchers by coordinates
res <- 0.3   # cell size in degrees (e.g., 0.5°, 0.25°, 0.1°)

df_map <- df %>%
  mutate(
    lon_bin = floor(lon / res) * res,
    lat_bin = floor(lat / res) * res
  ) %>%
  group_by(lon_bin, lat_bin) %>%
  summarise(
    n_vouchers = n(),
    lon = mean(lon),
    lat = mean(lat),
    .groups = "drop"
  )

set.seed(123)  # reproducibility

df_map_jitter <- df_map %>%
  mutate(
    lon_j = lon + runif(n(), -0.5, 0.5),
    lat_j = lat + runif(n(), -0.5, 0.5)
  )

# Convert to sf
pontos_sf <- st_as_sf(
  df_map_jitter,
  coords = c("lon_j", "lat_j"),
  crs    = 4326)

# Base world map
mundo <- ne_countries(
  scale = "medium",
  returnclass = "sf")

# Bounding box with margin
bbox <- st_bbox(pontos_sf)
xlim <- c(bbox$xmin - 25, bbox$xmax + 5)
ylim <- c(bbox$ymin - 10, bbox$ymax + 10)

# Plot
mapa_pontos <- ggplot() +
  geom_sf(data = mundo, fill = "gray95", color = "gray70") +
  geom_sf(
    data = pontos_sf,
    aes(fill = n_vouchers),
    shape  = 21,
    size   = 3.5,          # fixed size
    color  = "black",
    stroke = 0.3,
    alpha  = 0.8
  ) +
  scale_fill_gradientn(
    name = "Number of sampled specimens",
    colours = c(
      "#fde725",
      "orange",
      "darkorange",
      "red","darkred"
    ),
    limits = c(1, max(df_map$n_vouchers))
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude"
  )

# Map with variable point size
mapa_pontos_size <- ggplot() +
  geom_sf(data = mundo, fill = "gray95", color = "gray70") +
  geom_sf(
    data = pontos_sf,
    aes(fill = n_vouchers, size = n_vouchers),
    shape  = 21,
    color  = "black",
    stroke = 0.3,
    alpha  = 0.8
  ) +
  scale_fill_gradientn(
    name = "Number of sampled specimens",
    colours = c(
      "#fde725",
      "orange",
      "darkorange",
      "red", "darkred"
    ),
    limits = c(1, max(df_map$n_vouchers))
  ) +
  scale_size_continuous(
    name   = "Number of sampled specimens",
    range  = c(2, 8),
    limits = c(1, max(df_map$n_vouchers))
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude"
  )

# Display map
print(mapa_pontos)
print(mapa_pontos_size)

# Save map as SVG
ggsave(
  filename = "sampling_points_map.svg",
  plot     = mapa_pontos,
  device   = "svg",
  width    = 8,
  height   = 6,
  units    = "in"
)

