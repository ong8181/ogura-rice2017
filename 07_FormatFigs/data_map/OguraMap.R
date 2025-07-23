####
#### Ogura field map
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0
library(patchwork); packageVersion("patchwork") # 1.3.0
library(ggmap); packageVersion("ggmap") # 4.0.0


# ------------------------------------------ #
# Visualize sampling locations
# ------------------------------------------ #
# Register AIP
ggmap::register_stadiamaps("xxxxxxxxxxxxxxxxxxxxxxxx")
# Draw Hong Kong Map
ogura_lonlat1 <- c(left = 135.76, right = 135.80, bottom = 34.88, top = 34.92)
ogura_lonlat2 <- c(left = 135.4, right = 136.52, bottom = 34.7, top = 35.6)
jp_lonlat <- c(left = 127, right = 148, bottom = 29, top = 47)
# Maptype "alidade_smooth", "outdoors", "stamen_terrain"
ogura_map1 <- get_stadiamap(ogura_lonlat, zoom = 14, maptype = "stamen_terrain")
ogura_map2 <- get_stadiamap(ogura_lonlat2, zoom = 10, maptype = "stamen_terrain")
jp_map <- get_stadiamap(jp_lonlat, zoom = 6, maptype = "stamen_terrain")
ggmap(ogura_map1)
ggmap(ogura_map2)
ggmap(jp_map)

# Detailed map
ogura_lonlat3 <- c(left = 135.773, right = 135.777, bottom = 34.901, top = 34.9035)
ogura_map3 <- get_stadiamap(ogura_lonlat3, zoom = 18, maptype = "outdoors")
g3 <- ggmap(ogura_map3) + theme_bw()

#ogura_lonlat4 <- c(long = 135.7745, lat = 34.90225)
#ogura_map4 <- get_googlemap(center = ogura_lonlat4, zoom = 10, maptype = "satellite")


# Sampling locations
g1 <- ggmap(jp_map) +
  xlab(expression(paste("Longitude (", degree, "E)"))) +
  ylab(expression(paste("Latitude (", degree, "N)")))
g2 <- ggmap(ogura_map2) +
  xlab(expression(paste("Longitude (", degree, "E)"))) +
  ylab(expression(paste("Latitude (", degree, "N)")))


# Save plot
#ggsave("Ogura_map.jpg", (g1 + g2), width = 10, height = 8)
#ggsave("Ogura_map.pdf", (g1 + g2), width = 10, height = 8)
#ggsave("Ogura_Survey_map.pdf", g3, width = 10, height = 8)

