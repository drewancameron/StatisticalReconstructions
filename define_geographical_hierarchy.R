### Define geographical hierarchy

geographical_hierarchy <- list() # describes the relationship between unit levels

## <<< Perth example
library(rgdal)
MB_polygons <- readOGR("MB_PerthGCC.shp")
SA1_polygons <- readOGR("SA1_PerthGCC.shp")
SA2_polygons <- readOGR("SA2_PerthGCC.shp")
lowest_level <- list()
lowest_level$N <- length(MB_polygons)
geographical_hierarchy[[1]] <- lowest_level
next_level <- list()
next_level$N <- length(SA1_polygons)
next_level$aggregation_vector <- match(MB_polygons$SA1_CODE21,SA1_polygons$SA1_CODE21)
geographical_hierarchy[[2]] <- next_level
next_level <- list()
next_level$N <- length(SA2_polygons)
next_level$aggregation_vector <- match(MB_polygons$SA2_CODE21,SA2_polygons$SA2_CODE21)
geographical_hierarchy[[3]] <- next_level
## Perth example >>>

M <- length(geographical_hierarchy) # number of levels of the geographical hierarchy
