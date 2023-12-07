summarise_results <- function(start_i,end_i) {
  
  ran <- sample(1:length(start_i:end_i),2,replace=FALSE)
  
  results_list <- list()
  results_list$count_summaries <- list()
  for (k in 1:length(geographical_hierarchy)) {
    results_list$count_summaries[[k]] <- list()
    cell_count_array <- list()
    for (i in start_i:end_i) {
      load(paste0("./outputs/current_cell_count",sprintf("%03i",i),".dat"))
      cell_count_array[[length(cell_count_array)+1]] <- current_cell_count[[k]]$classes[,2]
    }
    cell_count_array <- do.call(rbind,cell_count_array)
    cell_count_array_mean <- apply(cell_count_array,2,mean)
    cell_count_array_median <- apply(cell_count_array,2,median)
    upperq <- function(x){quantile(x,0.975)}
    lowerq <- function(x){quantile(x,0.025)}
    cell_count_array_lower <- apply(cell_count_array,2,lowerq)
    cell_count_array_upper <- apply(cell_count_array,2,upperq)
    results_list$count_summaries[[k]]$mean <- cell_count_array_mean
    results_list$count_summaries[[k]]$median <- cell_count_array_median
    results_list$count_summaries[[k]]$upper <- cell_count_array_upper
    results_list$count_summaries[[k]]$lower <- cell_count_array_lower
    results_list$count_summaries[[k]]$draw_one <- cell_count_array[ran[1],]
    results_list$count_summaries[[k]]$draw_two <- cell_count_array[ran[2],]
  }
  save(results_list,file="./outputs/results_summarised.dat")
  return(results_list)
}

map_results <- function(results_list) {
  
  library(rgeos)
  library(rgdal)
  MB_polygons$MEA <- results_list$count_summaries[[1]]$mean
  MB_polygons$MED <- results_list$count_summaries[[1]]$median
  MB_polygons$HIG <- results_list$count_summaries[[1]]$upper
  MB_polygons$LOW <- results_list$count_summaries[[1]]$lower
  MB_polygons$RANO <- results_list$count_summaries[[1]]$draw_one
  MB_polygons$RANT <- results_list$count_summaries[[1]]$draw_two
  
  SA1_polygons$MEA <- results_list$count_summaries[[2]]$mean
  SA1_polygons$MED <- results_list$count_summaries[[2]]$median
  SA1_polygons$HIG <- results_list$count_summaries[[2]]$upper
  SA1_polygons$LOW <- results_list$count_summaries[[2]]$lower
  SA1_polygons$RANO <- results_list$count_summaries[[2]]$draw_one
  SA1_polygons$RANT <- results_list$count_summaries[[2]]$draw_two
  
  SA2_polygons$MEA <- results_list$count_summaries[[3]]$mean
  SA2_polygons$MED <- results_list$count_summaries[[3]]$median
  SA2_polygons$HIG <- results_list$count_summaries[[3]]$upper
  SA2_polygons$LOW <- results_list$count_summaries[[3]]$lower
  SA2_polygons$RANO <- results_list$count_summaries[[3]]$draw_one
  SA2_polygons$RANT <- results_list$count_summaries[[3]]$draw_two
  
  writeOGR(MB_polygons,".","MB_output",driver="ESRI Shapefile",overwrite=TRUE)
  writeOGR(SA1_polygons,".","SA1_output",driver="ESRI Shapefile",overwrite=TRUE)
  writeOGR(SA2_polygons,".","SA2_output",driver="ESRI Shapefile",overwrite=TRUE)
}

summarise_results_mock_data <- function(results_list,fiducial) {
  coverage_MB <- as.integer(fiducial[[1]]$classes[,2] <= results_list$count_summaries[[1]]$upper & 
                              fiducial[[1]]$classes[,2] >= results_list$count_summaries[[1]]$lower)
  coverage_SA1 <- as.integer(fiducial[[2]]$classes[,2] <= results_list$count_summaries[[2]]$upper & 
                               fiducial[[2]]$classes[,2] >= results_list$count_summaries[[2]]$lower)
  coverage_SA2 <- as.integer(fiducial[[3]]$classes[,2] <= results_list$count_summaries[[3]]$upper & 
                               fiducial[[3]]$classes[,2] >= results_list$count_summaries[[3]]$lower)
  return(list('coverage_MB'=coverage_MB,'coverage_SA1'=coverage_SA1,'coverage_SA2'=coverage_SA2))
}

add_coverage_MB <- function(coverage_MB) {
  MB_polygons <- readOGR("MB_output.shp")
  MB_polygons$COV <- coverage_MB
  writeOGR(MB_polygons,".","MB_output",driver="ESRI Shapefile",overwrite=TRUE)
}
