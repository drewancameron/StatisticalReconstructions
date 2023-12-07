build_mock_dataset <- function() {
load("./outputs/current_cell_count100.dat") 
fiducial <- current_cell_count
perturbed <- fiducial
for (i in 1:3) {
  nonzero <- which(perturbed[[i]]$classes > 0,arr.ind=TRUE)
  for (j in 1:length(nonzero[,1])) {
    req_mean <- perturbed[[i]]$classes[nonzero[j,1],nonzero[j,2]]
    dummy <- round(rnorm(1,mean=req_mean,sd=2.0))
    perturbed[[i]]$classes[nonzero[j,1],nonzero[j,2]] <- dummy*as.integer(dummy>2)
  }
  nonzero <- which(perturbed[[i]]$totals > 0)
  for (j in 1:length(nonzero)) {
    req_mean <- perturbed[[i]]$totals[nonzero[j]]
    dummy <- round(rnorm(1,mean=req_mean,sd=2.0))
    perturbed[[i]]$totals[nonzero[j]] <- dummy*as.integer(dummy>2)
  }
}
return(list('perturbed'=perturbed,'fiducial'=fiducial))
}
