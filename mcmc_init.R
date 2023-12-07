### Create starting state for MCMC sampler: forced to meet nested geography constraints

true_cell_count_initial <- list()
true_cell_count_initial[[1]] <- list()
true_cell_count_initial[[1]]$classes <- tablebuilder_outputs[[1]]$classes
true_cell_count_initial[[1]]$totals <- rowSums(true_cell_count_initial[[1]]$classes)
check_rowSums <- which(true_cell_count_initial[[1]]$totals==0 & tablebuilder_outputs[[1]]$totals>0)
if (length(check_rowSums)>0) {
  for (i in 1:length(check_rowSums)) {
    true_cell_count_initial[[1]]$classes[check_rowSums[i],sample(1:K,1)] <- 1
  }
}
true_cell_count_initial[[1]]$totals <- rowSums(true_cell_count_initial[[1]]$classes)

for (j in 2:M) {
  true_cell_count_initial[[j]] <- list()
  true_cell_count_initial[[j]]$classes <- aggregate(true_cell_count_initial[[1]]$classes,list(geographical_hierarchy[[j]]$aggregation_vector),sum)[,-1]
  check_aggregates <- which(true_cell_count_initial[[j]]$classes==0 & tablebuilder_outputs[[j]]$classes>0,arr.ind = TRUE)
  if (length(check_aggregates)>0) {
    for (i in 1:(dim(check_aggregates)[1])) {
      true_cell_count_initial[[1]]$classes[sample(which(geographical_hierarchy[[j]]$aggregation_vector==check_aggregates[i,1]),1),check_aggregates[i,2]] <- 1
    }
  }
  true_cell_count_initial[[j]]$totals <- aggregate(rowSums(true_cell_count_initial[[1]]$classes),list(geographical_hierarchy[[j]]$aggregation_vector),sum)[,-1]
  check_rowSums <- which(true_cell_count_initial[[j]]$totals==0 & tablebuilder_outputs[[j]]$totals>0)
  if (length(check_rowSums)>0) {
    for (i in 1:length(check_rowSums)) {
      true_cell_count_initial[[1]]$classes[sample(which(geographical_hierarchy[[j]]$aggregation_vector==check_rowSums[i]),1),sample(1:K,1)] <- 1
    }
  }
}

true_cell_count_initial[[1]]$totals <- rowSums(true_cell_count_initial[[1]]$classes)
for (j in 2:M) {
  true_cell_count_initial[[j]]$classes <- aggregate(true_cell_count_initial[[1]]$classes,list(geographical_hierarchy[[j]]$aggregation_vector),sum)[,-1]
  true_cell_count_initial[[j]]$totals <- aggregate(rowSums(true_cell_count_initial[[1]]$classes),list(geographical_hierarchy[[j]]$aggregation_vector),sum)[,-1]
}
