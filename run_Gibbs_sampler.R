### Blocked Gibbs sampler

Nsteps_total <- 50
N_cell_count_adjustment_attempts_per_step <- M*K*geographical_hierarchy[[1]]$N # suggested
N_warmup_noprior <- N_cell_count_adjustment_attempts_per_step # suggested
N_report <- floor(N_cell_count_adjustment_attempts_per_step/10)
ycount <- zcount <- 0

current_cell_count <- true_cell_count_initial
acceptance_rate_tracker <- matrix(0,nrow=nrow(tablebuilder_outputs[[1]]$classes),ncol=ncol(tablebuilder_outputs[[1]]$classes))

for (z in 1:(N_warmup_noprior+Nsteps_total*N_cell_count_adjustment_attempts_per_step)) {
  
  proposal <- list()
  proposal$location <- sample(1:geographical_hierarchy[[1]]$N,1)
  proposal$category <- sample(1:K,1)
  proposal$value <- sample(c(-1,1),1)
  log_likelihood_diff <- evaluate_log_likelihood_diff_across_geographies(current_cell_count,proposal)
  if (!(is.na(log_likelihood_diff))) {
    log_likelihood_diff <- log_likelihood_diff+evaluate_log_prior_diff(current_cell_count,proposal)
  }
  if (!(is.na(log_likelihood_diff)) & (log_likelihood_diff > log(runif(1)))) {
    current_cell_count[[1]]$classes[proposal$location,proposal$category] <- current_cell_count[[1]]$classes[proposal$location,proposal$category] + proposal$value
    current_cell_count[[1]]$totals[proposal$location] <- current_cell_count[[1]]$totals[proposal$location] + proposal$value
    for (j in 2:M) {
      current_cell_count[[j]]$classes[geographical_hierarchy[[j]]$aggregation_vector[proposal$location],proposal$category] <- current_cell_count[[j]]$classes[geographical_hierarchy[[j]]$aggregation_vector[proposal$location],proposal$category] + proposal$value
      current_cell_count[[j]]$totals[geographical_hierarchy[[j]]$aggregation_vector[proposal$location]] <- current_cell_count[[j]]$totals[geographical_hierarchy[[j]]$aggregation_vector[proposal$location]] + proposal$value
    }
    acceptance_rate_tracker[proposal$location,proposal$category] <- acceptance_rate_tracker[proposal$location,proposal$category]+1
  }
  
  if ((z%%N_report)==1) {
    cat(colMeans(acceptance_rate_tracker>1)," : ",colSums(current_cell_count[[1]]$classes)/colSums(tablebuilder_outputs[[length(tablebuilder_outputs)]]$classes),"\n")
  }
  
  if ((z%%N_cell_count_adjustment_attempts_per_step)==0) {
    cat("updating prior ...\n")
    ycount <- ycount + 1
    prior_update <- update_prior()
    log_type_expect <- prior_update$log_type_expect
    log_count_expect <- prior_update$log_count_expect
    acceptance_rate_tracker <- acceptance_rate_tracker*0
  }
  if (((z%%(N_cell_count_adjustment_attempts_per_step/12)))==0) {
    zcount <- zcount + 1
    save(current_cell_count,file=paste0("./outputs/current_cell_count",sprintf("%03i",zcount),".dat"))
  }
  
}
