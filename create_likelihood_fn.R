### Define likelihood function for model of TableBuilder perturbation process

error_sd <- 2.0
logdiffexp <- function(y,x) { # supposes x > y
  x+log(1-exp(y-x))
}
logsumexp <- function(y,x) {
  if (length(x)==1) {
    max(c(x,y))+log(1+exp(min(c(x,y))-max(c(x,y))))} else {
      mmax <- apply(cbind(x,y),1,max)
      mmin <- apply(cbind(x,y),1,min)
      mmax+log(1+exp(mmin-mmax))
    }
}

log_likelihood_diff_fn <- function(proposed_value,current_value,observed_value) {
  
  if (proposed_value==0 & observed_value>0) {
    proposed_likelihood <- NA
  } else if (proposed_value==0 & observed_value==0) {
    proposed_likelihood <- 0
  } else if (proposed_value>0 & observed_value==0) {
    proposed_likelihood <- pnorm(2.5,proposed_value,error_sd,log.p=TRUE)
  } else if (abs(proposed_value-observed_value) < 5 & abs(current_value-observed_value) < 5) {
    proposed_likelihood <- logdiffexp(pnorm(observed_value-0.5,proposed_value,error_sd,log.p = TRUE),pnorm(observed_value+0.5,proposed_value,error_sd,log.p = TRUE))
  } else {
    proposed_likelihood <- dnorm(observed_value,proposed_value,error_sd,log = TRUE)
  }
  
  if (current_value==0 & observed_value==0) {
    current_likelihood <- 0
  } else if (current_value>0 & observed_value==0) {
    current_likelihood <- pnorm(2.5,current_value,error_sd,log.p=TRUE)
  } else if (abs(proposed_value-observed_value) < 5 & abs(current_value-observed_value) < 5) {
    current_likelihood <- logdiffexp(pnorm(observed_value-0.5,current_value,error_sd,log.p = TRUE),pnorm(observed_value+0.5,current_value,error_sd,log.p = TRUE))
  } else {
    current_likelihood <- dnorm(observed_value,current_value,error_sd,log = TRUE)
  }
  
  return(proposed_likelihood - current_likelihood)
}

evaluate_log_likelihood_diff_across_geographies <- function(current_cell_count,proposal) {
  log_likelihood_diff <- 0
  if (current_cell_count[[1]]$classes[proposal$location,proposal$category]==0 & proposal$value==-1) {
    log_likelihood_diff <- NA
  } else {
    current_value <- current_cell_count[[1]]$classes[proposal$location,proposal$category]
    proposed_value <- current_value + proposal$value
    observed_value <- tablebuilder_outputs[[1]]$classes[proposal$location,proposal$category]
    log_likelihood_diff <- log_likelihood_diff + log_likelihood_diff_fn(proposed_value,current_value,observed_value)
    current_value <- current_cell_count[[1]]$totals[proposal$location]
    proposed_value <- current_value + proposal$value
    observed_value <- tablebuilder_outputs[[1]]$totals[proposal$location]
    log_likelihood_diff <- log_likelihood_diff + log_likelihood_diff_fn(proposed_value,current_value,observed_value)
    for (j in 2:M) {
      current_value <- current_cell_count[[j]]$classes[geographical_hierarchy[[j]]$aggregation_vector[proposal$location],proposal$category]
      proposed_value <- current_value + proposal$value
      observed_value <- tablebuilder_outputs[[j]]$classes[geographical_hierarchy[[j]]$aggregation_vector[proposal$location],proposal$category]
      log_likelihood_diff <- log_likelihood_diff + log_likelihood_diff_fn(proposed_value,current_value,observed_value)
      current_value <- current_cell_count[[j]]$totals[geographical_hierarchy[[j]]$aggregation_vector[proposal$location]]
      proposed_value <- current_value + proposal$value
      observed_value <- tablebuilder_outputs[[j]]$totals[geographical_hierarchy[[j]]$aggregation_vector[proposal$location]]
      log_likelihood_diff <- log_likelihood_diff + log_likelihood_diff_fn(proposed_value,current_value,observed_value)
    }
  }
  return(log_likelihood_diff)
}
