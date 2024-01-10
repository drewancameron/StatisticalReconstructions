### Define priors (and read in any ancillary files to define those: e.g. spatial covariates)

evaluate_log_prior_diff <- function(current_cell_count,proposal) { # dummy function definition
  log_prob_current <- 0
  log_prob_proposed <- 0
  log_prior_difference <- log_prob_proposed-log_prob_current
  return(log_prior_difference)
}

update_prior <- function(current_cell_count) {
}

## Perth example <<<
filename_SA1_covariate <- "SA1_IRSADs.csv"
xdata <- read.csv(filename_SA1_covariate,skip=11,header = FALSE)
SA1_COVARIATES <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("SA1_COVARIATES <- cbind(SA1_COVARIATES,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(filename_SA1_covariate,skip=9,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA1",ydata[-1][1:(length(SA1_COVARIATES[1,])-1)]))
colnames(SA1_COVARIATES) <- xcolnames
SA1_COVARIATES <- as.data.frame(SA1_COVARIATES)
SA1_COVARIATES <- SA1_COVARIATES[-which(is.na(SA1_COVARIATES$SA1)),]
SA1_COVARIATES <- SA1_COVARIATES[match(as.numeric(unique(SA1_polygons$SA1_CODE21)),SA1_COVARIATES$SA1),]
SA1_COVARIATES <- as.matrix(SA1_COVARIATES)
SA1_COVARIATES <- SA1_COVARIATES[,-c(1,12)]
SA1_COVARIATES[SA1_COVARIATES>0] <- 1
SA1_COVARIATES <- apply(SA1_COVARIATES,1,which.max)*as.integer(rowSums(SA1_COVARIATES)>0)
SA1_COVARIATES[SA1_COVARIATES==0] <- 11
library(INLA)
coords_MB <- coordinates(MB_polygons)
perth.mesh <-  inla.mesh.2d(coords_MB,cutoff = 0.01,max.n=300)
perth_A <- inla.mesh.project(perth.mesh,coords_MB)$A
perth_spde <- (inla.spde2.matern(perth.mesh,alpha=2)$param.inla)[c("M0","M1","M2")]
IRSAD_MB <- SA1_COVARIATES[geographical_hierarchy[[2]]$aggregation_vector]
MB_type <- match(MB_polygons$MB_CAT21,unique(MB_polygons$MB_CAT21))
library(TMB)
compile("perth_household_type_model.cpp")
dyn.load(dynlib("perth_household_type_model"))
compile("perth_household_count_model.cpp")
dyn.load(dynlib("perth_household_count_model"))
update_prior <- function() {
  input.data <- list('N_MB'=geographical_hierarchy[[1]]$N,
                     'spde_perth'=perth_spde,
                     'A_perth'=perth_A,
                     'irsads'=IRSAD_MB,
                     'N_MB_type'=max(MB_type),
                     'MB_types'=MB_type,
                     'N_pos'=current_cell_count[[1]]$classes[,2],
                     'N_neg'=current_cell_count[[1]]$totals-current_cell_count[[1]]$classes[,2])
  if (sum(as.integer((input.data$N_neg+input.data$N_pos)==0)) > 0) {
    input.data$N_neg[(input.data$N_neg+input.data$N_pos)==0] <- 1}
  parameters <- list('log_sd_field'=0,
                     'log_range_field'=0,
                     'log_sd_iid_residential'=0,
                     'log_sd_iid_other'=0,
                     'intercept'=0,
                     'irsad_offsets'=rep(-2,11),
                     'type_offsets'=rep(0,input.data$N_MB_type),
                     'field'=rep(0,perth.mesh$n),
                     'iid_errors'=rep(0,geographical_hierarchy[[1]]$N))
  obj <- MakeADFun(input.data,parameters,DLL='perth_household_type_model',random=c('irsad_offsets','type_offsets','field','iid_errors'))
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  rep <- sdreport(obj, getJointPrecision = TRUE)
  input.data <- list('N_MB'=geographical_hierarchy[[1]]$N,
                     'spde_perth'=perth_spde,
                     'A_perth'=perth_A,
                     'irsads'=IRSAD_MB,
                     'N_MB_type'=max(MB_type),
                     'MB_types'=MB_type,
                     'N_tot'=rowSums(current_cell_count[[1]]$classes))
  parameters <- list('log_sd_field'=0,
                     'log_range_field'=0,
                     'log_sd_iid_residential'=0,
                     'log_sd_iid_other'=0,
                     'intercept'=0,
                     'irsad_offsets'=rep(-2,11),
                     'type_offsets'=rep(0,input.data$N_MB_type),
                     'field'=rep(0,perth.mesh$n),
                     'iid_errors'=rep(0,geographical_hierarchy[[1]]$N))
  objx <- MakeADFun(input.data,parameters,DLL='perth_household_count_model',random=c('irsad_offsets','type_offsets','field','iid_errors'))
  optx <- nlminb(objx$par, objx$fn, objx$gr)
  repx <- sdreport(objx, getJointPrecision = TRUE)
  library(sparseMVN)
  xsample <- rmvn.sparse(10, c(repx$par.fixed, repx$par.random), Cholesky(repx$jointPrecision), prec = TRUE)
  colnames(xsample) <- names(c(repx$par.fixed,repx$par.random))
  save(xsample,file=paste0("./outputs/count_regression_perth_",sprintf("%03i",ycount),".dat"))
  ysample <- rmvn.sparse(10, c(rep$par.fixed, rep$par.random), Cholesky(rep$jointPrecision), prec = TRUE)
  colnames(ysample) <- names(c(rep$par.fixed,rep$par.random))
  save(ysample,file=paste0("./outputs/type_regression_perth_",sprintf("%03i",ycount),".dat"))
  xsample <- xsample[1,]
  ysample <- ysample[1,]
  log_count_expect <- log(objx$report(xsample)$predicted_surface_hcfmd_perth)
  log_type_expect <- log(obj$report(ysample)$predicted_surface_hcfmd_perth)
  cat("quantiles: Types: ",quantile(exp(log_type_expect),c(0.001,0.05,0.25,0.5,0.75,0.95,0.999)),"\n")
  cat("quantiles: Counts: ",quantile(exp(log_count_expect),c(0.001,0.05,0.25,0.5,0.75,0.95,0.999)),"\n")
  return(list('log_count_expect'=log_count_expect,'log_type_expect'=log_type_expect))
}
log_count_expect <- log(rep(5,geographical_hierarchy[[1]]$N))
log_type_expect <- log(rep(0.05,geographical_hierarchy[[1]]$N))
evaluate_log_prior_diff <- function(current_cell_count,proposal) { # dummy function definition
  log_prob_current <- 0
  log_prob_proposed <- 0
  log_prob_current <- log_prob_current + dbinom(current_cell_count[[1]]$classes[proposal$location,2],current_cell_count[[1]]$totals[proposal$location],exp(log_type_expect[proposal$location]),log=TRUE)
  log_prob_proposed <- log_prob_proposed + dbinom(current_cell_count[[1]]$classes[proposal$location,2]+as.integer(proposal$category==2)*proposal$value,current_cell_count[[1]]$totals[proposal$location]+proposal$value,exp(log_type_expect[proposal$location]),log=TRUE)
  log_prob_current <- log_prob_current + dpois(current_cell_count[[1]]$totals[proposal$location],exp(log_count_expect[proposal$location]),log=TRUE)
  log_prob_proposed <- log_prob_proposed + dpois(current_cell_count[[1]]$totals[proposal$location]+proposal$value,exp(log_count_expect[proposal$location]),log=TRUE)
  log_prior_difference <- log_prob_proposed-log_prob_current
  return(log_prior_difference)
}
## Perth example >>>
