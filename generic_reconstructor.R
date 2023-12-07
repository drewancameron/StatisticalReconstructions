### Generic code for reconstructions of Census cell counts from nested geographical TableBuilder outputs 
### Data structures are made concrete with a Perth GCC example

setwd("~/VIRTUAL_WA/STATISTICAL_RECONSTRUCTIONS/") # delete from final


## Analysis on actual TableBuilder outputs: TableBuilder files not provided in repo

source("define_geographical_hierarchy.R")

source("read_tablebuilder_outputs.R")

source("mcmc_init.R")

source("create_likelihood_fn.R")

source("define_priors.R")

source("run_Gibbs_sampler.R")

source("summarise_results.R")

results_list <- summarise_results(61,610)
map_results(results_list)

## Analysis on mock dataset with minimal perturbation model

source("define_geographical_hierarchy.R")

source("build_mock_dataset_minimal.R")

mockdataset <- build_mock_dataset()
save(mockdataset,file="mock_minimal.dat")
tablebuilder_outputs <- mockdataset$perturbed

source("mcmc_init.R")

source("create_likelihood_fn.R")

source("define_priors.R")

source("run_Gibbs_sampler.R")

source("summarise_results.R")

results_list <- summarise_results(61,610)
coverage <- summarise_results_mock_data(results_list,fiducial)
map_results(results_list)
add_coverage_MB(coverage[[1]])

## Analysis on mock dataset with alternative perturbation model

source("define_geographical_hierarchy.R")

source("build_mock_dataset_alternative.R")

mockdataset <- build_mock_dataset()
save(mockdataset,file="mock_alternative.dat")

tablebuilder_outputs <- mockdataset$perturbed

source("mcmc_init.R")

source("create_likelihood_fn.R")

source("define_priors.R")

source("run_Gibbs_sampler.R")

source("summarise_results.R")

results_list <- summarise_results(61,610)
coverage <- summarise_results_mock_data(results_list,fiducial)
map_results(results_list)
add_coverage_MB(coverage[[1]])

