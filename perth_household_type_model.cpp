#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  // Input Data
  
  DATA_INTEGER(N_MB);
  DATA_STRUCT(spde_perth,spde_t); // INLA SPDE object (components of precision matrix)
  DATA_SPARSE_MATRIX(A_perth); // INLA SPDE projection matrix: N_mesh to N_MB
  
  DATA_IVECTOR(irsads);
  
  DATA_INTEGER(N_MB_type);
  DATA_IVECTOR(MB_types);
  
  DATA_VECTOR(N_pos); 
  DATA_VECTOR(N_neg); 
  
  // Parameters
  
  PARAMETER(log_sd_field);
  PARAMETER(log_range_field);
  PARAMETER(log_sd_iid_residential);
  PARAMETER(log_sd_iid_other);
  
  PARAMETER(intercept);
  
  PARAMETER_VECTOR(irsad_offsets); //11
  
  PARAMETER_VECTOR(type_offsets); //N_MB_type
  
  PARAMETER_VECTOR(field);
  PARAMETER_VECTOR(iid_errors);
  
  // Parameter Transforms
  
  Type range = exp(log_range_field);
  Type kappa = 2.8284/range;
  Type sd = exp(log_sd_field);
  Type sd_iid_residential = exp(log_sd_iid_residential);
  Type sd_iid_other = exp(log_sd_iid_other);
  
  // Priors
  
  Type nll = 0.0;
  
  SparseMatrix<Type> Q_perth = Q_spde(spde_perth,kappa);
  nll += SCALE(GMRF(Q_perth),sd)(field);
  
  nll -= dnorm(log_sd_field,Type(-1),Type(0.5),true);
  nll -= dnorm(log_range_field,Type(1),Type(0.5),true);
  
  nll -= dnorm(log_sd_iid_residential,Type(-1),Type(0.5),true);
  nll -= dnorm(log_sd_iid_other,Type(-1),Type(0.5),true);
  
  nll -= dnorm(irsad_offsets,Type(0),Type(1.0),true).sum();
  nll -= dnorm(type_offsets,Type(0),Type(1.0),true).sum();
  
  for (int i=0; i<N_MB; i++) {
    if (MB_types[i]==1) {
      nll -= dnorm(iid_errors[i],Type(0),sd_iid_residential,true);
      } else {
      nll -= dnorm(iid_errors[i],Type(0),sd_iid_other,true);
    }
  }
  
  // Algebra
  
  vector<Type> baseline_field_perth(N_MB);
  baseline_field_perth = A_perth*field;
  
  vector<Type> irsad_effects(N_MB);
  for (int i=0; i<N_MB; i++) {irsad_effects[i] = irsad_offsets[irsads[i]-1];}
  vector<Type> type_effects(N_MB);
  for (int i=0; i<N_MB; i++) {type_effects[i] = type_offsets[MB_types[i]-1];}
  
  vector<Type> predicted_surface_hcfmd_perth(N_MB);
  predicted_surface_hcfmd_perth = intercept  + baseline_field_perth.array() + irsad_effects.array() + type_effects.array() + iid_errors.array();
  
  predicted_surface_hcfmd_perth = invlogit(predicted_surface_hcfmd_perth).array();
  
  // Likelihood
  
  for (int i=0; i<N_MB; i++) {
    nll -= dbinom(N_pos[i],N_pos[i]+N_neg[i],predicted_surface_hcfmd_perth[i],true);
  }
  
  // Reporting

  REPORT(irsad_effects);
  REPORT(type_effects);
  REPORT(field);
  
  REPORT(predicted_surface_hcfmd_perth);
  
  return nll;
}
