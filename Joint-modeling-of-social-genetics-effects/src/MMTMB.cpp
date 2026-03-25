
#include <TMB.hpp>
#include "MMTMB.hpp"
/**
   @file MMTMB.hpp
   \brief Global function to handle different model types in TMB. This function selects and runs the appropriate model based on the input string.
 
   The function supports. See the include file MMTMB.hpp for details
   - `MMTMB`: Multivariate  Model
 
   @author Jemay Salomon, PhD-Student
   @date 2023-2026
*/

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(model);
  if (model == "MMTMB") {
    return MMTMB(this);
  } else if(model=="lmm"){
    return lmm(this);
  }else{Rf_error("Unknown model type"); }
  return 0;
}


