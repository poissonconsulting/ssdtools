// Copyright 2015-2023 Province of British Columbia
// Copyright 2021 Environment and Climate Change Canada
// Copyright 2023-2024 Australian Government Department of Climate Change, 
// Energy, the Environment and Water
//
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
//
//       https://www.apache.org/licenses/LICENSE-2.0
//
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

// Compute the negative log-likelihood of the log-logistic distribution
// If Y ~ log-logistic( locationlog, scalelog) then log(Y) ~ logisitic(locationlog, scalelog)
// 
// Input data are left(1...n) right(1...n) weight(1...n)
// where 
//    n = sample size (inferred from the vectors)
//    left(i) right(i) specify the uncensored or censored data as noted below
//    weight(i)  - relative weight to be given to each observation's log-likelihood. Use values of 1 for ordinary likelihood
//
//  left(i) and right(i) can take the following forms
//     left(i) == right(i)  - non-censored data
//     left(i) <  right(i)  - interval censored data
//  left(i) must be non-negative (all concentrations must be non-negative)
//  right(i) can take the value Inf for no upper limit
// 
//  E.g.  left(i) right(i)
//          3       3       non-censored values
//          0       3       0 < concentration < 3
//          3       Inf     3 < concentration.
//
// Parameters are
//    locationlog - location on the log(Concentration) scale
//    log_scalelog    - log(scalelog)    on the log(Concentration) scale, i.e. scalelog=exp(log_scalelog)

// Refer to http://kaskr.github.io/adcomp/matrix_arrays_8cpp-example.html for help in coding the log-likelihood function

// Refer to https://github.com/kaskr/adcomp/wiki/Development
// on instructions for including TMB code in an R package

/// @file ll_llogis.hpp

#ifndef ll_llogis_hpp
#define ll_llogis_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type ll_llogis(objective_function<Type>* obj) {
  // Data
  DATA_VECTOR( left  );  // left and right values
  DATA_VECTOR( right );
  DATA_VECTOR( weight);  // weight

  // The order of these parameter statements determines the order of the estimates in the vector of parameters
  PARAMETER( locationlog );
  PARAMETER( log_scalelog    );

  Type scalelog;
  scalelog = exp(log_scalelog);  // convert to [0,Inf] scale
  
  Type nll = 0;  // negative log-likelihood
  int n_data    = left.size(); // number of data values
  Type pleft;    // probability that concentration < left(i)  used for censored data
  Type pright;   // probability that concentration < right(i) used for censored data
 
  // Probability of data conditional on parameter values for uncensored data
  // pdf of log(normal) obtained from pdf(normal) using the standard transformation theory
  for( int i=0; i<n_data; i++){
     if(left(i) == right(i)){   // uncensored values
        nll -= weight(i)*(dlogis( log(left(i)), locationlog, scalelog, true ) - log(left(i)));   // log likelihood for uncensored values
     };
     if(left(i) < right(i)){    // censored values; no builtin function so we code the cdf directly
        pleft = 0;
        if(left(i)>0){ pleft=1/(1+exp(-(log(left(i))-locationlog)/scalelog));};
        pright = 1;
        using std::isfinite;
        if(isfinite(right(i))){ pright=1/(1+exp(-(log(right(i))-locationlog)/scalelog));};
        nll -= weight(i)*log(pright-pleft);
     };
     
  };

  ADREPORT(scalelog);
  REPORT  (scalelog);
  
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
