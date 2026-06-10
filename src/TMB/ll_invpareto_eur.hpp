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

// Compute the negative log-likelihood of the European inverse Pareto distribution
// If Y ~ pareto (shape, scale) then 1/Y is inverse pareto (shape, scale)

// We are using the European Pareto distribution (the version implemented in the
// actuar package). Unlike the North American version (see ll_invpareto.hpp) it
// is unbounded above and so the scale parameter is estimated as a free parameter.

// The European inverse Pareto distribution has
//    cdf F(y) = (y / (y + scale))^shape, for y, shape, scale > 0
//    pdf f(y) = shape * scale * y^(shape - 1) / (y + scale)^(shape + 1)

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
// Parameters are
//    log_shape  - log(shape)
//    log_scale  - log(scale)

#ifndef ll_invpareto_eur_hpp
#define ll_invpareto_eur_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type ll_invpareto_eur(objective_function<Type>* obj)
{
  // Data
  DATA_VECTOR( left  );  // left and right values
  DATA_VECTOR( right );
  DATA_VECTOR( weight);  // weight

  // The order of these parameter statements determines the order of the estimates in the vector of parameters
   // Parameters
  PARAMETER( log_shape );
  PARAMETER( log_scale );

  Type shape = exp(log_shape);
  Type scale = exp(log_scale);

  Type nll = 0;
  int n_data    = left.size(); // number of data values
  Type pleft;    // probability that concentration < left(i)  used for censored data
  Type pright;   // probability that concentration < right(i) used for censored data
  Type y;

  for( int i=0; i<n_data; i++){
     if(left(i) == right(i)){  // uncensored data
        y = left(i);
        nll -= weight(i)*(log(shape) + log(scale) + (shape - Type(1))*log(y) - (shape + Type(1))*log(y + scale));
     };
     if(left(i) < right(i)){   // censored data
        pleft = 0;
        if(left(i) > 0){ pleft = pow(left(i)/(left(i)+scale), shape); };  // F(left)
        pright = 1;
        using std::isfinite;
        if(isfinite(right(i))){ pright = pow(right(i)/(right(i)+scale), shape); };  // F(right)
        nll -= weight(i)*log(pright-pleft);  // contribution to log-likelihood for censored values
     };

  };

  ADREPORT(shape);
  REPORT  (shape);
  ADREPORT(scale);
  REPORT  (scale);

  return nll;
};

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
