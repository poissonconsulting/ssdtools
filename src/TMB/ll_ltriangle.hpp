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

// Compute the negative log-likelihood of the log-triangular distribution
// If Y ~ log-triangular(locationlog, scalelog) then log(Y) ~ symmetric
// triangular with mode locationlog and half-width scalelog, i.e. support
// [locationlog - scalelog, locationlog + scalelog].
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
// Parameters are
//    locationlog  - mode on the log(Concentration) scale
//    log_scalelog - log(scalelog) on the log(Concentration) scale, i.e. scalelog=exp(log_scalelog)

/// @file ll_triangle.hpp

#ifndef ll_ltriangle_hpp
#define ll_ltriangle_hpp

// Cumulative distribution function of the symmetric triangular distribution
// with mode `location` and half-width `scale`, evaluated at x.
template<class Type>
Type ptri(Type x, Type location, Type scale) {
  Type z = (x - location) / scale;
  Type lower = (z + Type(1.0)) * (z + Type(1.0)) / Type(2.0); // -1 < z <= 0
  Type upper = Type(1.0) - (Type(1.0) - z) * (Type(1.0) - z) / Type(2.0); // 0 < z < 1
  Type mid = CppAD::CondExpLe(z, Type(0.0), lower, upper);
  return CppAD::CondExpLe(
    z, Type(-1.0), Type(0.0),
    CppAD::CondExpGe(z, Type(1.0), Type(1.0), mid));
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type ll_ltriangle(objective_function<Type>* obj) {
  // Data
  DATA_VECTOR( left  );  // left and right values
  DATA_VECTOR( right );
  DATA_VECTOR( weight);  // weight

  // The order of these parameter statements determines the order of the estimates in the vector of parameters
  PARAMETER( locationlog );
  PARAMETER( log_scalelog );

  Type scalelog;
  scalelog = exp(log_scalelog);  // convert to [0,Inf] scale

  Type nll = 0;  // negative log-likelihood
  int n_data = left.size(); // number of data values
  Type pleft;    // probability that concentration < left(i)  used for censored data
  Type pright;   // probability that concentration < right(i) used for censored data

  for( int i=0; i<n_data; i++){
     if(left(i) == right(i)){   // uncensored values
        // pdf of the symmetric triangular on the log scale is
        //   (scalelog - |log(y) - locationlog|) / scalelog^2
        // transformed to the concentration scale via the 1/y Jacobian (- log(left)).
        Type dev = log(left(i)) - locationlog;
        Type absdev = CppAD::CondExpGe(dev, Type(0.0), dev, -dev);
        Type inside = scalelog - absdev;
        // floor the density argument to keep the likelihood finite (and act as a
        // soft barrier) if a candidate support fails to cover a data point. The
        // floor is scaled to `scalelog` (the peak density argument) so the soft
        // barrier behaves consistently regardless of the magnitude of the data.
        Type floor = Type(1e-8) * scalelog;
        Type floored = CppAD::CondExpGt(inside, floor, inside, floor);
        nll -= weight(i) * (log(floored) - Type(2.0) * log(scalelog) - log(left(i)));
     };
     if(left(i) < right(i)){    // censored values
        pleft = 0;
        if(left(i)>0){ pleft = ptri(log(left(i)), locationlog, scalelog); };
        pright = 1;
        using std::isfinite;
        if(isfinite(right(i))){ pright = ptri(log(right(i)), locationlog, scalelog); };
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
