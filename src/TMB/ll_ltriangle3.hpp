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

// Compute the negative log-likelihood of the log-asymmetric-triangular
// distribution. If Y ~ ltriangle3(locationlog, scalelog, skew) then log(Y) has
// a triangular distribution with mode locationlog and support
//   [locationlog - (1 - skew) * scalelog, locationlog + (1 + skew) * scalelog]
// so the total width is 2 * scalelog and skew = 0 recovers the symmetric
// ltriangle distribution. The density on the log scale is, with a and b the
// support limits and c the mode,
//   2 (y - a) / ((b - a) (c - a))  for a <= y <= c
//   2 (b - y) / ((b - a) (b - c))  for c <  y <= b
// The lower and upper limbs have widths wl = (1 - skew) scalelog and
// wu = (1 + skew) scalelog, so in terms of the dimensionless height
//   t = (y - a) / wl  below the mode, t = (b - y) / wu above it,
// the density is t / scalelog on both limbs.
//
// See ll_ltriangle.hpp for the reasoning behind the soft-log barrier: the
// density is only piecewise smooth and is exactly zero outside the support, so
// log() is replaced by a C1 linear extension below a small threshold to keep
// the objective finite and make it unprofitable to exclude an observation.
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
//    skew         - skewness in [-1, 1], bounded by the optimizer rather than transformed

/// @file ll_ltriangle3.hpp

#ifndef ll_ltriangle3_hpp
#define ll_ltriangle3_hpp

// Cumulative distribution function of the asymmetric triangular distribution
// with support [a, b] and mode c, evaluated at x. `wl` and `wu` are the
// (floored) widths of the lower and upper limbs, c - a and b - c.
template<class Type>
Type ptri_ltriangle3(Type x, Type a, Type b, Type c, Type wl, Type wu) {
  Type w = b - a;
  Type lower = (x - a) * (x - a) / (w * wl); // a < x <= c
  Type upper = Type(1.0) - (b - x) * (b - x) / (w * wu); // c < x < b
  Type mid = CppAD::CondExpLe(x, c, lower, upper);
  return CppAD::CondExpLe(
    x, a, Type(0.0),
    CppAD::CondExpGe(x, b, Type(1.0), mid));
}

// Softened logarithm: log(x) above `eps`, and the C1 linear extension of log()
// below it (see softlog_ltriangle in ll_ltriangle.hpp).
template<class Type>
Type softlog_ltriangle3(Type x, Type eps) {
  Type floored = CppAD::CondExpGt(x, eps, x, eps);
  return CppAD::CondExpGt(x, eps, log(floored), log(eps) + (x - eps) / eps);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type ll_ltriangle3(objective_function<Type>* obj) {
  // Data
  DATA_VECTOR( left  );  // left and right values
  DATA_VECTOR( right );
  DATA_VECTOR( weight);  // weight

  // The order of these parameter statements determines the order of the estimates in the vector of parameters
  PARAMETER( locationlog );
  PARAMETER( log_scalelog );
  PARAMETER( skew );

  Type scalelog;
  scalelog = exp(log_scalelog);  // convert to [0,Inf] scale

  // `scalelog` is unbounded below, so guard against it underflowing to zero.
  // The limb widths reach exactly zero at skew = +/-1 (a right triangle), so
  // floor them as well to keep every branch of the CondExp expressions finite.
  Type tiny = Type(1e-300);
  Type scale = CppAD::CondExpGt(scalelog, tiny, scalelog, tiny);
  Type wl = (Type(1.0) - skew) * scale;  // width of the lower limb, c - a
  Type wu = (Type(1.0) + skew) * scale;  // width of the upper limb, b - c
  wl = CppAD::CondExpGt(wl, tiny, wl, tiny);
  wu = CppAD::CondExpGt(wu, tiny, wu, tiny);
  Type a = locationlog - wl;  // lower support limit
  Type b = locationlog + wu;  // upper support limit

  Type eps_dens = Type(1e-6);    // softlog threshold for the density
  Type eps_mass = Type(1e-10);   // softlog threshold for censored interval mass
  Type tmin = Type(-1e9);        // clamp keeping the density extension finite

  Type nll = 0;  // negative log-likelihood
  int n_data = left.size(); // number of data values
  Type pleft;    // probability that concentration < left(i)  used for censored data
  Type pright;   // probability that concentration < right(i) used for censored data

  for( int i=0; i<n_data; i++){
     if(left(i) == right(i)){   // uncensored values
        // dimensionless height of the density at log(y) on whichever limb the
        // observation falls, negative outside the support so that the soft
        // logarithm acts as a scale-free barrier
        Type y = log(left(i));
        Type dev = y - locationlog;
        Type tlower = Type(1.0) + dev / wl;
        Type tupper = Type(1.0) - dev / wu;
        Type t = CppAD::CondExpLe(dev, Type(0.0), tlower, tupper);
        t = CppAD::CondExpGt(t, tmin, t, tmin);
        nll -= weight(i) *
          (softlog_ltriangle3(t, eps_dens) - log(scale) - y);
     };
     if(left(i) < right(i)){    // censored values
        // `tc` is the smaller of the dimensionless distances of the censoring
        // interval's limits inside the support, so it is positive iff the
        // interval overlaps the support and is capped at one when it contains
        // the mode.
        Type tc = Type(1.0);
        pleft = 0;
        if(left(i)>0){
           Type logleft = log(left(i));
           pleft = ptri_ltriangle3(logleft, a, b, locationlog, wl, wu);
           Type tleft = (b - logleft) / wu;
           tc = CppAD::CondExpLt(tleft, tc, tleft, tc);
        };
        pright = 1;
        using std::isfinite;
        if(isfinite(right(i))){
           Type logright = log(right(i));
           pright = ptri_ltriangle3(logright, a, b, locationlog, wl, wu);
           Type tright = (logright - a) / wl;
           tc = CppAD::CondExpLt(tright, tc, tright, tc);
        };
        // The interval mass is exactly zero, with zero gradient, whenever the
        // interval lies outside the candidate support, so soften log() and add
        // the same linear barrier as the uncensored branch (see
        // ll_ltriangle.hpp).
        tc = CppAD::CondExpGt(tc, tmin, tc, tmin);
        Type barrier = CppAD::CondExpLt(tc, Type(0.0), tc / eps_dens, Type(0.0));
        nll -= weight(i) *
          (softlog_ltriangle3(pright - pleft, eps_mass) + barrier);
     };

  };

  ADREPORT(scalelog);
  REPORT  (scalelog);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
