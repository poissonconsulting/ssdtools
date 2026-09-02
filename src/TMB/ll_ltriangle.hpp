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
// The log-likelihood is only piecewise smooth: the |log(y) - locationlog| term
// puts a kink at every observation. This is inherent to the triangular family
// rather than an artefact of this implementation, and the kink set has measure
// zero in the parameter space, so L-BFGS-B copes, but a future reader should
// not assume the smoothness guarantees of the other distributions.
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

/// @file ll_ltriangle.hpp

#ifndef ll_ltriangle_hpp
#define ll_ltriangle_hpp

// Cumulative distribution function of the symmetric triangular distribution
// with mode `location` and half-width `scale`, evaluated at x.
template<class Type>
Type ptri_ltriangle(Type x, Type location, Type scale) {
  Type z = (x - location) / scale;
  Type lower = (z + Type(1.0)) * (z + Type(1.0)) / Type(2.0); // -1 < z <= 0
  Type upper = Type(1.0) - (Type(1.0) - z) * (Type(1.0) - z) / Type(2.0); // 0 < z < 1
  Type mid = CppAD::CondExpLe(z, Type(0.0), lower, upper);
  return CppAD::CondExpLe(
    z, Type(-1.0), Type(0.0),
    CppAD::CondExpGe(z, Type(1.0), Type(1.0), mid));
}

// Softened logarithm: log(x) above `eps`, and the C1 linear extension of log()
// below it. Continuous with a continuous first derivative at `eps`, finite for
// every finite argument, and decreasing without bound as the argument falls, so
// it acts as a barrier rather than as a floor. `x` must already be clamped from
// below so that the linear extension stays finite.
template<class Type>
Type softlog_ltriangle(Type x, Type eps) {
  Type floored = CppAD::CondExpGt(x, eps, x, eps);
  return CppAD::CondExpGt(x, eps, log(floored), log(eps) + (x - eps) / eps);
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

  // `scalelog` is unbounded below, so guard against it underflowing to zero,
  // which would make both the dimensionless deviation and log(scalelog)
  // non-finite.
  Type tiny = Type(1e-300);
  Type scale = CppAD::CondExpGt(scalelog, tiny, scalelog, tiny);

  Type eps_dens = Type(1e-6);    // softlog threshold for the density
  Type eps_mass = Type(1e-10);   // softlog threshold for censored interval mass
  Type tmin = Type(-1e9);        // clamp keeping the density extension finite

  Type nll = 0;  // negative log-likelihood
  int n_data = left.size(); // number of data values
  Type pleft;    // probability that concentration < left(i)  used for censored data
  Type pright;   // probability that concentration < right(i) used for censored data

  for( int i=0; i<n_data; i++){
     if(left(i) == right(i)){   // uncensored values
        // pdf of the symmetric triangular on the log scale is
        //   (scalelog - |log(y) - locationlog|) / scalelog^2
        // transformed to the concentration scale via the 1/y Jacobian
        // (- log(left)). Writing this in terms of the dimensionless
        //   t = 1 - |log(y) - locationlog| / scalelog
        // gives log(t) - log(scalelog) - log(y), so replacing log() by
        // `softlog_ltriangle` makes the out-of-support penalty scale free: it
        // grows without bound as the candidate support shrinks away from a data
        // point, rather than flattening at a constant, which would make
        // excluding the point profitable.
        Type dev = log(left(i)) - locationlog;
        Type absdev = CppAD::CondExpGe(dev, Type(0.0), dev, -dev);
        Type t = Type(1.0) - absdev / scale;
        t = CppAD::CondExpGt(t, tmin, t, tmin);
        nll -= weight(i) *
          (softlog_ltriangle(t, eps_dens) - log(scale) - log(left(i)));
     };
     if(left(i) < right(i)){    // censored values
        // `dist` is the distance on the log scale from the mode to the nearest
        // point of the censoring interval, zero when the mode lies inside it.
        Type dist = Type(0.0);
        pleft = 0;
        if(left(i)>0){
           Type logleft = log(left(i));
           pleft = ptri_ltriangle(logleft, locationlog, scale);
           Type dleft = logleft - locationlog;
           dist = CppAD::CondExpGt(dleft, dist, dleft, dist);
        };
        pright = 1;
        using std::isfinite;
        if(isfinite(right(i))){
           Type logright = log(right(i));
           pright = ptri_ltriangle(logright, locationlog, scale);
           Type dright = locationlog - logright;
           dist = CppAD::CondExpGt(dright, dist, dright, dist);
        };
        // Unlike the unbounded distributions, `ptri_ltriangle` returns exactly
        // 0 and exactly 1 outside the support, so the interval mass reaches
        // exactly 0, with zero gradient, whenever a censoring interval lies
        // entirely outside the candidate support. Soften log() so the
        // objective stays finite, and add the same scale-free linear barrier
        // as the uncensored branch on the interval's distance outside the
        // support, so that excluding a censored observation is never
        // profitable. `tc` is positive iff the interval overlaps the support.
        Type tc = Type(1.0) - dist / scale;
        tc = CppAD::CondExpGt(tc, tmin, tc, tmin);
        Type barrier = CppAD::CondExpLt(tc, Type(0.0), tc / eps_dens, Type(0.0));
        nll -= weight(i) *
          (softlog_ltriangle(pright - pleft, eps_mass) + barrier);
     };

  };

  ADREPORT(scalelog);
  REPORT  (scalelog);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
