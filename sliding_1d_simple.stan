data {
  int<lower=1> N_windows;       // the number of windows.
  int<lower=1> N_obs;           // the number of observations
  int<lower=0> Y[N_obs];        // observed cases in windows.
  int<lower=1> S;               // length of serial interval
  vector[S] W;                  // serial interval
}

parameters {
  // time and region specific R values, in log space
  real logR[N_windows  - 1];
  real xbeta[N_windows - 1];
  real<lower=0> xsigma;
}

transformed parameters {

  // expected value of cases in each WINDOW
  real<lower=0> M_window[N_windows] = rep_array(0.001, N_windows);

  // expected value of cases in each TIMESTEP
  real<lower=0> M[N_obs] = rep_array(0.001, N_obs);

  // no matter what N is, R will always be 1 less
  // than the number of Y
  real<lower=0> R[N_windows - 1] = rep_array(0.1, N_windows - 1);

  // ------ CALCULATE R(t) and INITIALIZE M[1] -------------
  // get R in exp() space for each window
  for(w in 1:(N_windows - 1)) {
    R[w] = exp(logR[w]);
  }

  // Initialize
  M_window[1] = Y[1] * 1.0;

  // ------ CALCULATE M(t) -------------
  /// NOW FOR EACH BLOCK
  for(t in 2:N_windows) {

    ///
    real inner_vec = 0;
    int S_loop_max = min(S, t + 1);
    int  forward_vec[S_loop_max];
    int  rev_vec[S_loop_max];

    ///
    for(si in 1:S_loop_max) {
      forward_vec[si] = si;
      rev_vec[si] = t + 1 - si;
    }

    ///
    real mx = 0;

    ///
    for(si in 1:S_loop_max) {
      int rev_i = rev_vec[si];

      if(rev_i < 1) {
        mx = 0;
      }

      if(rev_i >= 1){
        mx = M_window[rev_i];
      }

      inner_vec = inner_vec + W[si] * mx;
    }

    ///
    int wi = t - 1;
    M_window[t] = R[wi] * inner_vec;
  }

  // ------ REVERSE OUT M(t) -------------
  for(n in 1:N_obs) {
    M[n] = M_window[n];
  }

}


// Keep this simple
model {

  // ------ SIGMA, BETA, and LogR ------
  // priors and sample
  xsigma ~ inv_gamma(2, 1);
  xbeta ~ normal(0, 1);

  for(w in 1:(N_windows - 1)) {
    // WINDOW SPECIFIC LogR
    logR[w] ~ normal(xbeta[w], xsigma);
  }

  // target
  for(n in 1:N_obs) {
    target += poisson_lpmf(Y[n] | M[n]);
  }

}

