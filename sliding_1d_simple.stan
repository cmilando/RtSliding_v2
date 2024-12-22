data {
  int<lower=1> N;          // the number of observations.
  int<lower=0> Y[N];       // observed cases.
  int<lower=1> S;          // length of serial interval
  vector[S] W;             // serial interval
}

parameters {
  // time and region specific R values, in log space
  real logR[N-1];
  real xbeta[N-1];
  real<lower=0> xsigma;
}

transformed parameters {

  // expected value of cases, **J** each gets its own window, lots of zeros
  real<lower=0> M[N] = rep_array(0.001, N);
  real<lower=0> R[N-1] = rep_array(0.1, N-1);

  // ------ CALCULATE R(t) and INITIALIZE M[1] -------------
  // get R in exp() space for each window
  for(w in 1:(N-1)) {
    R[w] = exp(logR[w]);
    //

  }

  // Initialize
  M[1] = Y[1] * 1.0;

  // ------ CALCULATE M(t) -------------
  /// NOW FOR EACH BLOCK
  for(t in 2:N) {

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
        mx = M[rev_i];
      }

      inner_vec = inner_vec + W[si] * mx;
    }

    ///
    int wi = t - 1;
    M[t] = R[wi] * inner_vec;
  }

}


// THIS IS GOOD DON'T CHANGE IT
model {

  // ------ SIGMA, BETA, and LogR ------
  // priors and sample
  xsigma ~ inv_gamma(2, 1);
  xbeta ~ normal(0, 1);

  for(w in 1:(N-1)) {
    // WINDOW SPECIFIC LogR
    logR[w] ~ normal(xbeta[w], xsigma);
  }

  // target
  for(n in 1:N) {
    target += poisson_lpmf(Y[n] | M[n]);
  }

}

