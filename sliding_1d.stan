data {
  int<lower=1> N;          // the number of observations.
  int<lower=1> tau;        // sliding window size
  int<lower=1> t1;          // the number of observations.
  int<lower=1> max_ww;     // sw_row
  int SW[max_ww,t1];      // sliding window matrix
  int<lower=0> Y[N];       // observed cases.
  int<lower=1> S;          // length of serial interval
  vector[S] W;             // serial interval
}

parameters {
  // time and region specific R values, in log space
  real logR[max_ww];
}

transformed parameters {

  // expected value of cases, **J** each gets its own window, lots of zeros
  real<lower=0.00001> M[N, max_ww] = rep_array(0.1, N, max_ww);
  real<lower=0.00001> R[max_ww] = rep_array(1, max_ww);

  // ------ CALCULATE R(t) and INITIALIZE M[1] -------------
  // get R in exp() space for each window
  for(ww in 1:max_ww) {
    // First, get the R of this window
    R[ww] = exp(logR[ww]);
    //R[ww] = exp(xbeta[ww]);
    // initialize
    M[1, ww] = Y[1]; // **J**
  }

  // ------ CALCULATE M(t) -------------
  /// NOW FOR EACH BLOCK
  for(ww in 1:max_ww) {

      // next calculate the starting and ending n
      int startN = SW[ww, 1];
      int endN   = SW[ww, tau];

      // ok now for each of these get the Ms that occur only for these windows
      // the question is, what to do about
      for(m_i in startN:endN) {
        real inner_vec = 0;
        int sip_i = min(S, startN);
        for(si in 1:sip_i) {
          inner_vec = inner_vec + W[si] * M[(m_i + 1) - si, ww];
        }
        M[(m_i + 1), ww] = R[ww] * inner_vec;
      }

      // carry across
      for(wremain in min(ww + 1, max_ww):max_ww) {
        M[(startN + 1), wremain] = M[(startN + 1),ww];
      }


  }

}


model {

  // ------ SIGMA, BETA, and LogR ------
  // priors and sample
  //xsigma ~ inv_gamma(2, 1);
  logR ~ normal(0, 1);

  // WINDOW SPECIFIC LogR
  //for(ww in 1:max_ww) {
  //  logR[ww] ~ normal(xbeta[ww], xsigma);
  //}

  for(ww in 1:max_ww) {

      // Offset by 1 because its the end of the window
      int startN = SW[ww, 2];
      int endN   = SW[ww, t1];

      // SO THIS ENFORCES THAT THE OBSERVED CASES
      // ARE DRAWN FROM ALL OF THE WINDOWS THAT THIS DAY FALLS INTO
      // This is also what forces the betas of adjacent windows to be
      // related. and you have to add each target separately so its 1:1
      for (n in startN:endN) {
        Y[n] ~ poisson(M[n, ww]);
      }

  }

}


