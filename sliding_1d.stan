data {
  int<lower=1> N;          // the number of observations.
  int<lower=1> tau;        // sliding window size
  int<lower=1> max_ww;     // sw_row
  int SW[max_ww,tau];      // sliding window matrix
  int<lower=0> Y[N];       // observed cases.
  int<lower=1> S;          // length of serial interval
  vector[S] W;             // serial interval
  real init_cases;          // initial cases
}

parameters {
  // time and region specific R values, in log space
  real logR[max_ww];
}

transformed parameters {

  // expected value of cases, **J** each gets its own window, lots of zeros
  real<lower=0> M[N, max_ww] = rep_array(0.1, N, max_ww);
  real<lower=0> R[max_ww] = rep_array(1, max_ww);

  // ------ CALCULATE R(t) and INITIALIZE M[1] -------------
  // get R in exp() space for each window
  for(ww in 1:max_ww) {
    R[ww] = exp(logR[ww]);
  }

  // ------ CALCULATE M(t) -------------
  /// NOW FOR EACH BLOCK
  for(ww in 1:max_ww) {

      // next calculate the starting and ending n
      int startN = SW[ww, 1];
      int endN   = SW[ww, tau];

      // ok now for each of these get the Ms that occur only for these windows
      // the question is, what to do about
      for(t in startN:endN) {

        real inner_vec = 0;
        int S_loop_max = min(S, t + 1);
        int  forward_vec[S_loop_max];
        int  rev_vec[S_loop_max];
        for(si in 1:S_loop_max) {
          forward_vec[si] = si;
          rev_vec[si] = t + 1 - si;
        }
        real mx = 0;

        for(si in 1:S_loop_max) {
          int rev_i = rev_vec[si];

          if(rev_i < 0) {
            mx = 0;
          }

          if(rev_i == 0) {
            mx = init_cases;
          }

          if(rev_i > 0){
            mx = M[rev_i, ww];
          }

          inner_vec = inner_vec + W[si] * mx;
        }
        M[t, ww] = R[ww] * inner_vec;
      }

      // carry across
      if(startN < N) {
        for(wremain in min(ww + 1, max_ww):max_ww) {
          M[startN , wremain] = M[startN, ww];
        }
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
      int startN = SW[ww, 1];
      int endN   = SW[ww, tau];

      // SO THIS ENFORCES THAT THE OBSERVED CASES
      // ARE DRAWN FROM ALL OF THE WINDOWS THAT THIS DAY FALLS INTO
      // This is also what forces the betas of adjacent windows to be
      // related. and you have to add each target separately so its 1:1
      for (n in startN:endN) {
        Y[n] ~ poisson(M[n, ww]);
      }

  }

}


