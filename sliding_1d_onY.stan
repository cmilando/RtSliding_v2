data {
  int<lower=1> N;          // the number of observations.
  int<lower=1> tau;        // sliding window size
  int<lower=1> max_ww;     // sw_row
  int SW[max_ww,2];        // sliding window matrix
  int<lower=0> Y[N];       // observed cases.
  int<lower=1> S;          // length of serial interval
  vector[S] W;             // serial interval
}

parameters {
  // time and region specific R values, in log space
  real logR[max_ww];
  real xbeta[max_ww];
  real<lower=0> xsigma;
}

transformed parameters {

  // expected value of cases, **J** each gets its own window, lots of zeros
  real<lower=0> M[N] = rep_array(0.1, N);
  real<lower=0> R[max_ww] = rep_array(1, max_ww);

  // ------ CALCULATE R(t) and INITIALIZE M[1] -------------
  // get R in exp() space for each window
  for(ww in 1:max_ww) {
    R[ww] = exp(logR[ww]);
    //
    M[1] = Y[1] * 1.0;
  }

  // ------ CALCULATE M(t) -------------
  /// NOW FOR EACH BLOCK
  for(ww in 1:max_ww) {

      // next calculate the starting and ending n
      int startN = SW[ww, 1];
      int endN   = SW[ww, 2];  // THIS IS THE R(t) value you are setting

      // ok now for each of these get the Ms that occur only for these windows
      for(t in startN:endN) {

        // SETUP
        real inner_vec = 0;
        real mx = 0;
        int S_loop_max = min(S, t + 1);
        int  forward_vec[S_loop_max];
        int  rev_vec[S_loop_max];
        for(si in 1:S_loop_max) {
          forward_vec[si] = si;
          rev_vec[si] = t + 1 - si;
        }

        // LAMBDA from Cori 2013
        // CWM added a -1 so that you never get to W[1] = 0
        for(si in 1:(S_loop_max - 1)) {
          int rev_i = rev_vec[si];
          if(rev_i < 1) {
            mx = 0;
          }
          if(rev_i >= 1){
            mx = M[rev_i];
          }
          inner_vec = inner_vec + W[si] * mx;
        }

        // Finally, set (AND OVERWRITE) M(t)
        // Over writing is how you ensure they are connected to eachother
        M[t] = R[ww] * inner_vec;
        // so for a sliding window of 2

      }

  }

}


model {

  // ------ SIGMA, BETA, and LogR ------
  // priors and sample
  xsigma ~ inv_gamma(2, 1);
  xbeta ~ normal(0, 1);

  // WINDOW SPECIFIC LogR
  for(ww in 1:max_ww) {
    logR[ww] ~ normal(xbeta[ww], xsigma);
  }

  // SO THIS ENFORCES THAT THE OBSERVED CASES
  // ARE DRAWN FROM ALL OF THE WINDOWS THAT THIS DAY FALLS INTO
  // This is also what forces the betas of adjacent windows to be
  // related. and you have to add each target separately so its 1:1
  // AND YOU DONT HAVE TO DO all ww because they *should* be all the same
  for (n in 1:N) {
    Y[n] ~ poisson(M[n]);
  }

}


