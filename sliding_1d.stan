data {
  int<lower=1> N;          // the number of observations.
  int<lower=1> tau;        // sliding window size
  int<lower=1> max_ww;     // the number of sliding windows
  int SW[max_ww,2];        // sliding window matrix
  int<lower=0> Y[N];       // observed cases.
  int<lower=1> S;          // length of serial interval
  vector[S] W;             // serial interval
}

parameters {
  real logR[max_ww];       // time and region specific R values, in log space
  real xbeta[max_ww];
  real<lower=0> xsigma;
  real logInitCases;       // a guess at the initial cases
}

transformed parameters {

  // expected value of cases, **J** each gets its own window, lots of zeros
  real<lower=0> M[N, max_ww] = rep_array(0.1, N, max_ww);
  real<lower=0> R[max_ww]    = rep_array(1, max_ww);

  // ------ CALCULATE R(t) and INITIALIZE M[1] -------------
  // get R in exp() space for each window
  for(ww in 1:max_ww) {
    R[ww] = exp(logR[ww]);
    ///////////////////////////////////////
    // THIS IS WHERE THINGS NEED TO CHANGE for back imputations
    M[1, ww] = Y[1] * 1.0;
    ///////////////////////////////////////
  }

  // ------ CALCULATE M(t) -------------
  /// NOW FOR EACH BLOCK
  for(ww in 1:max_ww) {

      // the starting and ending n
      int startN = SW[ww, 1];
      int endN   = SW[ww, 2];

      // calculate Ms that occur only for these windows
      for(t in startN:endN) {

        // Setting up the forward and backwards iterators
        real inner_vec = 0;
        int  S_loop_max = min(S, t + 1); // This S and t + 1 because ...
        int  forward_vec[S_loop_max];
        int  rev_vec[S_loop_max];
        for(si in 1:S_loop_max) {
          forward_vec[si] = si;
          rev_vec[si] = t + 1 - si;
        }

        real mx = 0;

        // LAMBDA from Cori 2013.
        for(si in 1:S_loop_max) {

          int rev_i = rev_vec[si];

          // CASE 1:
          if(rev_i < 0) {
            mx = 0.;
          }
          ///////////////////////////////////////
          // CASE 2: THIS IS WHERE THINGS NEED TO CHANGE for back imputations
          if(rev_i == 0) {
            mx = exp(logInitCases);
          }
          ///////////////////////////////////////

          // CASE 3: past M values do exist
          if(rev_i > 0){
            mx = M[rev_i, ww];
          }

          // Sum up the inner vector
          // remember that W[1] will ALWAYS = 0
          inner_vec = inner_vec + W[si] * mx;
        }
        M[t, ww] = R[ww] * inner_vec;
      }

      // ------ CARRY FORWARDS -------------
      // carry your current guess forwards in time
      // this is poor man's recursion
      // you have to do this because otherwise when you do
      // M[rev_i, ww] in the lambda loop above, the next ww
      // won't have any values in it
      if(startN < N) {
        for(wremain in min(ww + 1, max_ww):max_ww) {
          M[startN , wremain] = M[startN, ww];
        }
      }
   }
}


model {

  // ------ SIGMA, BETA, and LogR and Guess ------
  //
  logInitCases ~ normal(0, 1);

  // priors and sample
  xsigma ~ inv_gamma(2, 1);  // this gets the variance across the region
  xbeta ~ normal(0, 1);      // this gets the value

  // WINDOW SPECIFIC LogR
  for(ww in 1:max_ww) {
    logR[ww] ~ normal(xbeta[ww], xsigma);
  }

  // ------ TARGET ------
  for(ww in 1:max_ww) {

      // Get sliding window
      int startN = SW[ww, 1];
      int endN   = SW[ww, 2];

      // SO THIS ENFORCES THAT THE OBSERVED CASES
      // ARE DRAWN FROM ALL OF THE WINDOWS THAT THIS DAY FALLS INTO
      // This is also what forces the betas of adjacent windows to be
      // related. and you have to add each target separately so its 1:1
      // Wow, alarming that that ran without needing mw
      for (n in startN:endN) {
        target += poisson_lpmf(Y[n] | M[n, ww]);
      }

  }

}


