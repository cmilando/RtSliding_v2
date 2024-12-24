data {
  int<lower=1> N;          // the number of observations.
  int<lower=1> J;          // the number of regions
  real P[J, J];            // transfer matrix for P
  int<lower=1> tau;        // sliding window size
  int<lower=1> max_ww;     // the number of sliding windows
  int SW[max_ww,2];        // sliding window matrix
  int<lower=0> Y[N, J];    // observed cases.
  int<lower=1> S;          // length of serial interval
  vector[S] W;             // serial interval
}

parameters {
  real logR[max_ww, J];       // time and region specific R values, in log space
  real xbeta[max_ww, J];
  real<lower=0> xsigma[J];
  real logInitCases[J];       // a guess at the initial cases
}

transformed parameters {

  // expected value of cases, **J** each gets its own window, lots of zeros
  real<lower=0> M[N, max_ww, J] = rep_array(0.1, N, max_ww, J);
  real<lower=0> R[max_ww, J]    = rep_array(1, max_ww, J);

  // ------ CALCULATE R(t) and INITIALIZE M[1] -------------
  // get R in exp() space for each window
  for(j in 1:J) {
    for(ww in 1:max_ww) {
      R[ww, j] = exp(logR[ww, j]);
      M[1, ww, j] = Y[1, j] * 1.0;
    }
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
        int  S_loop_max = min(S, t + 1);
        // This S and t + 1 so you get back to M[1] when t=2, because w[1] = 0 always
        int  forward_vec[S_loop_max];
        int  rev_vec[S_loop_max];
        for(si in 1:S_loop_max) {
          forward_vec[si] = si;
          rev_vec[si] = t + 1 - si;
        }

        // LAMBDA from Cori 2013, modified to include J
        for(j in 1:J) {

          // these get reset for each j
          real inner_vec = 0;
          real mx = 0;

          // Loop through
          for(si in 1:S_loop_max) {

            int rev_i = rev_vec[si];

            ///////////////////////////////////////
            // BACK IMPUTATIONS -- just one step for now
            if(rev_i < 0) {
              mx = 0.;
            }
            if(rev_i == 0) {
              mx = exp(logInitCases[j]);
            }
            ///////////////////////////////////////

            // CASE 2: past M values do exist
            if(rev_i > 0){
              mx = M[rev_i, ww, j];
            }

            // Sum up the inner vector
            // remember that W[1] will ALWAYS = 0
            inner_vec = inner_vec + W[si] * mx;
          }
          // Note that since there is transfer - the calculation of R
          // comes BEFORE the transfer occurrs, because we want the
          // REGION-SPECIFIC R
          // although maybe this isn't such a bad idea
          M[t, ww, j] = R[ww, j] * inner_vec;
        }

        // And now you need 1 more step to include the influence of transfer
        // between states
        real MM[J];
        for(j in 1:J) {
          MM[j] = 0;
          for(jp in 1:J) {
            MM[j] += P[jp, j] * M[t, ww, jp];
            // P[jp, j] does NOT sum to 1 !
            // if P = [ 0.8 , 0.2
            //          0.4,  0.6 ]
            // represting that 1 -> 2 transfer is 0.2
            // and 2 -> 1 transfer is 0.4
            // then MM[1] = 0.8 * M[1] + 0.4 * M[2]
            // so the row is where the value starts
            // and the column is where it ends up
          }
        }
        for(j in 1:J) {
          M[t, ww, j] = MM[j];
        }
      }

      // ------ CARRY FORWARDS -------------
      // carry your current guess forwards in time
      // this is poor man's recursion
      // you have to do this because otherwise when you do
      // M[rev_i, ww] in the lambda loop above, the next ww
      // won't have any values in it
      if(startN < N) {
        for(j in 1:J) {
          for(wremain in min(ww + 1, max_ww):max_ww) {
            M[startN , wremain, j] = M[startN, ww, j];
          }
        }
      } //
   } // ww
}


model {

  // ------ SIGMA, BETA, and LogR and Guess ------
  // priors and sample
  // WINDOW SPECIFIC LogR
  logInitCases ~ normal(0, 1); // starting cases in each region
  xsigma ~ inv_gamma(2, 1);  // this gets the variance across the region
  for(j in 1:J) {
    xbeta[1:max_ww, j] ~ normal(0, 1);      // this gets the value for each window
    for(ww in 1:max_ww) {
      logR[ww, j] ~ normal(xbeta[ww, j], xsigma[j]);
    }
  }

  // ------ TARGET ------
  for(j in 1:J) {
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
          target += poisson_lpmf(Y[n, j] | M[n, ww, j]);
        }

    }
  }

}


