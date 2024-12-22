data {
  int<lower=1> N;          // the number of observations.
  int<lower=1> tau;        // sliding window size
  int<lower=1> max_ww;     // sw_row
  int SW[max_ww,2];        // sliding window matrix
  int SWT[N,tau];        // sliding window matrix
  int<lower=0> Y[N];       // observed cases.
  int<lower=1> S;          // length of serial interval
  vector[S] W;             // serial interval
}

parameters {
  // time and region specific R values, in log space
  real logR[max_ww];
  real xbeta[max_ww];
  real<lower=0> xsigma;
  real<lower=0, upper=1> Umatrix[N,tau];
}

transformed parameters {

  // expected value of cases, **J** each gets its own window, lots of zeros
  real<lower=0> M[N] = rep_array(0.1, N);
  real<lower=0> R[max_ww] = rep_array(1, max_ww);

  // ------ CALCULATE R(t) and INITIALIZE M[1] -------------
  // get R in exp() space for each window
  for(ww in 1:max_ww) {
    R[ww] = exp(logR[ww]);
  }

  M[1] = Y[1] * 1.0;

  // ------ USE RECURSION TO CALCULATE M(t) -------------
  // YOU PROBABLY NEED TO WRITE A FUNCTION AND THEN CALL IT IN THE LOOP

  for(t in 2:N) {

    // *******************
    // maybe the first step is
    // which Rs reach back to this point
    // sort of an inversion of the SW matrix
    // count down from the top
    // you could just choose one at random to have an influence
    // and then over time it would converge?
    // Honestly, not a terrible idea
    real non_empty[tau];
    int tx = 0;

    for(tau_i in 1:tau) {
      if(SWT[t, tau_i] > 0) {
        tx = 1;
      } else {
        tx = 0;
      }
      non_empty[tau_i] = tx * Umatrix[t, tau_i];
    }

    int max_U = 0;
    real trueMax = -1;

    for(tau_i in 1:tau) {
      if(non_empty[tau_i] > trueMax) {
        trueMax = non_empty[tau_i];
        max_U = tau_i;
      }
    }

    //print("t = ", t);
    //print("max_U = ", max_U);
    //print("non empty = ", non_empty);
    //print("Umatrix = ", Umatrix[t, 1:tau]);

    int random_w = SWT[t, max_U];
    real Rstar = R[random_w];
    // *******************

    // Now, you can continue as normal
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

    // apply this randomly chosen Rstar from the list of potential Rwws
    M[t] = Rstar * inner_vec;

  }

}


// THIS IS GOOD DON'T CHANGE IT
model {

  // ------ SIGMA, BETA, and LogR ------
  // priors and sample
  xsigma ~ inv_gamma(2, 1);
  xbeta ~ normal(0, 1);

  // WINDOW SPECIFIC LogR
  for(ww in 1:max_ww) {
    logR[ww] ~ normal(xbeta[ww], xsigma);
  }

  for(tau_i in 1:tau) {
    for(n in 1:N) {
      Umatrix[n, tau_i] ~ uniform(0.0, 1.0);
    }
  }

  // SO THIS ENFORCES THAT THE OBSERVED CASES
  // ARE DRAWN FROM ALL OF THE WINDOWS THAT THIS DAY FALLS INTO
  // This is also what forces the betas of adjacent windows to be
  // related. and you have to add each target separately so its 1:1
  for (n in 1:N) {
    Y[n] ~ poisson(M[n]);
  }

}


