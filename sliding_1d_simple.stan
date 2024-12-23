data {
  int<lower=1> N_windows;       // the number of windows.
  int<lower=1> N_obs;           // the number of observations
  int<lower=0> Y[N_obs];        // observed cases in windows.
  int<lower=1> S;               // length of serial interval
  vector[S] W;                  // serial interval
  int<lower=0> SW[N_windows,2]; // sliding window matrix
  matrix[N_obs,N_obs] A;        // sliding window value
  int<lower=1> tau;
}

parameters {
  // time and region specific R values, in log space
  real logR[N_windows  - 1];
  real xbeta[N_windows - 1];
  real<lower=0> xsigma;
}

transformed parameters {

  // expected value of cases in each WINDOW
  vector<lower=0>[N_windows] M_window = rep_vector(0.001, N_windows);

  // the back half of M_window for matrix math step at the end
  vector<lower=0>[N_obs] M_for_matrix   = rep_vector(0.001, N_obs);

  // expected value of cases in each TIMESTEP
  vector<lower=0>[N_obs] M = rep_vector(0.01, N_obs);

  // no matter what N is, R will always be 1 less
  // than the number of Y
  // ALWAYS
  vector<lower=0>[N_windows - 1] R = rep_vector(0.1, N_windows - 1);

  // ------ CALCULATE R(t) and INITIALIZE M[1] -------------
  // get R in exp() space for each window
  for(w in 1:(N_windows - 1)) {
    R[w] = exp(logR[w]);
  }

  // Initialize with Y_w !! Not Y
  M_window[1] = 0;
  for(n in SW[1,1]:SW[1,2]) {
    M_window[1] += Y[n];
  }

  //print("-----");
  //print("R b4: ", R);
  //print("M_window b4: ", M_window);

  // ------ CALCULATE M(t) -------------
  /// NOW FOR EACH BLOCK
  for(t in 2:N_windows) {

    ///
    real inner_vec = 0;
    int  S_loop_max = min(S, t + 1);
    int  forward_vec[S_loop_max];
    int  rev_vec[S_loop_max];

    ///
    for(si in 1:S_loop_max) {
      forward_vec[si] = si;
      rev_vec[si] = t + 1 - si;
    }

    ///
    real mx = 0;

    //print("****");
    //print("t: ", t);
    //print("forward_vec: ", forward_vec);
    //print("rev_vec: ", rev_vec);

    ///
    for(si in 1:S_loop_max) {
      //print(">> si = ", si);
      int rev_i = rev_vec[si];

      if(rev_i < 1) {
        mx = 0;
      } else {
        mx = M_window[rev_i];
      }

      //print("W[si]: ", W[si]);
      //print("mx: ", mx);
      //print("W[si] * mx: ", W[si] * mx);

      inner_vec = inner_vec + W[si] * mx;
      //print("inner_vec: ", inner_vec);
    }

    ///
    int wi = t - 1;
    M_window[t] = R[wi] * inner_vec;
    //print("M_window: ", M_window);

  }
  //print("@@");
  //print("M_window after: ", M_window);

  // ------ REVERSE OUT M(t) -------------

  // FIRST, YOU NEED TO ADD UP ANY CASES YOU DONT HAVE ALREADY
  // TO THE FRONT OF M_WINDOW
  if(tau > 1) {
    for(n in 1:(tau-1)) {
      M_for_matrix[n] = M[n]; // This needs to be expectation or it solves
    }
    for(n in tau:N_obs) {
      M_for_matrix[n] = M_window[n - tau + 1];
    }
  } else {
    for(n in 1:N_obs) {
      M_for_matrix[n] = M_window[n];
    }
  }

  //print("M_for_matrix ", M_for_matrix);

  // Solve the linear system A[N x N] * M_matrix[N x 1] = M[N x 1]
  M = mdivide_left(A, M_for_matrix);

  //print("M ", M);


}


// BELOW HERE IS GOOD
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
  // This is what drives the solution
  // and fixes the beginning part ... lolol
  for(n in 1:N_obs) {
    target += poisson_lpmf(Y[n] | M[n]);
  }

}

