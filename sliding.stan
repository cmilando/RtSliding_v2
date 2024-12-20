data {
  int<lower=1> N;          // the number of observations.
  int<lower=1> SWcol;      // sw_col -- aka sliding window size
  int<lower=1> SWrow;      // sw_row
  int SW[SWrow,SWcol];     // sliding window matrix
  int<lower=1> J;          // the number of regions
  int<lower=0> Y[N,J];     // observed cases.
  matrix[N*J,J] P;         // transfer matrix, changes with time
  int<lower=1> S;          // length of serial interval
  vector[S] W;             // serial interval
  vector[J] init_cases;    // initial cases
}

parameters {
  vector<lower=0.00001>[J] xsigma;  // region-specific st-dev
  matrix[SWrow,J] xbeta;            // time-region specific beta, BY WINDOW
  matrix[N,J] logR;                // time-region specific R values, in Log space, BY DAY
}

transformed parameters {
  matrix[N+1,J] M;           // expected value of cases
  matrix[N,J] R;           // time and region specific R values, expressed normally
  matrix[J, J] RR;         // Diag R matrix

  // ------ CALCULATE R(t) -------------
  // get R in exp() space
  for(j in 1:J) {
    for(n in 1:N) {
      R[n, j] = exp(logR[n, j]);
    }
  }

  // ------ EQUATION (11b) -------------
  // initialize M at t = 0
  for(j in 1:J) {
    M[1, j] = init_cases[j];
  }

  // calculate m[t,j] for t >= 1
  for(n in 1:N) {

    // M needs its own iterator
    int mn = n + 1;

    // calculate the inner part first
    int tau_end = min(S, n);

    // Create diagonal matrix RR
    RR = diag_matrix(to_vector(R[n, ]));

    // Create MM based on t value
    // MM is m(t - 1), ..., m(1)
    // where rows are regions
    // and columns are time
    matrix[J, tau_end] MM;
    for(tt in 1:tau_end) {
      MM[, tt] = to_vector(M[mn - tt, ]);
    }

    // Create WW vector
    matrix[tau_end, 1] WW;
    WW[1:tau_end, 1] = to_vector(W[1:tau_end]);

    // Calculate result
    int start_P = J*(n - 1) + 1;
    int end_P   = J*(n - 1) + J;
    M[mn, ] = to_row_vector(P[start_P:end_P, ]' * RR * MM * WW);

  }

}

model {


  // ------ EQUATION (11a) -------
  // MODIFIED BY CWM TO APPLY SLIDING WINDOW TO R(t) INSTEAD OF M
  // priors and sample
  xsigma ~ inv_gamma(2, 1);

  for(j in 1:J) {

      // weak prior on beta
      xbeta[, j] ~ normal(0, 1);

      // sample logR
      for(ww in 1:SWrow) {
        // beginning of the sliding window
        int start_W = SW[ww, 1];
        // end of the sliding window
        int end_W   = SW[ww, SWcol];
        // each block of R(t)
        // so each block of R(t) belongs to a different xbeta
        // each xbeta is for each window
        logR[start_W:end_W, j] ~ normal(xbeta[ww, j], xsigma[j]);
      }

  }

  // ------ EQUATION (11c) -------
  for(j in 1:J) {
      for(n in 1:N) {
        int mn = n+1;
        Y[n, j] ~ poisson(M[mn, j]);
      }
  }

}

