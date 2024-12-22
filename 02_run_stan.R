#' ============================================================================
#' PHASE 2
#' ============================================================================

# **** CHANGE THIS TO SEE THE IMPACT ***** #
tau        = 2
max_ww     = maxt - tau + 1

# --- set up sliding window ---
# STARTN = 2 so you avoid any weird beginning stuff
# You can only do this because you know what init_cases is, this is just
# used for
sliding_windows = get_SW(maxt, tau)

dim(sliding_windows)
head(sliding_windows)
tail(sliding_windows)


Y_w = vector("numeric", length = max_ww)

for(w in 1:max_ww) {
  startN = sliding_windows[w, 1]
  endN = sliding_windows[w, 2]
  Y_w[w] = sum(NOBS[startN:endN])
}
head(Y_w)
tail(Y_w)
tail(NOBS)

# Can you recover Y from Y_w and max_ww
# yes because you know what M[1] is
Y_recovered = vector("numeric", length = length(NOBS))

A <- matrix(0, nrow = max_ww, ncol = maxt)

# Fill the design matrix
for (w in 1:max_ww) {
  startN <- sliding_windows[w, 1]
  endN <- sliding_windows[w, 2]
  A[w, startN:endN] <- 1
}

# ok now add your knowns
A0 = matrix(0, nrow = 1, ncol = maxt)
A0[1, 1] = 1
Y0 = NOBS[1]

Afull = rbind(A0, A)
Yfull = c(Y0, Y_w)

NOBS_recovered <- solve(Afull, Yfull)

NOBS_recovered

#' ============================================================================
stan_data <- list(
  N_windows = N_w,        # number of windows
  N_obs = maxt,           # number of observations
  Y = Y_w,                # cases
  S = length(sip),        # serial interval length
  W = sip                 # serial interval vector,
)

# if this fails, its mostly in the initialization it seems
m_hier_onY <- rstan::stan(file = 'sliding_1d_simple.stan',
                            data = stan_data,
                            iter = 2000,
                            cores = 1,
                            chains = 1)


## wow that seemed to work
## so now, just move the summing inside STAN so you can get an M prediction
## and you're good
