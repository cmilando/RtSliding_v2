#' ============================================================================
#' PHASE 2
#' ============================================================================

# **** CHANGE THIS TO SEE THE IMPACT ***** #
tau        = 1
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

N_w = max_ww

N_w

#' ============================================================================
stan_data <- list(
  N = N_w,         # number of days
  Y = Y_w,                 # cases
  S = length(sip),          # serial interval length
  W = sip                  # serial interval vector
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
