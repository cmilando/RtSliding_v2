#' ============================================================================
#' PHASE 2
#' ============================================================================

# **** CHANGE THIS TO SEE THE IMPACT ***** #
tau        = 1
max_ww     = maxt - tau
M = matrix(data = NA, nrow = maxt, ncol = max_ww)

# --- set up sliding window ---
# STARTN = 2 so you avoid any weird beginning stuff
# You can only do this because you know what init_cases is, this is just
# used for
sliding_windows = get_SW(maxt, tau)
dim(sliding_windows)
head(sliding_windows)
tail(sliding_windows)

# -- Also get reveresed sliding windows --
SWT = get_SWT(sliding_windows, tau, maxt)
dim(SWT)
head(SWT)

#' ============================================================================
stan_data_onY <- list(
  N = maxt,                  # number of days
  tau = tau,                 #
  max_ww = max_ww,           #
  SW = sliding_windows,      #
  SWT = SWT,
  Y = NOBS,                  # cases
  S = length(sip),           # serial interval length
  W = sip                   # serial interval vector
)

# if this fails, its mostly in the initialization it seems
m_hier_onY <- rstan::stan(file = 'sliding_1d_simple.stan',
                            data = stan_data_onY,
                            iter = 2000,
                            cores = 1,
                            chains = 1)

