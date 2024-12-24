#' ============================================================================
#' PHASE 2
#' ============================================================================

source("01_simulate_data.R")

# **** CHANGE THIS TO SEE THE IMPACT ***** #
tau        = 5
max_ww     = maxt - tau
M = matrix(data = NA, nrow = maxt, ncol = max_ww)

# --- set up sliding window ---
# STARTN = 2 so you avoid any weird beginning stuff
# You can only do this because you know what init_cases is, this is just
# used for
sliding_windows = get_SW(maxt, tau)

# adding small windows up front
if(tau > 1) {
  second_col <- 2:(tau)

  sliding_windows_beta <- matrix(data = c(rep(2, length(second_col)),
                                          second_col),
                          nrow = length(second_col),
                          ncol = 2)

  sliding_windows <- rbind(sliding_windows_beta, sliding_windows)
  max_ww <- nrow(sliding_windows)
}

# --- OK NOW< RUN IN 1 D ---
stan_data <- list(
  N = maxt,                  # number of days
  tau = tau,                 #
  max_ww = max_ww,           #
  SW = sliding_windows,      #
  Y = NOBS,                  # cases
  S = length(sip),           # serial interval length
  W = sip                   # serial interval vector
)

## // if tau = 1 and GuessM = 0, EpiEstim and our estimate line up

initf1 <- function() {
  list(logR = rep(0,times = max_ww))
}

# if this fails, its mostly in the initialization it seems

m_hier <- rstan::stan(file = 'sliding_1d.stan',
                      data = stan_data,
                      iter = 2000,
                      init = initf1,
                      cores = 1,
                      chains = 1)

source("03_extract_output.R")
source("04_plot.R")
