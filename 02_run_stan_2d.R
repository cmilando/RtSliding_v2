#' ============================================================================
#' PHASE 2
#' ============================================================================

source("01_simulate_data_2d.R")

# **** CHANGE THIS TO SEE THE IMPACT ***** #
tau        = 5
max_ww     = maxt - tau

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

## ----------
J = 2

NOBS_mat <- cbind(x1$NOBS, x2$NOBS)

head(NOBS_mat)

P_base <- matrix(c(c(1.0, 0.0), c(0.0, 1.0)), nrow = 2, byrow = T)

p_mod <- expand.grid(j1 = seq(0, 1, by = 0.1),
            j2 = seq(0, 1, by = 0.1))

xpi = 2

P <- P_base

P[1, 1] <- P_base[1, 1] - p_mod[xpi, 1]
P[1, 2] <- P_base[1, 2] + p_mod[xpi, 1]
P[2, 1] <- P_base[2, 1] - p_mod[xpi, 2]
P[2, 2] <- P_base[2, 2] + p_mod[xpi, 2]

P <- matrix(c(c(0.8, 0.2), c(0.4, 0.6)), nrow = 2, byrow = T)
P <- matrix(c(c(0.9, 0.1), c(0.1, 0.9)), nrow = 2, byrow = T)
P <- matrix(c(c(0.7, 0.3), c(0.1, 0.9)), nrow = 2, byrow = T)
P <- matrix(c(c(0.7, 0.3), c(0.3, 0.7)), nrow = 2, byrow = T)

matrix_to_mac_filename(P)

stopifnot(all(rowSums(P) == 1))

# --- OK NOW< RUN IN 1 D ---
stan_data <- list(
  N = maxt,                  # number of days
  J = J,                     # number of zones
  P = P,                     # transfer matrix
  tau = tau,                 # sliding window size
  max_ww = max_ww,           # max number of windows
  SW = sliding_windows,      #
  Y = NOBS_mat,              # cases
  S = length(sip),           # serial interval length
  W = sip                   # serial interval vector
)

## // if tau = 1 and GuessM = 0, EpiEstim and our estimate line up

# if this fails, its mostly in the initialization it seems
initf1 <- function() {
  list(logR = matrix(0, nrow = max_ww, ncol = J))
}

# Run STAN
m_hier <- rstan::stan(file = 'sliding_2d.stan',
                      data = stan_data,
                      iter = 3000,
                      init = initf1,
                      cores = 1,
                      chains = 1)

# interestingly, just changing to matrices made it slower
# but results seem the same

source("03_extract_output_2d.R")

