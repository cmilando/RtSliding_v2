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
first_break = mean(sliding_windows[1, ])
second_break = mean(sliding_windows[max_ww, ])

# adding small windows up front and the back
if(tau > 1) {
  second_col <- 2:(tau)
  first_col <- (maxt-tau+2):maxt
  stopifnot(length(first_col) == (length(second_col)))

  sliding_windows_top <- matrix(data = c(rep(2, length(second_col)),
                                          second_col),
                          nrow = length(second_col),
                          ncol = 2)

  sliding_windows_bottom <- matrix(data = c(first_col,
                                         rep(maxt, length(first_col))),
                                nrow = length(first_col),
                                ncol = 2)

  sliding_windows <- rbind(sliding_windows_top,
                           sliding_windows,
                           sliding_windows_bottom)

  max_ww <- nrow(sliding_windows)
}

sliding_windows

## ----------
J = 2

NOBS_mat <- cbind(x1$NOBS, x2$NOBS)

head(NOBS_mat)

P <- matrix(c(c(1.0, 0.0), c(0.0, 1.0)), nrow = 2, byrow = T)
P

stopifnot(all(rowSums(P) == 1))

## ---------
# Include a reporting delay at the tail
# problems here if this is different from tau
# or not problems but just changes where the discontinuity is
# but I believe this works as epxected
# because things are changing but you have less information
# this is probably where `nowcasting` comes into play
NDELAY = 2
stopifnot(NDELAY<tau)
report_delay <- c(0.6, 0.3)

# ok now create 100 different y tails and pass that in
set.seed(123)
NSIM <- 5
ar <- array(rep(NaN, NSIM*(NDELAY)*J), c(NSIM, NDELAY, J))
for(j in 1:J) {
  observed_Y <- NOBS_mat[(maxt - (NDELAY-1):0), j]
  expected_Y <- observed_Y / report_delay
  residual_Y = expected_Y - observed_Y
  sample_Y <- sapply(residual_Y, function(yy) rpois(NSIM, yy))
  ar[, , j] = as.integer(observed_Y + sample_Y)
}


# the solution here is to do like EpiNow2 and separate this section out

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
  W = sip,                   # serial interval vector
  NDELAY = NDELAY,
  NSIM = NSIM,
  Ymod = ar
)

## // if tau = 1 and GuessM = 0, EpiEstim and our estimate line up

# if this fails, its mostly in the initialization it seems
initf1 <- function() {
  list(logR = matrix(0, nrow = max_ww, ncol = J))
}

# Run STAN
m_hier <- rstan::stan(file = 'sliding_2d_plus.stan',
                      data = stan_data,
                      iter = 3000,
                      init = initf1,
                      cores = 1,
                      chains = 1)

# interestingly, just changing to matrices made it slower
# but results seem the same

source("03_extract_output_2d.R")


