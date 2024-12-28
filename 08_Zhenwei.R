library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("00_functions.R")
set.seed(123)

#' ============================================================================
#' Data
#' ============================================================================
NOBS_mat <- read.csv("example_data.csv")
head(NOBS_mat)

# remove the first row, as these are the initial cases
NOBS_mat <- NOBS_mat[2:nrow(NOBS_mat), ]
plot(NOBS_mat[,1], type = 'l')
lines(NOBS_mat[,2], col = 'red')
lines(NOBS_mat[,3], col = 'green')

maxt = nrow(NOBS_mat)
J = ncol(NOBS_mat)

#' ============================================================================
#' Serial Interval
#' ============================================================================
si_shape <- 2 ## shape parameter for serical interval assuming gamma distribution
si_rate <- 0.5 ## rate parameter for serical interval assuming gamma distribution
si_t_max <- 14 ## maximum number of days with non-zero probability for serial interval
w <- sapply(1:si_t_max, function(x){
  pgamma(x, si_shape, si_rate) - pgamma(x-1, si_shape, si_rate)
})

w <- c(0, w/sum(w))
sip <- w
#' ============================================================================
#' Transfer Matrix
#' ============================================================================
P <- matrix(c(0.8, 0.2, 0.1,  ## mobility matrix
              0.15, 0.6, 0.3,
              0.05, 0.2, 0.6), 3)
all(rowSums(P) == 1)

#' ============================================================================
#' Sliding Windows
#' ============================================================================
tau        = 40
max_ww     = maxt - tau

sliding_windows = get_SW(maxt, tau)
first_break = mean(sliding_windows[1, ])
second_break = mean(sliding_windows[max_ww, ])

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

#' ============================================================================
#' Run Stan
#' ============================================================================
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
                      cores = 4,
                      chains = 1)
