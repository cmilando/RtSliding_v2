#' ============================================================================
#' PHASE 2
#' ============================================================================

# **** CHANGE THIS TO SEE THE IMPACT ***** #
tau        = 4
max_ww     = maxt - tau + 1

# --- set up sliding window ---
# STARTN = 2 so you avoid any weird beginning stuff
# You can only do this because you know what init_cases is, this is just
# used for
sliding_windows = get_SW(maxt, tau)
sliding_windows

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
#
# # Can you recover Y from Y_w and max_ww
# # yes because you know what M[1] is
# Y_recovered = vector("numeric", length = length(NOBS))

A <- matrix(0, nrow = max_ww, ncol = maxt)

# Fill the design matrix
for (w in 1:max_ww) {
  startN <- sliding_windows[w, 1]
  endN <- sliding_windows[w, 2]
  A[w, startN:endN] <- 1
}

# ok now add your knowns
if(tau > 1){
  A0 = matrix(0, nrow = (tau - 1), ncol = maxt)
  A = rbind(A0, A)
  # then just fill in as diagonal matrix
  for(tt in 1:(tau-1)) {
    A[tt, tt] = 1
  }
}
dim(A)
stopifnot(ncol(A) == nrow(A))
A
A

if(tau > 1) {
  Yfull <- Y_w
  for(tt in 1:(tau-1)) {
    Yfull = c(NOBS[tt], Yfull)
  }
} else {
  Yfull <- NOBS
}
NOBS
Yfull
#
NOBS_recovered <- solve(A, Yfull)
#
NOBS_recovered
NOBS

#solve(A, c(8,30,8.15485,18.5264))

#' ============================================================================
stan_data <- list(
  N_windows = max_ww,        # number of windows
  N_obs = maxt,           # number of observations
  Y = NOBS,                # cases
  S = length(sip),        # serial interval length
  W = sip,                # serial interval vector,
  SW = sliding_windows,
  A = A,
  tau = tau
)


# if this fails, its mostly in the initialization it seems
# just has to be high enough to not fail
initf1 = function() {
   list(logR = rep(4, times = (max_ww-1)))
}


m_hier_onY <- rstan::stan(file = 'sliding_1d_simple.stan',
                            data = stan_data,
                            iter = 2000,
                            init = initf1,
                            cores = 1,
                            chains = 1)

# m_hier_onY <- rstan::stan(file = 'sliding_1d_simple_fixedTop.stan',
#                           data = list(
#                             N = max_ww,
#                             Y = Y_w,
#                             S = length(sip),
#                             W = sip
#                           ),
#                           iter = 2000,
#                           cores = 1,
#                           chains = 1)

dim(out$M)

source("03_plot_output.R")
