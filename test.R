# --- Initialize ---
set.seed(123)
xsigma = 0.1
tau = 7
maxt = 80
sip = c(0.1, 0.2, 0.3, 0.2, 0.1, 0.05, 0.05)
S = length(sip)
max_ww = maxt - tau
M = matrix(data = NA, nrow = maxt, ncol = max_ww)
MTRUE = vector("numeric", maxt)
init_cases = 100

#' ============================================================================
#' PHASE 1
#' ============================================================================


getSPLINE <- function(peak_x, peak_val,
                    max_x, startR = 1, endR = 0.2,
                  returnLOG = T) {

  require(splines)
  x = 1:max_x
  p1 <- seq(from = startR, to = peak_val, length.out = peak_x)
  p2 <- seq(from = peak_val, to = endR, length.out = max_x - peak_x + 1)
  R <- c(p1, p2[2:length(p2)])

  b2 <- lm(R ~ bs(x, knots = peak_x, degree = 3, intercept = F))
  b2.pred <- predict(b2)

  if(returnLOG) {
    return(log(b2.pred))
  } else {
    return(b2.pred)
  }

}

R_TRUE_matrix = getSPLINE(20, 1.5, maxt - 1, returnLOG = F)
plot(R_TRUE_matrix)

# -- INITIALIZE
M[1,] = init_cases
MTRUE[1] = init_cases
NOBS <- vector("numeric", maxt)

# --- GET TRUE M VALUES ----
get_Mt <- function(Rt, init_cases, w) {

  # Rt = R_TRUE_matrix
  # init_cases = init_cases
  # w = sip

  Mout <- rep(NA, length(Rt) + 1)
  Mout[1] <- init_cases

  for(t in 2:maxt) {

    inner_vec = 0
    sip_i = min(length(w), t - 1)
    for(si in 1:sip_i) {
      #print(si)
      inner_vec = inner_vec + w[si] * Mout[t - si]
    }
    rt_offset = t - 1
    Mout[t] = Rt[rt_offset] * inner_vec

  }
  return(Mout)
}

MTRUE <- get_Mt(R_TRUE_matrix, init_cases, sip)
#MTRUE

plot(MTRUE, type = 'l')
NOBS <- sapply(MTRUE, function(x) rpois(1, x))
# OBS <- sapply(MTRUE, function(x) x)
points(NOBS)

# --- Reverse out the R(t) ---
get_Rt <- function(m, w) {

  R <- rep(NA, length(m) - 1)

  for(t in 2:length(m)) {
    tau_end = min(length(w), t - 1)
    c_mat <- w[1:tau_end] %*% m[t - 1:tau_end]
    rt_offset = t - 1
    R[rt_offset] <- m[t] / c_mat
  }
  R
  return(R)
}

RT_calc = get_Rt(m = NOBS, w = sip)
plot(R_TRUE_matrix, type = 'l')
points(RT_calc, col = 'red')

#' ============================================================================
#' PHASE 2
#' ============================================================================

# --- set up sliding window ---
sliding_windows <- matrix(data = NA, nrow = max_ww, ncol = tau + 1)
window_i = 1
startN = 1
endN = tau + 1
sliding_windows[window_i, ] = startN:endN

for(window_i in 2:max_ww) {
  startN = startN + 1
  endN = endN + 1
  sliding_windows[window_i, ] = startN:endN
}

sliding_windows

for(i in 1:nrow(sliding_windows)) {
  for(j in 1:ncol(sliding_windows)) {
    sliding_windows[i,j] <- as.integer(sliding_windows[i, j])
  }
}
dim(sliding_windows)
head(sliding_windows)

## OK NOW< RUN IN 1 D
stan_data <- list(
  N = maxt,                 # number of days
  tau = tau,
  t1 = tau + 1,
  max_ww = max_ww,
  SW = sliding_windows,
  Y = NOBS,       # cases
  S = length(sip),              # serial interval length
  W = sip                      # serial interval vector
)

library(rstan)

initf1 <- function() {
  # M[N, max_ww]
  list(logR = rep(0,times = max_ww))
}

m_hier <- rstan::stan(file = 'sliding_1d.stan',
                      data = stan_data,
                      iter = 2000,
                      init = initf1,
                      cores = 1,
                      chains = 1)

#' ============================================================================
#' PHASE 3
#' ============================================================================
out <- rstan::extract(m_hier)

head(out$M[1, , ])
dim(out$M[1, , ])

# -- OBSERVATIONS
# ok for each value of M, its a function of those specific winwos
Mlist = vector('list', maxt)
dim(out$M)
for(n in 1:maxt) {
  start_W = max(1,n - tau)
  end_W = min(n, max_ww)
  Mlist[[n]] <- out$M[,n, start_W:end_W]
}

Mlist_med <- sapply(Mlist, quantile, probs = 0.5)
Mlist_lb <- sapply(Mlist, quantile, probs = 0.025)
Mlist_ub <- sapply(Mlist, quantile, probs = 0.975)

#OFFSET ALREADY ACCOUNTED FOR WELL DONE, mn above
plot(Mlist_med, type = 'l')
lines(Mlist_lb, type = 'l')
lines(Mlist_ub, type = 'l')
points(NOBS, col = 'red')

## ***********
## FIXED TO HERE 12.21.2024 at 2:50am you monster go to sleep
## THERE IS JUST MORE OFFSET BY 1 STUFF TO FIX
## ***********

# -- RT
# ok for each value of M, its a function of those specific winwos
# REMEBER THAT THESE ARE FOR TAU
dim(out$R)
Rlist_med <-apply(out$R, 2, quantile, probs = 0.5)
Rlist_lb <- apply(out$R, 2, quantile, probs = 0.025)
Rlist_ub <- apply(out$R, 2, quantile, probs = 0.975)

# REMEMBER TO OFFSET
plot(x = (tau+1):maxt, y = Rlist_med, type = 'l')
lines(x = (tau+1):maxt, y = Rlist_lb, type = 'l')
lines(x = (tau+1):maxt, y = Rlist_ub, type = 'l')
lines(x = 2:maxt, y = RT_calc, col = 'red')

## --------------------------------------------------------------


