# --- Initialize ---

# Framework:
# -- assuming that M[0] is init_cases
# -- EpiEstim doesn't calculate windows until 2-2
# -- and windows can have the same start-end day
# -- SIP has to start with 0

set.seed(123)

xsigma     = 0.1
# **** CHANGE THIS TO SEE THE IMPACT ***** #
tau        = 1
maxt       = 80
# *** THIS HAS TO START WTIH 0
sip        = c(0, 0.1, 0.2, 0.3, 0.2, 0.1, 0.05, 0.05)
S          = length(sip)
max_ww     = maxt - tau + 1
init_cases = 100

M = matrix(data = NA, nrow = maxt, ncol = max_ww)

MTRUE = vector("numeric", maxt)

#' ============================================================================
#' PHASE 1
#' ============================================================================


getSPLINE <- function(peak_x, peak_val,
                    max_x, startR = 1, endR = 0.4,
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

R_TRUE_matrix = getSPLINE(20, 1.5, maxt, returnLOG = F)
plot(R_TRUE_matrix)

# --- GET TRUE M VALUES ----
get_Mt <- function(Rt, init_cases, w) {

  Mout <- rep(0, length(Rt))

  for(t in 1:maxt) {

    inner_vec = 0

    S_loop_max = min(length(w), t + 1)
    forward_vec = 1:S_loop_max
    rev_vec = t + 1 - forward_vec

    for(si in 1:S_loop_max) {
      rev_i = rev_vec[si]

      if(rev_i < 0) {
        mx = 0
      } else if(rev_i == 0) {
        mx = init_cases
      } else {
        mx = Mout[rev_i]
      }

      inner_vec = inner_vec + w[si] * mx
    }
    Mout[t] = Rt[t] * inner_vec

  }
  Mout
  return(Mout)
}

MTRUE <- get_Mt(R_TRUE_matrix, init_cases, sip)
MTRUE
stopifnot(!any(is.na(MTRUE)))
stopifnot(all((MTRUE != 0)))

plot(MTRUE, type = 'l')
NOBS <- sapply(MTRUE, function(x) rpois(1, x))
points(NOBS)
stopifnot(all((NOBS != 0)))

# --- Reverse out the R(t) ---
get_Rt <- function(m, init_cases, w) {

  Rt <- rep(0, length(m))

  for(t in 1:maxt) {

    inner_vec = 0

    S_loop_max = min(length(w), t + 1)
    forward_vec = 1:S_loop_max
    rev_vec = t + 1 - forward_vec

    for(si in 1:S_loop_max) {
      rev_i = rev_vec[si]

      if(rev_i < 0) {
        mx = 0
      } else if(rev_i == 0) {
        mx = init_cases
      } else {
        mx = m[rev_i]
      }

      inner_vec = inner_vec + w[si] * mx
    }
    Rt[t]  = m[t] / inner_vec

  }
  Mout
  return(Rt)
}

RT_calc = get_Rt(m = NOBS, init_cases = init_cases, w = sip)
length(RT_calc)
length(R_TRUE_matrix)
plot(RT_calc, col = 'red', type = 'l')
lines(R_TRUE_matrix, type = 'l')

#' ============================================================================
#' PHASE 2
#' ============================================================================

# --- set up sliding window ---
# STARTN = 2 so you avoid any weird beginning stuff
sliding_windows <- matrix(data = NA, nrow = max_ww, ncol = tau)
window_i = 1
startN = 1
endN = startN + (tau-1)
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
tail(sliding_windows)

## OK NOW< RUN IN 1 D
stan_data <- list(
  N = maxt,                  # number of days
  tau = tau,                 #
  max_ww = max_ww,           #
  SW = sliding_windows,      #
  Y = NOBS,                  # cases
  S = length(sip),           # serial interval length
  W = sip,                   # serial interval vector
  init_cases = init_cases
)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

initf1 <- function() {
  # M[N, max_ww]
  list(logR = rep(0,times = max_ww))
}

# if this fails, its mostly in the initialization it seems

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
tail(out$M[1, , ])
dim(out$M[1, , ])

# -- OBSERVATIONS
# ok for each value of M, its a function of those specific winwos
Mlist = vector('list', maxt)
dim(out$M)
for(n in 1:maxt) {
  start_W = max(1,n - tau+1)
  end_W = min(n, max_ww)
  #print(paste(n, start_W, end_W))
  Mlist[[n]] <- out$M[,n, unique(start_W:end_W)]
}

Mlist[[1]]

Mlist_med <- sapply(Mlist, quantile, probs = 0.5)
Mlist_lb <- sapply(Mlist, quantile, probs = 0.025)
Mlist_ub <- sapply(Mlist, quantile, probs = 0.975)

#OFFSET ALREADY ACCOUNTED FOR WELL DONE, mn above
plot(Mlist_med, type = 'l')
lines(Mlist_lb, type = 'l')
lines(Mlist_ub, type = 'l')
points(NOBS, col = 'red')

# -- RT
# ok for each value of M, its a function of those specific winwos
# REMEBER THAT THESE ARE FOR TAU
dim(out$R)
#out$R <- out$R[, 2:(dim(out$R)[2])]
Rlist_med <-apply(out$R, 2, quantile, probs = 0.5)
Rlist_lb <- apply(out$R, 2, quantile, probs = 0.025)
Rlist_ub <- apply(out$R, 2, quantile, probs = 0.975)

# REMEMBER TO OFFSET
plot(x = (tau):maxt, y = Rlist_med, type = 'l')
lines(x = (tau):maxt, y = Rlist_lb, type = 'l')
lines(x = (tau):maxt, y = Rlist_ub, type = 'l')
lines(x = 1:maxt, y = RT_calc, col = 'red')
RT_calc

## CANT GET TO LOW CASES

#' ============================================================================
#' Plots
#' ============================================================================

# Adjust graphics parameters for side-by-side plots
par(mfrow = c(2, 1))  # Two plots in one row

# Plot 1: M(t)
plot(Mlist_med, type = 'l', col = 'blue',
     ylab = 'M(t) Values', xlab = NA, lwd = 1.5,
     main = 'M(t) with CI')
lines(Mlist_lb, type = 'l', col = 'blue', lty = 2)
lines(Mlist_ub, type = 'l', col = 'blue', lty = 2)
points(NOBS, col = 'black')

# Plot 2: R(t)
plot(x = tau:maxt, y = Rlist_med, type = 'l',
     col = 'green', ylab = 'R(t) Values', xlab = 'Time',
     lwd = 2, main = 'R(t) with CI')
# lines(x = (tau):maxt, y = Rlist_lb, type = 'l', col = 'green', lty = 2)
# lines(x = (tau):maxt, y = Rlist_ub, type = 'l', col = 'green', lty = 2)
lines(x = 1:maxt, y = RT_calc, col = 'red', lwd = 1)


# Plot 2b: compared with epiEstim
# ok so you need to include init_cases
incid = c(NOBS)
this_si = sip

t_start <- seq(2, length(incid)-(tau - 1))
t_start
t_end <- t_start + (tau - 1)
t_end
library(EpiEstim)
res <- estimate_R(incid = incid,
                  method = "non_parametric_si",
                  config = make_config(list(
                    si_distr = this_si,
                    t_start = t_start,
                    t_end = t_end)),
                  backimputation_window = 10)
#plot(res)

epiE_x = res$R$t_end
epiE_med = res$R$`Median(R)`
epiE_lb = res$R$`Quantile.0.025(R)`
epiE_ub = res$R$`Quantile.0.975(R)`

lines(x = epiE_x, y = epiE_med, type = 'l',
     col = 'purple')

# lines(x = epiE_x, y = epiE_lb, type = 'l',
#       col = 'purple', lty = 2)
#
# lines(x = epiE_x, y = epiE_ub, type = 'l',
#       col = 'purple', lty = 2)

abline(v = S)

# Reset graphical parameters to default
par(mfrow = c(1, 1))


