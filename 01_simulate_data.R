#' ============================================================================
#' Framework:
#' -- assuming that M[0] is init_cases
#' -- EpiEstim doesn't calculate windows until 2-2
#' -- and windows can have the same start-end day
#' -- SIP has to start with 0
#' -- R(t) represents how cases go from m(t -1) to m(t)
#' ============================================================================

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("00_functions.R")
set.seed(123)

#' ============================================================================
#' PHASE 1
#' ============================================================================

# --- Initialize ---
maxt       = 40
init_cases = 100
sip        = c(0, 0.1, 0.2, 0.3, 0.2, 0.1, 0.05, 0.05)
S          = length(sip)
MTRUE       = vector("numeric", maxt)

##
stopifnot(sip[1] == 0)

## set graph params
par(mfrow = c(1, 2))  # Two plots in one row

##
R_TRUE_matrix = getSPLINE(ceiling(maxt/4), 1.5, maxt)

##
MTRUE <- get_Mt(R_TRUE_matrix, init_cases, sip)
NOBS <- sapply(MTRUE, function(x) rpois(1, x))
stopifnot(all((NOBS != 0)))

plot(MTRUE, type = 'l', ylab = 'Cases', xlab = 'Day',
     col = 'red', lwd = 1.5,
     ylim = range(MTRUE, NOBS),
     main = paste0('Cases (N=', length(R_TRUE_matrix), ")"))
points(NOBS, cex = 0.8)

legend("topright", legend = c("Function", "Observed"),
       col = c("red", "black"),
       lty = 1, lwd = 2, cex = 0.8)

##
RT_calc = get_Rt(m = NOBS, init_cases = init_cases, w = sip)
stopifnot(length(RT_calc) == length(R_TRUE_matrix))
plot(RT_calc, col = 'black', type = 'l', ylab = 'R(t)', xlab = 'Day',
     ylim = range(RT_calc, R_TRUE_matrix),
     main = paste0('Rt (N=', length(R_TRUE_matrix), ")"), lty = 3)
lines(R_TRUE_matrix, type = 'l', col = 'red')
abline(h = 1, lty = 2)
legend("topright", legend = c("Function", "Observed"),
       col = c("red", "black"),
       lty = 1, lwd = 2, cex = 0.8)

# Reset graphical parameters to default
par(mfrow = c(1, 1))



