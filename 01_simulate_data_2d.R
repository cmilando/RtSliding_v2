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

sip        = si(4, 4, 2)
S          = length(sip)
maxt       = 30
stopifnot(sip[1] == 0)

# --- Initialize ---
init_cases = 100
R_TRUE_matrix = getSPLINE(ceiling(maxt/6), 1.5, maxt)
MTRUE <- get_Mt(R_TRUE_matrix, init_cases, sip)
NOBS <- sapply(MTRUE, function(x) rpois(1, x))
stopifnot(all((NOBS != 0)))
RT_calc = get_Rt(m = NOBS, init_cases = init_cases, w = sip)
stopifnot(length(RT_calc) == length(R_TRUE_matrix))

x1 <- data.frame(R_TRUE_matrix, MTRUE, NOBS, RT_calc)
x1$J = 1

# --- Initialize ---
init_cases = 150
R_TRUE_matrix = getSPLINE(ceiling(3*maxt/4), 1.2, maxt)
MTRUE <- get_Mt(R_TRUE_matrix, init_cases, sip)
NOBS <- sapply(MTRUE, function(x) rpois(1, x))
stopifnot(all((NOBS != 0)))
RT_calc = get_Rt(m = NOBS, init_cases = init_cases, w = sip)
stopifnot(length(RT_calc) == length(R_TRUE_matrix))

x2 <- data.frame(R_TRUE_matrix, MTRUE, NOBS, RT_calc)
x2$J = 2

