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
maxt       = 30
init_cases = 100
sip        = si(4, 4, 2)
S          = length(sip)
MTRUE       = vector("numeric", maxt)

##
stopifnot(sip[1] == 0)

##
R_TRUE_matrix = getSPLINE(ceiling(maxt/4), 1.5, maxt)

##
MTRUE <- get_Mt(R_TRUE_matrix, init_cases, sip)

NOBS <- sapply(MTRUE, function(x) rpois(1, x))
stopifnot(all((NOBS != 0)))

RT_calc = get_Rt(m = NOBS, init_cases = init_cases, w = sip)
stopifnot(length(RT_calc) == length(R_TRUE_matrix))

## add the back-calculated cases
# back_imputed <- backimpute(NOBS, 4)
# NOBS <- c(back_imputed, NOBS)
# maxt <- length(NOBS)
