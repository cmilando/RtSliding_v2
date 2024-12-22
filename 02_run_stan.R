#' ============================================================================
#' PHASE 2
#' ============================================================================

# **** CHANGE THIS TO SEE THE IMPACT ***** #
tau        = 1
max_ww     = maxt - tau
M = matrix(data = NA, nrow = maxt, ncol = max_ww)

# --- set up sliding window ---
# STARTN = 2 so you avoid any weird beginning stuff
# You can only do this because you know what init_cases is, this is just
# used for
sliding_windows = get_SW(maxt, tau)
dim(sliding_windows)
head(sliding_windows)
tail(sliding_windows)

#' #' ============================================================================
#' ## OK NOW< RUN IN 1 D
#' stan_data <- list(
#'   N = maxt,                  # number of days
#'   tau = tau,                 #
#'   max_ww = max_ww,           #
#'   SW = sliding_windows,      #
#'   Y = NOBS,                  # cases
#'   S = length(sip),           # serial interval length
#'   W = sip                   # serial interval vector
#' )
#'
#' initf1 <- function() {
#'   list(logR = rep(0,times = max_ww))
#' }
#'
#' # if this fails, its mostly in the initialization it seems
#'
#' m_hier <- rstan::stan(file = 'sliding_1d.stan',
#'                       data = stan_data,
#'                       iter = 2000,
#'                       init = initf1,
#'                       cores = 1,
#'                       chains = 1)
#'
#' #' ============================================================================
#' stan_data_init <- list(
#'   N = maxt,                  # number of days
#'   tau = tau,                 #
#'   max_ww = max_ww,           #
#'   SW = sliding_windows,      #
#'   Y = NOBS,                  # cases
#'   S = length(sip),           # serial interval length
#'   W = sip,                   # serial interval vector
#'   init_cases = init_cases
#' )
#'
#' # if this fails, its mostly in the initialization it seems
#'
#' m_hier_init <- rstan::stan(file = 'sliding_1d_w_init.stan',
#'                       data = stan_data_init,
#'                       iter = 2000,
#'                       init = initf1,
#'                       cores = 1,
#'                       chains = 1)
#'
#' #' ============================================================================
#' stan_data_gamma <- list(
#'   N = maxt,                  # number of days
#'   tau = tau,                 #
#'   max_ww = max_ww,           #
#'   SW = sliding_windows,      #
#'   Y = NOBS,                  # cases
#'   S = length(sip),           # serial interval length
#'   W = sip                   # serial interval vector
#' )
#'
#' # if this fails, its mostly in the initialization it seems
#'
#' m_hier_gamma <- rstan::stan(file = 'sliding_1d_gamma.stan',
#'                            data = stan_data_gamma,
#'                            iter = 2000,
#'                            init = initf1,
#'                            cores = 1,
#'                            chains = 1)

#' ============================================================================
stan_data_onY <- list(
  N = maxt,                  # number of days
  tau = tau,                 #
  max_ww = max_ww,           #
  SW = sliding_windows,      #
  Y = NOBS,                  # cases
  S = length(sip),           # serial interval length
  W = sip                   # serial interval vector
)

# if this fails, its mostly in the initialization it seems

m_hier_onY <- rstan::stan(file = 'sliding_1d.stan',
                            data = stan_data_onY,
                            iter = 2000,
                            init = initf1,
                            cores = 1,
                            chains = 1)

