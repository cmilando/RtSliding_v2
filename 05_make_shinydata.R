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

P_base <- matrix(c(c(1.0, 0.0), c(0.0, 1.0)), nrow = 2, byrow = T)

p_mod <- expand.grid(j1 = seq(0, 0.9, by = 0.1),
                     j2 = seq(0, 0.9, by = 0.1))

#####

for(xpi in 1:nrow(p_mod)) {
  timestamp(suffix = paste(" >", paste(p_mod[xpi,], collapse = '-')))
  P <- P_base

  P[1, 1] <- P_base[1, 1] - p_mod[xpi, 1]
  P[1, 2] <- P_base[1, 2] + p_mod[xpi, 1]

  P[2, 1] <- P_base[2, 1] + p_mod[xpi, 2]
  P[2, 2] <- P_base[2, 2] - p_mod[xpi, 2]

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

  f <- paste0("shinydata/x",matrix_to_mac_filename(P),".RData")
  out <- rstan::extract(m_hier)

  ###
  MlistJ = vector("list", J)
  SWT <- get_SWT(sliding_windows, maxt)
  for(j in 1:J) {

    Mlist = vector('list', maxt)
    dim(out$M)

    for(n in 1:maxt) {
      these_windows <- SWT[[n]]
      Mlist[[n]] <- as.vector(out$M[ , n, these_windows, j])
    }

    Mlist_x <- 1:maxt
    Mlist_med <- sapply(Mlist, quantile, probs = 0.5)
    Mlist_lb <- sapply(Mlist, quantile, probs = 0.025)
    Mlist_ub <- sapply(Mlist, quantile, probs = 0.975)

    stopifnot(length(Mlist_med) == maxt)

    MlistJ[[j]] = data.frame(
      J = j, NOBS = NOBS_mat[, j],
      Mlist_x, Mlist_med, Mlist_lb, Mlist_ub
    )
  }

  Mlist_df  <- do.call(rbind, MlistJ)
  Mlist_df$J <- factor(Mlist_df$J)

  ###
  library(EpiEstim)
  RlistJ = vector("list", J)
  for(j in 1:J) {

    Rlist_x = apply(sliding_windows, 1, mean)

    ##
    Rlist_med <-apply(out$R[, ,j], 2, quantile, probs = 0.5)
    Rlist_lb <- apply(out$R[, ,j], 2, quantile, probs = 0.025)
    Rlist_ub <- apply(out$R[, ,j], 2, quantile, probs = 0.975)

    RlistJ[[j]] = data.frame(
      Model = 'STAN',
      J = j,
      x = Rlist_x,
      med = Rlist_med, lb = Rlist_lb, ub = Rlist_ub
    )

    #
    NOBS <- NOBS_mat[, j]
    t_start <- seq(2, length(NOBS)-(tau - 1))
    t_end   <- t_start + (tau - 1)

    if(tau > 1) {
      res <- estimate_R(incid = NOBS,
                        method = "non_parametric_si",
                        config = make_config(list(
                          si_distr = sip,
                          t_start = t_start,
                          t_end = t_end)),
                        backimputation_window = S * 2)
    } else {
      res <- estimate_R(incid = NOBS,
                        method = "non_parametric_si",
                        config = make_config(list(
                          si_distr = sip,
                          t_start = t_start,
                          t_end = t_end)))
    }
    # made it a mid window so its apples-apples
    epiE_x = (res$R$t_end -  res$R$t_start)/2 + res$R$t_start
    epiE_med = res$R$`Median(R)`
    epiE_lb = res$R$`Quantile.0.025(R)`
    epiE_ub = res$R$`Quantile.0.975(R)`
    EpiDf <- data.frame(
      Model = 'EpiEstim',
      J = j,
      x = epiE_x,
      med = epiE_med, lb = epiE_lb, ub = epiE_ub
    )

    RlistJ[[j]] <- rbind(RlistJ[[j]], EpiDf)

  }

  Rlist_df  <- do.call(rbind, RlistJ)
  Rlist_df$J <- factor(Rlist_df$J)

  ##
  xM = Mlist_df
  xR = Rlist_df
  save(xM, xR, P, file =f)
}
