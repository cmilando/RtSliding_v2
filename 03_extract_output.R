#' ============================================================================
#' PHASE 3
#' ============================================================================
out <- rstan::extract(m_hier)

# -- OBSERVATIONS
# ok for each value of M, its a function of those specific winwos
Mlist = vector('list', maxt)
dim(out$M)
for(n in (tau+1):maxt) {
  start_W = max(1,n - tau)
  end_W = min(n, max_ww)
  Mlist[[n]] <- out$M[,n, unique(start_W:end_W)]
}

Mlist_x <- 1:maxt
Mlist_med <- sapply(Mlist, quantile, probs = 0.5)
Mlist_lb <- sapply(Mlist, quantile, probs = 0.025)
Mlist_ub <- sapply(Mlist, quantile, probs = 0.975)

stopifnot(length(Mlist_med) == maxt)

# -- RT
# ok for each value of M, its a function of those specific winwos

# -- EpiEstim
library(EpiEstim)
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
epiE_x = res$R$t_end
epiE_med = res$R$`Median(R)`
epiE_lb = res$R$`Quantile.0.025(R)`
epiE_ub = res$R$`Quantile.0.975(R)`

# -- RT ANALYTICAL
# ok for each value of M, its a function of those specific winwos
analytical_R <- get_analytical_R(NOBS, sliding_windows, sip)

#
# REMEBER THAT THE X HERE IS FOR WINDOW,
# So first, get the windows that EpiEstim are using
dim(out$R)

# match windows
first_w = -1
i = 1
while(first_w < 0) {
  if(res$R[1, 1] == sliding_windows[i, 1] &
     res$R[1, 2] == sliding_windows[i, 2]) {
    first_w = i
  }
  i = i + 1
}
# first_w = 1

Rlist_x = sliding_windows[first_w:max_ww, 2]
stopifnot(identical(Rlist_x, epiE_x))

##
Rlist_med <-apply(out$R, 2, quantile, probs = 0.5)[first_w:max_ww]
Rlist_lb <- apply(out$R, 2, quantile, probs = 0.025)[first_w:max_ww]
Rlist_ub <- apply(out$R, 2, quantile, probs = 0.975)[first_w:max_ww]


##
Ryrange = range(Rlist_med, Rlist_lb, Rlist_ub,
                RT_calc,
                epiE_med, epiE_lb, epiE_ub)
Ryrange

