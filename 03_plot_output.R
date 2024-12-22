#' ============================================================================
#' PHASE 3
#' ============================================================================
out <- rstan::extract(m_hier_onY)

# -- OBSERVATIONS
Mlist_med <- apply(out$M, 2, quantile, probs = 0.5)
Mlist_lb <- apply(out$M, 2, quantile, probs = 0.025)
Mlist_ub <- apply(out$M, 2, quantile, probs = 0.975)

stopifnot(length(Mlist_med) == maxt)

# -- RT
# ok for each value of M, its a function of those specific winwos

# -- EpiEstim
t_start <- seq(2, length(NOBS)-(tau - 1))
t_end   <- t_start + (tau - 1)

library(EpiEstim)
res <- estimate_R(incid = NOBS,
                  method = "non_parametric_si",
                  config = make_config(list(
                    si_distr = sip,
                    t_start = t_start,
                    t_end = t_end)))

epiE_x = res$R$t_end
epiE_med = res$R$`Median(R)`
epiE_lb = res$R$`Quantile.0.025(R)`
epiE_ub = res$R$`Quantile.0.975(R)`

#
# REMEBER THAT THE X HERE IS FOR WINDOW,
# So first, get the windows that EpiEstim are using
dim(out$R)
first_w = -1
i = 1
while(first_w < 0) {
  if(res$R[1, 1] == sliding_windows[i, 1] &
     res$R[1, 2] == sliding_windows[i, 2]) {
    first_w = i
  }
  i = i + 1
}
first_w

Rlist_x = sliding_windows[first_w:max_ww, 2]
stopifnot(identical(Rlist_x, epiE_x))

##
Rlist_med <-apply(out$R, 2, quantile, probs = 0.5)[first_w:max_ww]
Rlist_lb <- apply(out$R, 2, quantile, probs = 0.025)[first_w:max_ww]
Rlist_ub <- apply(out$R, 2, quantile, probs = 0.975)[first_w:max_ww]
##

##
Ryrange = range(Rlist_med, Rlist_lb, Rlist_ub,
                RT_calc,
                epiE_med, epiE_lb, epiE_ub)
Ryrange
#' ============================================================================
#' Plots
#' ============================================================================

# --------------------
# Adjust graphics parameters for side-by-side plots
par(mfrow = c(1, 2))  # Two plots in one row

# --------------------
# Plot 1: M(t)
# --------------------
plot(NOBS, col = 'black', cex = 0.8, type = 'p',
     ylab = 'M(t) Values', xlab = 'Day',
     main = paste0('M(t) for tau = ',tau))
polygon(c(1:maxt, rev(1:maxt)), c(Mlist_lb, rev(Mlist_ub)),
        col = rgb(0, 0, 1, 0.2), border = NA)
lines(Mlist_med, type = 'l', col = 'blue',lwd = 1.5)
legend("topright", legend = c("Observed", "STAN"),
       col = c("black", "blue"),
       lty = 1, lwd = 2, cex = 0.8)

# --------------------
# Plot 2: R(t)
# --------------------
plot(x = 1:maxt, y = RT_calc, col = 'black', lwd = 1, type = 'l',
     ylab = 'R(t) Values', xlab = 'Day', lty = 3,
     ylim = Ryrange,
     #ylim = c(0, 2),
     main = paste0('R(t) for tau = ',tau))

plot_a = 0.05

##
lines(x = Rlist_x, y = Rlist_med, type = 'l', col = 'blue')
polygon(c(Rlist_x, rev(Rlist_x)),
        c(Rlist_lb, rev(Rlist_ub)),
        col = rgb(0, 0, 1, plot_a), border = NA)

##
lines(x = epiE_x, y = epiE_med, type = 'l', col = 'green')
polygon(c(epiE_x, rev(epiE_x)),
        c(epiE_lb, rev(epiE_ub)),
        col = rgb(0, 1, 0, plot_a), border = NA)

abline(h = 1, lty = 2)

legend("topright", legend = c("Function", "STAN","EpiEstim"),
       col = c("black", "blue", "green"),
       lty = 1, lwd = 2, cex = 0.8)

# --------------------
# Reset graphical parameters to default
par(mfrow = c(1, 1))
