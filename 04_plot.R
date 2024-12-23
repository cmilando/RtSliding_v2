#' ============================================================================
#' Plots
#' ============================================================================

# Adjust graphics parameters for side-by-side plots

png("plot.png", units = 'in',
    width = 9, height = 7.5, res = 600) # 300 DPI
## set graph params
par(mfrow = c(2, 2))  # Two plots in one row

###
plot(MTRUE, type = 'l', ylab = 'Cases', xlab = 'Day',
     col = 'red', lwd = 1.5,
     ylim = range(MTRUE, NOBS),
     main = paste0('Cases'))
points(NOBS, cex = 0.8)

legend("topright", legend = c("Function", "Observed"),
       col = c("red", "black"),
       lty = 1, lwd = 2, cex = 0.8)

###
plot(RT_calc, col = 'black', type = 'l', ylab = 'R(t)', xlab = 'Day',
     ylim = range(RT_calc, R_TRUE_matrix),
     main = paste0('R(t)'), lty = 3)
lines(R_TRUE_matrix, type = 'l', col = 'red')
abline(h = 1, lty = 2)
legend("topright", legend = c("Function", "Observed"),
       col = c("red", "black"),
       lty = 1, lwd = 2, cex = 0.8)


# Plot 1: M(t)
plot(NOBS, col = 'black', cex = 0.8, type = 'p',
     ylab = 'M(t) Values', xlab = 'Day',
     main = paste0('M(t) for tau = ',tau))
polygon(c(1:maxt, rev(1:maxt)), c(Mlist_lb, rev(Mlist_ub)),
        col = rgb(0, 0, 1, 0.2), border = NA)
lines(Mlist_med, type = 'l', col = 'blue',lwd = 1.5)
legend("topright", legend = c("Observed", "STAN"),
       col = c("black", "blue"),
       lty = 1, lwd = 2, cex = 0.8)

# Plot 2: R(t)
plot(x = 1:maxt, y = RT_calc, col = 'black', lwd = 1, type = 'l',
     ylab = 'R(t) Values', xlab = 'Day', lty = 3,
     ylim = Ryrange, # ylim = c(0, 2),
     main = paste0('R(t) for tau = ',tau))

head(Rlist_med)
length(Rlist_med)
length(epiE_med)

lines(x = Rlist_x, y = Rlist_med, type = 'l', col = 'blue')
polygon(c(Rlist_x, rev(Rlist_x)),
        c(Rlist_lb, rev(Rlist_ub)),
        col = rgb(0, 0, 1, 0.2), border = NA)

lines(x = analytical_R$t_end, lty = 2,
      y = analytical_R$med, type = 'l', col = 'red')
polygon(c(analytical_R$t_end, rev(analytical_R$t_end)),
        c(analytical_R$lb, rev(analytical_R$ub)),
        col = rgb(1, 0, 0, 0.2), border = NA)


lines(x = epiE_x, y = epiE_med, type = 'l', col = 'green')
polygon(c(epiE_x, rev(epiE_x)),
        c(epiE_lb, rev(epiE_ub)),
        col = rgb(0, 1, 0, 0.2), border = NA)

abline(h = 1, lty = 2)

legend("topright", legend = c("Function", "STAN","EpiEstim",
                              "AnalyticalR"),
       col = c("black", "blue", "green", "red"),
       lty = 1, lwd = 2, cex = 0.8)



# Reset graphical parameters to default
par(mfrow = c(1, 1))

dev.off()
