#' ============================================================================
#' PHASE 3
#' ============================================================================
out <- rstan::extract(m_hier)

# -- OBSERVATIONS
# ok for each value of M, its a function of those specific winwos
# don't oversample, just take the last row
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
library(ggpubr)
library(tidyverse)
p1 <- ggplot(Mlist_df) + theme_classic2() +
  geom_point(aes(x = Mlist_x, y = NOBS, color = J, group = J),
             shape = 21) +
  geom_ribbon(aes(x = Mlist_x, ymin = Mlist_lb,
                  ymax = Mlist_ub, fill = J, group = J),
              alpha = 0.1) +
  geom_line(aes(x = Mlist_x, y = Mlist_med, color = J, group = J)) +
  geom_vline(xintercept = c(first_break, second_break),
             color = 'white', linewidth = 1) +
  geom_vline(xintercept = c(first_break, second_break),
           color = 'grey', linewidth = 0.5, linetype = 'dotted')

# -- RT
# ok for each value of M, its a function of those specific winwos
RlistJ = vector("list", J)
for(j in 1:J) {

  Rlist_x = apply(sliding_windows, 1, mean)

  ##
  Rlist_med <-apply(out$R[, ,j], 2, quantile, probs = 0.5)
  Rlist_lb <- apply(out$R[, ,j], 2, quantile, probs = 0.025)
  Rlist_ub <- apply(out$R[, ,j], 2, quantile, probs = 0.975)

  RlistJ[[j]] = data.frame(
    J = j,
    Rlist_x,
    Rlist_med, Rlist_lb, Rlist_ub
  )

}

Rlist_df  <- do.call(rbind, RlistJ)
Rlist_df$J <- factor(Rlist_df$J)
library(ggpubr)
library(tidyverse)
p2 <- ggplot(Rlist_df) + theme_classic2() +
  coord_cartesian(xlim = range(Mlist_x)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_ribbon(aes(x = Rlist_x, ymin = Rlist_lb,
                  ymax = Rlist_ub, fill = J, group = J),
              alpha = 0.1) +
  geom_line(aes(x = Rlist_x, y = Rlist_med, color = J, group = J)) +
  geom_vline(xintercept = c(first_break, second_break),
             color = 'white', linewidth = 1) +
  geom_vline(xintercept = c(first_break, second_break),
             color = 'grey', linewidth = 0.5, linetype = 'dotted')

library(patchwork)
print(p1 + p2 + plot_layout(nrow = 1))
