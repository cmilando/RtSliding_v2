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
library(lemon)
p1 <- ggplot(Mlist_df) + theme_classic2() +
  facet_rep_wrap(~J,nrow = J) +
  geom_point(aes(x = Mlist_x, y = NOBS, group = J),
             shape = 21) +
  geom_ribbon(aes(x = Mlist_x, ymin = Mlist_lb,
                  ymax = Mlist_ub, group = J),
              alpha = 0.1) +
  geom_line(aes(x = Mlist_x, y = Mlist_med, group = J)) +
  geom_vline(xintercept = c(first_break, second_break),
             color = 'white', linewidth = 1) +
  geom_vline(xintercept = c(first_break, second_break),
           color = 'grey', linewidth = 0.5, linetype = 'dotted') +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  ggtitle(paste("Cases and M(t) with τ =", tau)) +
  ylab('# cases') + xlab("Day")

# -- RT
# ok for each value of M, its a function of those specific winwos
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

  if(tau == 1) {
    RlistJ[[j]] <- RlistJ[[j]][2:nrow(RlistJ[[j]]), ]
    warning('first window omitted since tau == 1')
  }

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
library(ggpubr)
library(tidyverse)
p2 <- ggplot(Rlist_df) + theme_classic2() +
  facet_rep_wrap(~J,nrow = J) +
  coord_cartesian(xlim = range(Mlist_x)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_ribbon(aes(x = x, ymin = lb,
                  ymax = ub, fill = Model, group = Model),
              alpha = 0.2) +
  geom_line(aes(x = x, y = med, color = Model, group = Model)) +
  geom_vline(xintercept = c(first_break, second_break),
             color = 'white', linewidth = 1) +
  geom_vline(xintercept = c(first_break, second_break),
             color = 'grey', linewidth = 0.5, linetype = 'dotted') +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  ggtitle(paste("R(t) with τ =", tau)) + ylab('R(t)') + xlab("Day")

library(patchwork)
print(p1 + p2 + plot_layout(nrow = 1))
ggsave('plot.png', units = 'in',
            width = 9, height = 6)
