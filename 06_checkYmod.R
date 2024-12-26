

RlistJ = vector("list", J)

for(j in 1:J) {

  Rlist_x = 1:(length(report_delay))
  #stopifnot(identical(Rlist_x, epiE_x))

  ##
  Rlist_med <-apply(ar[, ,j], 2, quantile, probs = 0.5)
  Rlist_lb <- apply(ar[, ,j], 2, quantile, probs = 0.025)
  Rlist_ub <- apply(ar[, ,j], 2, quantile, probs = 0.975)

  RlistJ[[j]] = data.frame(
    J = j,
    Rlist_x, Rlist_med, Rlist_lb, Rlist_ub
  )

}

Rlist_df  <- do.call(rbind, RlistJ)
Rlist_df$J <- factor(Rlist_df$J)
library(ggpubr)
library(tidyverse)
ggplot(Rlist_df) + theme_classic2() +
  geom_ribbon(aes(x = Rlist_x, ymin = Rlist_lb,
                  ymax = Rlist_ub, fill = J, group = J),
              alpha = 0.1) +
  geom_point(aes(x = Rlist_x, y = Rlist_med, color = J, group = J)) +
geom_line(aes(x = Rlist_x, y = Rlist_med, color = J, group = J))

