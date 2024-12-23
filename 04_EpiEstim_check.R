


library(EpiEstim)
incid = NOBS
incid[1] <- NOBS[1]
#incid[1] <- 1000
method = "non_parametric_si"
t_start <- seq(2, length(NOBS)-(tau - 1))
t_end   <- t_start + (tau - 1)
si_distr <- sip

res <- estimate_R(incid = incid,
                  method = "non_parametric_si",
                  config = make_config(list(
                    si_distr = sip,
                    t_start = t_start,
                    t_end = t_end)))


##
mean_prior = 5
std_prior = 5
a_prior <- (mean_prior / std_prior)^2
b_prior <- std_prior^2 / mean_prior
a_prior
b_prior

#########################################################
overall_infectivity <- function(incid, si_distr) {

  maxt <- length(incid)
  lambda <- vector()
  lambda[1] <- NA
  for (t in seq(2, maxt))
    lambda[t] <- sum(si_distr[seq_len(t)] * incid[seq(t, 1)], na.rm = TRUE)
  return(lambda)
}

vnapply <- function(X, FUN, ...) {
  vapply(X, FUN, numeric(1), ...)
}

#########################################################

nb_time_periods <- length(t_start)
seq_len(nb_time_periods)

## <<< THIS IS A FUNCTION THAT INCLUDES M[0] <<<<<
lambda <- overall_infectivity(incid, si_distr)

a_posterior <-
  vnapply(seq_len(nb_time_periods), function(t)
      a_prior + sum(incid[seq(t_start[t], t_end[t])])  ## SUM INCID, DOESN'T CHANGE WITH M[1]
)

b_posterior <-
  vnapply(seq_len(nb_time_periods), function(t)
      1 / (1 / b_prior + sum(lambda[seq(t_start[t], t_end[t])])) ## SUM LAMBDA, DOES CHANGE WITH M[1]
)

set.seed(123)
median_posterior <- qgamma(0.5,
                           shape = a_posterior,
                           scale = b_posterior, lower.tail = TRUE, log.p = FALSE
)
median_posterior
res$R$`Median(R)`
