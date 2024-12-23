#' Title
#'
#' @param peak_x
#' @param peak_val
#' @param max_x
#' @param startR
#' @param endR
#'
#' @return
#' @export
#'
#' @examples
getSPLINE <- function(peak_x, peak_val,
                      max_x, startR = 1, endR = 0.4) {

  require(splines)
  x = 1:max_x
  p1 <- seq(from = startR, to = peak_val, length.out = peak_x)
  p2 <- seq(from = peak_val, to = endR, length.out = max_x - peak_x + 1)
  R <- c(p1, p2[2:length(p2)])

  b2 <- lm(R ~ bs(x, knots = peak_x, degree = 3, intercept = F))
  b2.pred <- predict(b2)

  return(b2.pred)

}


#' Get expected value of cases from Rt, init_cases, and w
#'
#' @param Rt
#' @param init_cases
#' @param w
#'
#' @return
#' @export
#'
#' @examples
get_Mt <- function(Rt, init_cases, w) {

  Mout <- rep(0, length(Rt))

  for(t in 1:maxt) {

    inner_vec = 0

    S_loop_max = min(length(w), t + 1)
    forward_vec = 1:S_loop_max
    rev_vec = t + 1 - forward_vec

    for(si in 1:S_loop_max) {
      rev_i = rev_vec[si]

      if(rev_i < 0) {
        mx = 0
      } else if(rev_i == 0) {
        mx = init_cases
      } else {
        mx = Mout[rev_i]
      }

      inner_vec = inner_vec + w[si] * mx
    }
    Mout[t] = Rt[t] * inner_vec

  }

  stopifnot(!any(is.na(Mout)))
  stopifnot(all((Mout != 0)))

  return(Mout)
}


#' Get Rt from expected value of cases, init_cases, and w
#'
#' @param m
#' @param init_cases
#' @param w
#'
#' @return
#' @export
#'
#' @examples
get_Rt <- function(m, init_cases, w) {

  Rt <- rep(0, length(m))

  for(t in 1:maxt) {

    inner_vec = 0

    S_loop_max = min(length(w), t + 1)
    forward_vec = 1:S_loop_max
    rev_vec = t + 1 - forward_vec

    for(si in 1:S_loop_max) {
      rev_i = rev_vec[si]

      if(rev_i < 0) {
        mx = 0
      } else if(rev_i == 0) {
        mx = init_cases
      } else {
        mx = m[rev_i]
      }

      inner_vec = inner_vec + w[si] * mx
    }
    Rt[t]  = m[t] / inner_vec

  }
  return(Rt)
}


#' Get sliding window matrix based on maxt and tau
#'
#' @param maxt
#' @param tau
#'
#' @return
#' @export
#'
#' @examples
get_SW <- function(maxt, tau) {
  max_ww     = maxt - tau
  sliding_windows <- matrix(data = NA, nrow = max_ww, ncol = 2)
  window_i = 1
  startN = 2
  endN = startN + (tau-1)
  sliding_windows[window_i, 1] = startN
  sliding_windows[window_i, 2] = endN

  for(window_i in 2:max_ww) {
    startN = startN + 1
    endN = endN + 1
    stopifnot(endN <= maxt)
    sliding_windows[window_i, 1] = startN
    sliding_windows[window_i, 2] = endN
  }

  for(i in 1:nrow(sliding_windows)) {
    for(j in 1:ncol(sliding_windows)) {
      sliding_windows[i,j] <- as.integer(sliding_windows[i, j])
    }
  }
  return(sliding_windows)

}


#' Get analytical R, the guts of EpiEstim
#'
#' @param incid
#' @param sliding_windows
#' @param sip
#' @param mean_prior
#' @param std_prior
#'
#' @return
#' @export
#'
#' @examples
get_analytical_R <- function(incid, sliding_windows, sip,
                             mean_prior = 5, std_prior = 5) {

  t_start <- sliding_windows[, 1]
  t_end <- sliding_windows[, 2]
  si_distr <- sip

  ## SINCE R has a GAMMA PRIOR, you can skip to this
  a_prior <- (mean_prior / std_prior)^2
  b_prior <- std_prior^2 / mean_prior

  #########################################################
  overall_infectivity <- function(incid, si_distr) {
    maxt <- length(incid)
    lambda <- vector()
    lambda[1] <- NA
    for (t in seq(2, maxt)) {
      lambda[t] <- sum(si_distr[seq_len(t)] * incid[seq(t, 1)], na.rm = TRUE)
    }
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

  ## SUM INCID, DOESN'T CHANGE WITH M[1]
  a_posterior <-
    vnapply(seq_len(nb_time_periods), function(t) {
      a_prior + sum(incid[seq(t_start[t], t_end[t])])
    })

  ## SUM LAMBDA, DOES CHANGE WITH M[1]
  b_posterior <-
    vnapply(seq_len(nb_time_periods), function(t) {
      1 / (1 / b_prior + sum(lambda[seq(t_start[t], t_end[t])]))
    })

  ## output
  get_posterior <- function(q) {
    qgamma(q,shape = a_posterior,
           scale = b_posterior,
           lower.tail = TRUE, log.p = FALSE)
  }

  out_df <- data.frame(
    t_start, t_end, lb = get_posterior(0.025),
    med = get_posterior(0.5),
    ub = get_posterior(0.975)
  )

  return(out_df)
}


#' Get Gamma SI
#'
#' @param ndays
#' @param shape
#' @param rate
#' @param leading0
#'
#' @return
#' @export
#'
#' @examples
si <- function(ndays, shape, rate, leading0 = TRUE) {

  prob <- numeric(ndays)

  for (i in 1:ndays){
    prob[i] <- pgamma(i, shape = shape, rate = rate) -
      pgamma(i - 1, shape = shape, rate = rate)
  }

  result <- prob/sum(prob)

  # Add a leading 0 to indicate no cases
  if(leading0) result <- c(0, result)

  return(result)
}
