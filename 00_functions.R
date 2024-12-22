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
  Mout
  return(Rt)
}

#' Get initial R(t) based on tau
#'
#' @param m
#' @param w
#' @param tau
#'
#' @return
#' @export
#'
#' @examples
get_first_Rt <- function(m, w, tau) {

  Rt <- 0

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
