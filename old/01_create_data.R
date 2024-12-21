## code to prepare ALL datasets goes here
library(splines)

#' -----------------------------------------
## shape parameter for serial interval assuming gamma distribution
## rate parameter for serial interval assuming gamma distribution
si_shape <- 2
si_rate  <- 0.5

w <- sapply(1:14, function(x){
  pgamma(x, si_shape, si_rate) - pgamma(x-1, si_shape, si_rate)
})

# this should have a leading 0 right ?
# maybe implemented elsewhere
w <- c(w/sum(w))

plot(w)
#' -----------------------------------------

getR <- function(peak_x, peak_val, max_x, startR = 1, endR = 0.2) {
  # peak_x = 20
  # peak_val = 1.5
  # max_x = 80
  #
  x = 1:max_x
  p1 <- seq(from = startR, to = peak_val, length.out = peak_x)
  p2 <- seq(from = peak_val, to = endR, length.out = max_x - peak_x + 1)
  R <- c(p1, p2[2:length(p2)])

  b2 <- lm(R ~ bs(x, knots = peak_x, degree = 3, intercept = F))
  b2.pred <- predict(b2)
  # plot(x, R)
  # lines(x, b2.pred, col = 'red')

  return(b2.pred)

}

tmax = 80
Rmatrix = cbind(R(20, 2, 80), R(40, 1.8, 80))
head(Rmatrix)
plot(1:tmax, Rmatrix[,1], type = 'l')
lines(1:tmax, Rmatrix[,2], type = 'l', col = 'red')

M <- matrix(NA, nrow = tmax + 1, ncol = 2)
N <- matrix(NA, nrow = tmax + 1, ncol = 2)
R_rev <- matrix(NA, nrow = tmax, ncol = 2)
R_this <- matrix(NA, nrow = tmax, ncol = 2)

#' -----------------------------------------

# some transfer

# P <- matrix(c(0.8, 0.4,
#               0.2, 0.6), 2)


# no transfer

P <- matrix(c(1, 0,
              0, 1), 2)


P
#' -----------------------------------------

# initialize
t = 1
J <- 2
M[1,] <- rep(100, J)
N[1,] <- rep(100, J)
S <- length(w)
add_noise <- T
noise_p = 0.25
set.seed(123)

#' -----------------------------------------
for(t in 1:tmax) {

  # what are the boundaries of serial_interval
  w_i = min(S, t)

  mt = t+1

  # #
  RR <- diag(Rmatrix[t, ])
  RR

  ## MM is m(t - 1), ..., m(1)
  ## where rows are regions
  ## and columns are time
  if(t == 1) {
    MM <- as.matrix(M[1, ])
  } else {
    MM <- t(as.matrix(M[mt - 1:w_i, ]))
  }
  MM
  WW <-  as.matrix(w[1:w_i])
  inner_vec <- RR %*% MM %*% WW
  inner_vec

  outer_vec = t(P) %*% inner_vec

  M[mt, ] = outer_vec

  if(!add_noise) {
    N[mt, ] = sapply(outer_vec, function(x) rpois(1, x))
  } else {
    N[mt, ] = sapply(outer_vec, function(x) rpois(1, runif(1, (1-noise_p)*x,
                                                          (1+noise_p)*x)))
  }

  ## --------------------------------------------------
  ## reverse calc R
  sum_m_w_mat <- matrix(rep(MM %*% WW, times = J), byrow = T, nrow = J)

  # now the a, b, c ... matrix is really just P * rowwise the above
  c_mat <- t(P) * sum_m_w_mat
  c_mat

  R_rev[t, ] <- solve(t(c_mat) %*% c_mat) %*% t(c_mat) %*% M[mt, ]

  R_rev[t, ]
  Rmatrix[t, ]
  DIG = 12
  stopifnot(round(R_rev[t, ], DIG) == round(Rmatrix[t, ], DIG))

  ## --------------------------------------------------
  ## calc this R based on N not MM
  if(t == 1) {
    NN <- as.matrix(N[1, ])
  } else {
    NN <- t(as.matrix(N[mt - 1:w_i, ]))
  }

  sum_m_w_mat <- matrix(rep(NN %*% WW, times = J), byrow = T, nrow = J)

  # now the a, b, c ... matrix is really just P * rowwise the above
  c_mat <- t(P) * sum_m_w_mat
  c_mat

  R_this[t, ] <- solve(t(c_mat) %*% c_mat) %*% t(c_mat) %*% N[mt, ]
  head(R_this)
}

# Now, reset M and N to remove initial cases
M <- M[2:nrow(M),]
N <- N[2:nrow(N),]

plot(M[,1], type = 'l', col = 'red')
lines(N[, 1], col = 'red')
lines(M[,2], type = 'l', col = 'green')
lines(N[, 2], col = 'green')

plot(R_rev[,1], type = 'l', col = 'red', ylim = c(0, 5))
lines(R_this[, 1], col = 'red')
lines(R_rev[,2], type = 'l', col = 'green')
lines(R_this[, 2], col = 'green')
R_this[,2]
# ------------------------------------------------------------------------
## NOW EXPORT


## SAMPLE DATES
ref_date = as.Date('2020-01-01')
sample_dates = seq.Date(from = ref_date, to = ref_date + nrow(N) - 1, by = 'day')
sample_dates

##  SAMPLE MULTI SITE
sample_multi_site = data.frame(date = sample_dates, N[1:(nrow(N)),])
colnames(sample_multi_site)[2:3] <- c('Tatooine', 'Hoth')


## TRANSFER MATRIX
# make daily
P_list <- lapply(1:(nrow(N)), function(x) P)
transfer_matrix <- do.call(rbind, P_list)
dim(transfer_matrix)
colnames(transfer_matrix) <- c('Tatooine', 'Hoth')
rn = vector("character", nrow(N))
for(i in seq(2, (nrow(N)) * 2, by = 2)) {
  rn[i-1] = paste0(sample_dates[i/2], ':Tatooine')
  rn[i]   = paste0(sample_dates[i/2], ':Hoth')
}
rownames(transfer_matrix) <- rn[1:nrow(transfer_matrix)]
head(transfer_matrix)
tail(transfer_matrix)


