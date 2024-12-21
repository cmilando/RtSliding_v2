
source('01_create_data.R')

site_names <- colnames(sample_multi_site)[c(2, 3)]
plot(sample_multi_site$date, y = sample_multi_site$Hoth, col = 'red')
points(sample_multi_site$date, y = sample_multi_site$Tatooine, col = 'blue')

Y <- matrix(integer(1), nrow = nrow(sample_multi_site), ncol = 2)
for(i in 1:nrow(Y)) {
  for(j in c(2, 3)) {
    Y[i,j-1] <- as.integer(sample_multi_site[i,j])
  }
}

all(is.integer(Y))

si <- function (ndays, shape, rate, leading0 = TRUE)
{
  prob <- numeric(ndays)
  for (i in 1:ndays) {
    prob[i] <- pgamma(i, shape = shape, rate = rate) -
      pgamma(i - 1, shape = shape, rate = rate)
  }
  result <- prob/sum(prob)
  if (leading0)
    result <- c(0, result)
  return(result)
}

sip <- si(14, 4.29, 1.18, leading0 = FALSE)

report_dates = sample_multi_site$date
case_matrix = Y

# report dates
stopifnot(all(!is.na(as.Date(report_dates))))
for(i in 2:length(report_dates)) {
  if(as.numeric(report_dates[i] - report_dates[i-1]) > 1) stop()
}

# case days
stopifnot(all(!is.na(case_matrix)))
stopifnot(all(is.integer(case_matrix)))
stopifnot(nrow(case_matrix) == length(report_dates))

# transfer matrix
stopifnot(all(!is.na(transfer_matrix)))
stopifnot(all(is.numeric(transfer_matrix)))
stopifnot(nrow(transfer_matrix) == ncol(case_matrix) * length(report_dates))
stopifnot(ncol(transfer_matrix) == ncol(case_matrix))
for(i in 1:nrow(transfer_matrix)) {
  if(sum(transfer_matrix[i,]) != 1) stop()
}

# need to do validation on serial interval
stopifnot(all(!is.na(sip)))
stopifnot(all(is.numeric(sip)))
stopifnot(sip[1] != 0)

# ------------------------------------------------------------

NN <- as.integer(nrow(case_matrix))

tau = as.integer(7)
## NOTE: interesting, seems to work better for smaller tau?
## that probably is not the expected behavior
## so, it converges for tau = 3, but not for tau = 7


SWrow = as.integer(NN-tau+1)
sliding_windows <- matrix(data = NA, nrow = SWrow, ncol = tau)
window_i = 1
for(NN_i in tau:NN) {
  # NN_i = tau
  start_i = NN_i - tau + 1
  end_i = NN_i
  sliding_windows[window_i, ] = start_i:end_i
  window_i = window_i + 1
}
for(i in 1:nrow(sliding_windows)) {
  for(j in 1:ncol(sliding_windows)) {
    sliding_windows[i,j] <- as.integer(sliding_windows[i, j])
  }
}

sliding_windows

# ------------------------

# int<lower=1> N;          // the number of observations.
# int<lower=1> tau;      // sw_col -- aka sliding window size
# int<lower=1> max_ww;      // sw_row
# int SW[max_ww,tau];     // sliding window matrix
# int<lower=0> Y[N];       // observed cases. **J**
#   int<lower=1> S;          // length of serial interval
# vector[S] W;             // serial interval
# int init_cases;          // initial cases **J**

# Data list for Stan
# also need to do data validation here
stan_data <- list(
  N = 7,                 # number of days
  tau = tau,
  max_ww = SWrow,
  SW = matrix(sliding_windows[1,],nrow =SWrow),
  Y = case_matrix[1:7,1],        # cases
  S = 2,                          # serial interval length
  W = c(0.5, 0.5),                # serial interval vector
  init_cases = c(100)     # initial cases
)

library(rstan)

m_hier <- rstan::stan(file = 'sliding_1d.stan',
                      data = stan_data,
                      verbose = T,
                      iter = 2000,
                      cores = 1,
                      chains = 1)

out <- rstan::extract(m_hier)

dim(out$M)
dim(out$R)

summary(out$R)

summary(out$M)
M













# -------------------------
# Data list for Stan
# also need to do data validation here
stan_data <- list(
  N = NN,                 # number of days
  SWcol = tau,
  SWrow = SWrow,
  SW = sliding_windows,
  J = ncol(case_matrix),  # n regions
  Y = case_matrix,        # cases
  P = transfer_matrix,    # transfer matrix
  S = length(sip),        # serial interval length
  W = sip,                # serial interval vector
  init_cases = c(100, 100)     # initial cases
)

library(rstan)

m_hier <- rstan::stan(file = 'sliding.stan',
                      data = stan_data,
                      verbose = T,
                      iter = 2000,
                      cores = 4,
                      chains = 4)

out <- rstan::extract(m_hier)

### first step -- is it making betas smoother? Yes

xx <- t(out$xbeta[, , 1])
dim(xx)
beta_m <- apply(xx, 1, quantile, probs = 0.5)
xx_l <- apply(xx, 1, quantile, probs = 0.025)
xx_u <- apply(xx, 1, quantile, probs = 0.975)

plot(beta_m)
lines(xx_l)
lines(xx_u)

#### second test -- is that affecting logR? -- seems like no
## why not?

xx <- t(out$logR[, , 1])
Rt <- exp(apply(xx, 1, quantile, probs = 0.5))
xx_l <- exp(apply(xx, 1, quantile, probs = 0.025))
xx_u <- exp(apply(xx, 1, quantile, probs = 0.975))

plot(Rt)
lines(x = 1:80, y = xx_l)
lines(x = 1:80, y = xx_u)

lines(Rmatrix[,1], col = 'green')
lines(R_rev[,1], col = 'red')
lines(R_this[,1], col = 'red')

#### second test -- is that affecting logR? -- seems like no
## why not?

xx <- t(out$Y_out[, , 1])
Y_o <- apply(xx, 1, quantile, probs = 0.5)
xx_l <- apply(xx, 1, quantile, probs = 0.025)
xx_u <- apply(xx, 1, quantile, probs = 0.975)

plot(Y_o)
lines(x = 1:80, y = xx_l)
lines(x = 1:80, y = xx_u)

lines(Rmatrix[,1], col = 'green')
lines(R_rev[,1], col = 'red')
lines(R_this[,1], col = 'red')
