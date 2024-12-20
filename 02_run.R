
source('01_create_data.R')

site_names <- colnames(sample_multi_site)[c(2, 3)]
plot(sample_multi_site$date, y = sample_multi_site$Tatooine, col = 'red')
points(sample_multi_site$date, y = sample_multi_site$Hoth, col = 'blue')

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

tau = as.integer(3)
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
  init_cases = case_matrix[1, ]     # initial cases
)


initf1 <- function() {
  #
  # vector<lower=0.00005>[J] xsigma;  // region-specific st-dev
  # matrix[NW,J] xbeta;         // time-region specific beta
  # matrix[NW,J] logR;          // time-region specific R values, in Log space
  #
  list(xsigma = rep(0.1, ncol(case_matrix)),
       xbeta = matrix(0, nrow = SWrow, ncol = ncol(case_matrix)),
       logR = matrix(0.1, nrow = NN, ncol = ncol(case_matrix)))
}

library(rstan)

m_hier <- rstan::stan(file = 'sliding.stan',
                      data = stan_data,
                      init = initf1,
                      verbose = T,
                      cores = 4,
                      chains = 4)

out <- rstan::extract(m_hier)

dim(out$M)
dim(out$xsigma)

data_l <- lapply(1:dim(out$M)[3], function(i) {
  data.frame(
    x      = 1:dim(out$M)[2],
    y_real = Y[, i],
    y      = apply(out$M[, , i], 2, mean),
    yl     = apply(out$M[, , i], 2, quantile, probs = 0.025),
    yh     = apply(out$M[, , i], 2, quantile, probs = 0.975),
    Rt     = apply(out$R[, , i], 2, mean),
    Rtl    = apply(out$R[, , i], 2, quantile, probs = 0.025),
    Rth    = apply(out$R[, , i], 2, quantile, probs = 0.975),
    region = site_names[i] # must be a string
  )
})

data_all <- do.call(rbind, data_l)
head(data_all)

# Summarise data
data_all_summarise <- aggregate(
  cbind(y, y_real, yl, yh, Rt, Rtl, Rth) ~ x + region,
  data = data_all,
  FUN = mean
)

head(data_all_summarise)

# Generate a color palette
regions <- unique(data_all_summarise$region)
# Define a set of distinct colors
colors <- c("red", "blue", "green", "purple", "orange", "brown",
            "pink", "yellow", "cyan", "magenta")
colors <- colors[1:length(regions)]
names(colors) <- regions

# Plot expected cases
plot(
  x = as.integer(data_all_summarise$x),
  y = as.numeric(data_all_summarise$y),
  type = "n",
  xlab = "Days",
  ylab = "Cases",
  main = "Expected Cases"
)

for (region_i in regions) {
  region_data <- subset(data_all_summarise, region == region_i)
  points(x = region_data$x, region_data$y_real, col = colors[region_i])
  polygon(
    c(region_data$x, rev(region_data$x)),
    c(region_data$yl, rev(region_data$yh)),
    col = adjustcolor(colors[region_i], alpha.f = 0.3), border = NA
  )
  lines(region_data$x, region_data$y, col = colors[region_i], lwd = 0.5)
}


lines(M[,1], type = 'l', col = 'blue')
lines(N[, 1], col = 'blue')
lines(M[,2], type = 'l', col = 'red')
lines(N[, 2], col = 'red')

legend("topright",
       legend = regions,
       col = colors,
       lty = rep(1, length(regions)),
       cex = 0.8,
       pt.cex = 1.5 ) # Text size

# Plot R(t)
plot(
  data_all_summarise$x, data_all_summarise$Rt,
  xlab = "Days", ylab = "Reproduction Number",
  type = "n",
  main = "R(t)",
  ylim = c(0, 5)
)
for (region_i in regions) {
  region_data <- subset(data_all_summarise, region == region_i)
  polygon(
    c(region_data$x, rev(region_data$x)),
    c(region_data$Rtl, rev(region_data$Rth)),
    col = adjustcolor(colors[region_i], alpha.f = 0.3), border = NA
  )
  lines(region_data$x, region_data$Rt, col = colors[region_i], lwd = 0.5)
}
abline(h = 1, col = "black", lwd = 1, lty = 1)

lines(R_rev[,1], type = 'l', col = 'blue')
lines(R_this[, 1], col = 'blue')
lines(R_rev[,2], type = 'l', col = 'red')
lines(R_this[, 2], col = 'red')

legend("topright",
       legend = regions,
       col = colors,
       lty = c(1, 1),
       cex = 0.8,
       pt.cex = 1.5 ) # Text size
