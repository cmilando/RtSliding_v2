# /// TRANSFORMED PARAMS
Rtau = vector('numeric', max_ww)

# -- ww loop
for(ww in 1:max_ww) {

  xbeta = rnorm(1, 0, 1)
  logR = rnorm(1, BETAmatrix[ww], xsigma)
  Rtau[ww] = exp(logR)

  startN = sliding_windows[ww, 1]
  endN   = sliding_windows[ww, tau]
  xwindow = startN:endN

  # // renewal equation forwards
  # // but this is why you have to shoot all the values across
  # // 𝑚(𝑡)= 𝑅(𝑡) Σ𝜏<𝑡 𝜔(𝜏) 𝑚(𝑡−𝜏)
  for(m_i in startN:endN) {
    inner_vec = 0
    sip_i = min(S, startN)
    for(si in 1:sip_i) {
      inner_vec = inner_vec + sip[si] * M[m_i + 1 - si, ww];
    }
    M[m_i + 1, ww] = Rtau[ww] * inner_vec
  }

  # carry across
  if((startN + 1) < endN) {
    for(wremain in min(ww + 1, max_ww):max_ww) {
      M[startN + 1, wremain] = M[startN + 1,ww]
    }
  }

}

plot(Rtau)
lines(exp(BETAmatrix))

# /// MODEL BLOCK
Mlist = vector('list', maxt)
# -- get just the
# ok for each value of M, its a function of those specific winwos
for(n in 1:maxt) {
  start_W = max(1,n - tau + 1)
  end_W = min(n, max_ww)
  mn = n + 1
  Mlist[[n]] <- M[mn, start_W:end_W]
}

Mlist_med <- sapply(Mlist, quantile, probs = 0.5)
Mlist_lb <- sapply(Mlist, quantile, probs = 0.025)
Mlist_ub <- sapply(Mlist, quantile, probs = 0.975)

plot(Mlist_med, type = 'l')
lines(Mlist_lb, type = 'l')
lines(Mlist_ub, type = 'l')
