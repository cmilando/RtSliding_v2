library(EpiEstim)

data("Flu2009")

###
tau = 2
t_start <- seq(2, nrow(Flu2009$incidence)-(tau - 1))
t_end <- t_start + (tau - 1)
head(t_start)
head(t_end)
I = Flu2009$incidence$I
res <- estimate_R(incid = I,
                  method = "non_parametric_si",
                  config = make_config(list(
                    si_distr = Flu2009$si_distr,
                    t_start = t_start,
                    t_end = t_end)))
plot(res)

###3

incid = NOBS
this_si = sip
tau = 1

t_start <- seq(2, length(incid)-(tau - 1))
t_start
t_end <- t_start + (tau - 1)
t_end

head(t_start)
head(t_end)
tail(t_start)
tail(t_end)

method = "non_parametric_si",
config = make_config(list(
  si_distr = this_si,
  t_start = t_start,
  t_end = t_end))


  method <- match.arg(method)
  if (any(dt >= 2)) {
    msg <- "backimputation_window is currently not supported when dt > 1"
    if (backimputation_window > 0)
      stop(msg)
    out <- estimate_R_agg(incid, dt = dt, dt_out = dt_out,
                          iter = iter, config = config, method = method, grid = grid)
    return(out)
  }
  if (backimputation_window) {
    incid <- backimpute_I(incid, window_b = backimputation_window)
  }
  config <- make_config(incid = incid, method = method, config = config)
  config <- process_config(config)
  check_config(config, method)
  if (method == "si_from_data") {
    si_data <- process_si_data(si_data)
    config <- process_config_si_from_data(config, si_data)
    if (!is.null(config$mcmc_control$seed)) {
      cdt <- dic.fit.mcmc(dat = si_data, dist = config$si_parametric_distr,
                          burnin = config$mcmc_control$burnin, n.samples = config$n1 *
                            config$mcmc_control$thin, init.pars = config$mcmc_control$init_pars,
                          seed = config$mcmc_control$seed)
    }
    else {
      cdt <- dic.fit.mcmc(dat = si_data, dist = config$si_parametric_distr,
                          burnin = config$mcmc_control$burnin, n.samples = config$n1 *
                            config$mcmc_control$thin, init.pars = config$mcmc_control$init_pars)
    }
    MCMC_conv <- check_cdt_samples_convergence(cdt@samples)
    c2e <- coarse2estim(cdt, thin = config$mcmc_control$thin)
    cat(paste("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@",
              "\nEstimating the reproduction number for these serial interval",
              "estimates...\n", "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"))
    si_sample <- c2e$si_sample
  }
  if (!is.null(config$seed)) {
    set.seed(config$seed)
  }
  out <- estimate_R_func(incid = incid, method = method, si_sample = si_sample,
                         config = config)
  if (method == "si_from_data") {
    out[["MCMC_converged"]] <- MCMC_conv
  }

  out

  <environment: namespace:EpiEstim>
