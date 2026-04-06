source("packages.R") # attach necessary packages

#' *State space models*:
#' - process model:                    `mu_proc = b0 + b1 * x1 + ...`
#' - process output (states):          `Y_proc ~ MVN(mu_proc, s_proc)`
#' - process observations at time `t`: `O_t ~ MVN(Y_proc, s_obs)`
#' 
#' but estimating the process model requires us to work backwards:
#' - model the space of possible states (i.e., true outcomes, true responses)
#' - process observations at time `t`: `O_t ~ MVN(Y_proc, s_obs)`
#' - process output (states):          `Y_proc ~ MVN(mu_proc, s_proc)`
#' - process model:                    `mu_proc = b0 + b1 * x1 + ...`
#' uncertainty needs to be propagated accordingly across each step

# Forecasting from dynamic models ----

#' applications of state-space models:
#' Kalman filter and Apollo missions: `https://doi.org/10.1109/MCS.2010.936465`

air_passengers <-
  tibble(time = 1:length(AirPassengers),
         dec_date = time(AirPassengers) %>% as.numeric(),
         year = floor(dec_date),
         month = round((dec_date - year) * 12) + 1, # to ensure January is 1
         passengers = as.numeric(AirPassengers)) # in thousands

air_passengers <- mutate(air_passengers, lag_12_passengers = lag(passengers, 12))
data_train <- filter(air_passengers, year <= 1955) %>%
  filter(! is.na(lag_12_passengers)) # lagged value before 1st observation is NA
data_test <- filter(air_passengers, year > 1955)

m_gam <- mvgam(formula = passengers ~ 0, # no error in observation process
               trend_formula = ~
                 log(lag_12_passengers) + # since we are on the log link scale
                 s(year, k = 5, bs = "tp") +
                 s(month, k = 10, bs = "cc"),
               trend_model = "None",
               noncentred = TRUE,
               knots = list(month = c(0.5, 12.5)),
               family = poisson(link = "log"),
               data = data_train,
               newdata = data_test, # calculate forecast while fitting
               chains = 4,
               burnin = 750,
               samples = 500,
               control = list(max_treedepth = 20, adapt_delta = 0.9),
               parallel = TRUE,
               silent = 2)

summary(m_gam) # check diagnostics
plot(m_gam)

# fit a GAM with an AR(1) term
# not fixing the AR(1) coefficient causes issues with PSIS diagnostics
get_mvgam_priors(formula = passengers ~ 0, # no error in observation process
                 trend_formula = ~
                   log(lag_12_passengers) +
                   s(year, k = 5, bs = "tp") +
                   s(month, k = 10, bs = "cc"),
                 trend_model = AR(p = 1),
                 data = data_train)

m_gam_ar <- mvgam(formula = passengers ~ 0, # no error in observation process
                  trend_formula = ~
                    log(lag_12_passengers) +
                    s(year, k = 5, bs = "tp") +
                    s(month, k = 10, bs = "cc"),
                  trend_model = AR(p = 1), # AR(1) model
                  # adding bounds for AR(1) range to improve diagnostics
                  priors = prior(normal(0.4, 0.01), class = ar1,
                                 lb = 0.35, ub = 0.45),
                  noncentred = TRUE, # use a noncentered AR(1) model
                  knots = list(month = c(0.5, 12.5)),
                  family = poisson(link = "log"),
                  data = data_train,
                  newdata = data_test, # calculate forecast while fitting
                  chains = 4,
                  burnin = 750,
                  samples = 1000,
                  control = list(max_treedepth = 20, adapt_delta = 0.95),
                  parallel = TRUE, silent = 2)

plot(m_gam_ar)
summary(m_gam_ar) # check diagnostics

# plot diagnostics
# one point is slightly problematic
# not restricting the range of the AR(1) coef results in Pareto shapes > 5 due
# to very low effective sample sizes as the AR(1) fights with the smooth terms
layout(matrix(c(1, 1:3), ncol = 2, byrow = TRUE))
plot(m_gam_ar, type = "forecast")
plot(loo(m_gam_ar), diagnostic = "k")
abline(v = 50, col = "grey", lty = "dashed")
plot(loo(m_gam_ar), diagnostic = "ESS") #' same as `diagnostic = "n_eff"`
layout(1)

#' Dynamic `mvgam` models contain draws for many quantities, all stored as MCMC
#' draws in an object of class `stanfit` in the `model_output` slot:
#' - `β` coefficients for linear predictor terms (called `b`)
#' - Family-specific shape/scale parameters:
#'    - `ϕ` for Negative Binomial,
#'    - `σ_obs` for Normal / LogNormal
#' - Trend-specific parameters:
#'    - `α` and `ρ` for GP trends,
#'    - `σ` and `ar1` for AR trends
#' - In-sample posterior predictions: `ypred`
#' - In-sample posterior trend estimates: `trend`

class(m_gam_ar$model_output)
m_gam_ar$model_output@model_pars # names of model parameters
m_gam_ar$model_output@par_dims # dimensions of parameter vectors/matrices
m_gam_ar$model_output # summary table of parameter samples
m_gam_ar$model_output@sim$samples[[1]][1:10, 1:8] # df of samples for chain 1

# view posterior draws of the trend
plot(m_gam_ar, type = "forecast") # with base plot
plot(forecast(m_gam_ar)) # with ggplot2
plot(forecast(m_gam_ar), realisations = TRUE) # CIs = summaries of realizations

# random draws from the posterior (NOTE: x axis is time since first observation)
plot(m_gam_ar, type = "forecast", realisations = TRUE, n_realisations = 10)
plot(m_gam_ar, type = "trend", realisations = TRUE, n_realisations = 10) +
  geom_vline(xintercept = nrow(data_train), lty = "dashed")
plot(m_gam_ar, type = "smooths", realisations = TRUE, n_realisations = 10,
     trend_effects = TRUE)

# draws summarized to credible intervals (NOTE: x axis is time since first obs)
plot(m_gam_ar, type = "forecast")
plot(m_gam_ar, type = "trend") +
  geom_vline(xintercept = nrow(data_train), lty = "dashed")
plot(m_gam_ar, type = "smooths", trend_effects = TRUE)

# generate forecasts for up to and of 2026
# predicting later is useful if data are not available or too large to add
data_test
preds_2026 <- tibble(time = max(data_train$time) + 1:(12 * (2027 - 1963)),
                     dec_date = max(data_train$dec_date) + time/12,
                     year = floor(dec_date),
                     month = round((dec_date - year) * 12) + 1) %>%
  left_join(air_passengers %>% select(time, passengers, lag_12_passengers),
            by = "time")

#' need to predict with `for` loop since lag-12 values are not always avaiable
#' not the best way to predict: it does not include uncertainty in lagged values
for(i in which(is.na(preds_2026$passengers))) {
  if(is.na(preds_2026$lag_12_passengers[i])) {
    preds_2026$lag_12_passengers[i] <- preds_2026$passengers[i - 12]
  }
  
  preds_2026$passengers[i] <-
    predict(m_gam_ar, preds_2026[i, ], type = "expected")[, "Estimate"]
}

# to stop from calculating score (since values are predicted)
preds_2026 <- rename(preds_2026, estimate = passengers)

# model predicts that passengers will exceed 8 million by the end of 1972
print(preds_2026, n = 13) #' `lag_12_passengers` for row 13 is `estimate` for 1
max(preds_2026$dec_date)
plot_mvgam_fc(m_gam_ar, newdata = preds_2026, ylim = c(0, 8e6 / 1e3),
              realisations = TRUE)
plot_mvgam_fc(m_gam_ar, newdata = preds_2026, ylim = c(0, 8e6 / 1e3))
abline(v = 200, lty = "dashed")
filter(preds_2026, time == 200)

#' `plot(forecast(m_gam_ar, newdata = new_data))` fails with error:
#' `arguments imply differing number of rows: 792, 852`
#' function is adding the test data twice: `nrow(data_test)` is 60

#' **break**

# interpreting predictions ----
summary(m_gam_ar)

# - coefficients are often hard to interpret for GAMs, especially if non-gaussian
# - p-values are often too small for smooth terms bc they ignore uncertainty in
#   the smoothness parameter
# - predictions are more interpretable than coefficients

plot_preds <- function(.d, scale){
  .d %>%
    as.data.frame() %>%
    mutate(time = 1:n()) %>%
    ggplot() +
    geom_ribbon(aes(time, ymin = Q2.5, ymax = Q97.5), alpha = 0.3) +
    geom_line(aes(time, Estimate)) +
    labs(x = "Time", y = paste0("Estiamted effect (", scale, " scale)"))
}

# on link scale: for understanding coefficients of the linear predictor
predict(m_gam_ar, type = "link") %>% #' = `brms::posterior_linpred()`
  plot_preds(scale = "link") +
  ylim(c(0, 410))

# on expected scale: for understanding effects on the mean response
predict(m_gam_ar, type = "expected") %>% #' = `brms::posterior_epred()`
  plot_preds(scale = "expected") +
  ylim(c(0, 410))

# on response scale: for understanding effects on the individual observations
predict(m_gam_ar, type = "response") %>% #' = `brms::posterior_predict()`
  plot_preds(scale = "response") +
  ylim(c(0, 410))

#' can include process error with `process_error = TRUE`
predict(m_gam_ar, type = "response", process_error = TRUE) %>%
  plot_preds(scale = "response") +
  ylim(c(0, 410))

# example with count data:
# - link: values are all real numbers (+ or -); scale is additive
# - expected: values are > 0 (including decimals); scale is multiplicative
# - response: values are > 0 (including decimals); scale is multiplicative
# 
# but note that:
# - the mean can be any real number > 0 (including decimals)
# - the predicted response values can only be integers > 0

#' link-scale partial effects are centered around 0
#' allows to add intercept term to the partial effects
#' intercept = est. mean response, averaged across all smooths, on link scale
draw(m_gam_ar$trend_mgcv_model)
plot(forecast(m_gam_ar, type = "link"))
head(as.vector(forecast(m_gam_ar, type = "link")$forecasts$series1))

#' if `link = "log"` (does not apply if `link = "logit"`):
#' expected-scale partial effects are centered around 1
#' show the relative change in response with the predictor
#' allows to multiply intercept term by the partial effects
#' intercept = est. mean response, averaged across all smooths, on resp. scale
draw(m_gam_ar$trend_mgcv_model, fun = exp) # relative change in passengers
draw(m_gam_ar$trend_mgcv_model, # partial change in passengers
     fun = \(x) exp(coef(m_gam_ar$trend_mgcv_model)["(Intercept)"]) * exp(x))
plot(forecast(m_gam_ar, type = "expected"))
head(as.vector(forecast(m_gam_ar, type = "expected")$forecasts$series1))

#' *NOTE:* in `{mgcv}` and `{gratia}`, predictions on the "response" scale are
#'         actually on the expected scale, since they only include uncertainty
#'         in the mean
#'         in `{mvgam}`, predictions on the response scale include uncertainty
#'         at the observation level, rather than just the mean, and the values
#'         are always integers for count models
plot(forecast(m_gam_ar, type = "response"))
head(as.vector(forecast(m_gam_ar, type = "response")$forecasts$series1))

# useful for checking if the model has a good fit: can it simulate data well?
pp_check(m_gam_ar, type = "ribbon", ndraws = 100)
pp_check(m_gam_ar, type = "intervals", ndraws = 100)
pp_check(m_gam_ar, type = "scatter", ndraws = 9)
pp_check(m_gam_ar, type = "scatter_avg", ndraws = 100)
pp_check(m_gam_ar, type = "hist", ndraws = 8, bins = 10)
pp_check(m_gam_ar, type = "dens_overlay", ndraws = 100)
pp_check(m_gam_ar, type = "ecdf_overlay", ndraws = 100)

pp_check(m_gam_ar, type = "stat", ndraws = 10, stat = "mean", binwidth = 1)
pp_check(m_gam_ar, type = "stat", ndraws = 10, stat = "median", binwidth = 2.5)
pp_check(m_gam_ar, type = "stat", ndraws = 10, stat = "sd", binwidth = 2.5)
pp_check(m_gam_ar, type = "stat", ndraws = 10, stat = "var", binwidth = 250)

pp_check(m_gam_ar, type = "stat_2d", ndraws = 1000, stat = c("mean", "sd"))

#' **break**

# comparing models ----
layout(1:2)
acf(air_passengers$passengers)
pacf(air_passengers$passengers)
layout(1)

m_bad <- mvgam(formula = passengers ~ 0,
               trend_formula = ~ 1,
               trend_model = AR(p = 1), # AR(1) model
               noncentred = TRUE, # use a noncentered AR(1) model
               knots = list(month = c(0.5, 12.5)),
               family = poisson(link = "log"),
               data = data_train, # calculate forecast while fitting
               newdata = data_test,
               chains = 4,
               burnin = 750,
               samples = 500,
               parallel = TRUE,
               silent = 2)

layout(matrix(c(1, 1:3), ncol = 2, byrow = TRUE))
plot(m_bad, type = "forecast") # the forecast is an AR(1) around the mean
plot(loo(m_bad), diagnostic = "k")
plot(loo(m_bad), diagnostic = "ESS") #' same as `diagnostic = "n_eff"`
# points closer to the estimated mean have better scores
layout(1)

# both models predict decently well for past data
plot(hindcast(m_gam_ar)) / plot(hindcast(m_bad))

# but they do not model the data the same way
plot_predictions(m_gam_ar, by = "time") / # GAM model "understands" the trends
  plot_predictions(m_bad, by = "time") # the AR model always assumes stationarity

# the naive model forecasts quite badly (reverts to the long-term mean)
layout(1:2)
plot_mvgam_fc(m_bad)
plot_mvgam_fc(m_gam_ar)
layout(1)

loo_compare(m_gam_ar, m_bad) #' `m_bad` is clearly much worse
loo_compare(m_gam_ar, m_gam) #' `m_gam_ar` performs better: 12 / 1.6 = 7.5 SEs

# predicting from new data (no terms to predict for the bad model)
plot_predictions(m_gam_ar, condition = "month", points = 0.5)
plot_predictions(m_gam_ar, condition = "year", points = 0.5)

# rates of change: useful for link and expected scales, not for response scale
newd_slopes <- tibble(time = 1:100,
                      year = mean(data_train$year),
                      month = seq(0, 12, length.out = length(time)),
                      lag_12_passengers = mean(data_train$lag_12_passengers))

plot_grid(
  draw(m_gam_ar$trend_mgcv_model, select = "s(month)",
       data = tibble(lag_12_passengers = 0,
                     month = seq(0, 12, length.out = 400),
                     year = 0)),
  plot_slopes(m_gam_ar, variables = "month", by = "month", type = "link",
              newdata = newd_slopes) +
    geom_hline(yintercept = 0, linetype = "dashed"),
  ncol = 1)

plot_slopes(m_gam_ar, variables = "month", by = "month", type = "expected",
            newdata = newd_slopes) +
  geom_hline(yintercept = 0, linetype = "dashed")

# rates of change are too dramatic on the response scale because observations
# are too stochastic
plot_slopes(m_gam_ar, variables = "month", by = "month", type = "response",
            newdata = newd_slopes) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-1000, 1000))

# assessing model fits with forecasts ----
# measures for a good forecast:
# - reliability: is able to predict unobserved data well
# - sharpness: can produce precise predictions with low uncertainty
# - skill: can predict more accurately than other models or a baseline
# 
# we assess the fit of models by evaluating the fits for data that were not
# fit to the model, such as with (leave-one-out) cross-validation. however,
# we cannot simply leave out a subset of our data because data close in time
# are correlated, so data dropped randomly are not independent from those not
# used to fit the model.
# leave-future-out CV allows us to assess the model using future data that are
# less correlated on the observed data.

# point-based forecast evaluation
# 
#' given a forecast horizon, H:
#' - forecast error: `e = obs - pred` at a sufficiently distant future time
#' - estimate the point-based measure:
#'   - mean absolute error: `mean(abs(e))`
#'   - mean squared error: `mean((e)^2)`; similar to variance
#'   - root mean squared error: `sqrt(mean((e)^2))`; similar to SD
#'   - mean abs % error: `100 * mean(abs(e / k)`; scale independent if `k > 0`
#'     `k` can be observations or another benchmark. for more info:
#'     `https://www.youtube.com/watch?v=ek5xLEoQN3E`

# interval-based forecast evaluation
# the scaled interval score (SIS) evaluates forecasts based on deviation from
# an interval of y (e.g., a credible interval)
#' for more info, see: `https://doi.org/10.1371/journal.pcbi.1008618`
calculate_sis <- function(u, l, alpha, y) {
  sis <- case_when(l >= y & y <= u ~ u - l, # in the [l, u] interval
                   y < l ~ u - l + 2 / alpha * (l - y), # below the interval
                   y > u ~ u - l + 2 / alpha * (y - u)) # above the interval
}

# Probabilistic (i.e., distribution-based) forecast evaluation
# 
# rather than focusing on point-specific estimates of model performance, it's
# better to look at the full forecast distribution rather than just a subset of
# points
# this is what we've been doing when comparing forecasts with DRPS/CRPS!
# DRPS/CRPS: Discrete/Continuous Ranked Probability Score
plot(forecast(m_gam_ar))

# for example, we only observe some values from a distribution, but we can
# still estimate the model for the full estimated distribution!
# the probabilities are conditional on the observed data and model parameters
set.seed(20)
d_probs <- tibble(y = seq(-4, 4, by = 0.01),
                  dens = dnorm(y),
                  log_dens = dnorm(y, log = TRUE),
                  sampled = sample(c(TRUE, FALSE), size = length(y), replace = TRUE,
                                   prob = c(0.005, 0.995)))

d_obs <- filter(d_probs, sampled)

# probability densities
ggplot() +
  geom_area(aes(y, dens), d_probs, fill = "grey", color = "black") +
  geom_rug(aes(y), d_obs, lwd = 1, color = "red4") +
  geom_segment(aes(x = y, xend = y, y = 0, yend = dens),
               d_obs, color = "red4") +
  geom_segment(aes(x = y, xend = -Inf, y = dens, yend = dens),
               d_obs, color = "red4",
               arrow = arrow(angle = 15, type = "closed")) +
  ylab("Probability density")

# log probability density
ggplot() +
  geom_area(aes(y, log_dens), d_probs, fill = "grey", color = "black") +
  geom_rug(aes(y), d_obs, lwd = 1, color = "red4") +
  geom_segment(aes(x = y, xend = y, y = 0, yend = log_dens),
               d_obs, color = "red4") +
  geom_segment(aes(x = y, xend = -Inf, y = log_dens, yend = log_dens),
               d_obs, color = "red4",
               arrow = arrow(angle = 15, type = "closed")) +
  ylab("Log(probability density)")

# taking log(density):
# - helps keep numbers more manageable
# - moves operations to the additive scale
# - can be sensitive to outliers (i.e., very negative log(density))

# Continuous Ranked Probability Score
# "-1" forces a flip in the pdf, making the observation the MLE value
if (FALSE) {
  if_else(est < y,
          integrate(f = (pnorm(est)    )^2, lower = -Inf, upper = Inf),
          integrate(f = (pnorm(est) - 1)^2, lower = -Inf, upper = Inf))
}
# SIS converges to CRPS when evaluating many of equally spaced intervals for SIS

# Discrete Ranked Probability Score (CRPS for discrete random variables)
# "-1" forces a flip in the pdf, making the observation the MLE value
if (FALSE) {
  if_else(est < y,
          sum(f = (ppois(est)    )^2, lower = 0, upper = Inf),
          sum(f = (ppois(est) - 1)^2, lower = 0, upper = Inf))
}

#' `{mvgam}` can produce scores quickly and easily
score(forecast(m_gam_ar))$series1 %>% head() # CRPS by default; = DRPS in output
score(forecast(m_gam_ar))$all_series %>% head()

# DRPS produces same values since P(Y = x) = 0 at non-integer values of x
score(forecast(m_gam_ar), score = "drps")$series1 %>% head()
score(forecast(m_gam_ar), score = "drps")$all_series %>% head()

#' Expected log predictive density, as used in `loo_compare()`
score(forecast(m_gam_ar, type = "link"), score = "elpd")$series1 %>% head()

# averaging SIS values across many intervals will converge to DPRS/CPRS values
score(forecast(m_gam_ar), score = "sis")$series1 %>% head()
score(forecast(m_gam_ar), score = "sis", interval_width = 0.5)$series1 %>% head()
score(forecast(m_gam_ar), score = "sis")$all_series %>% head()
