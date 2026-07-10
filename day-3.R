source("packages.R") # attach necessary packages

#' recap:
#' - ARIMA models:
#'   - AR if *data* are correlated through time
#'   - I if data need to be detrended by taking differences
#'   - MA if *errors* are correlated through time
#' - CAR models are a continuous-time version of AR models
#' - GAM's smooth terms are a continuous version of random effects
#' - to estimate change, data should have 3+ observations per period of interest
#' - state space models separate the observation process from the latent process

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

#' today's topics:
#' - forecasting from dynamic models
#' - interpreting the different types of predictions
#' - comparing models and assessing them with forecasts:
#'   - Point-based forecast evaluation
#'   - Probabilistic forecast evaluation
#'   - Bayesian posterior predictive checks

# forecasting from dynamic models ----

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

ggplot(air_passengers, aes(dec_date, passengers, lty = year < 1955)) +
  geom_line() +
  geom_vline(xintercept = 1955, lty = "dashed") +
  labs(x = "Year CE", y = "International airline passengers (thousands)") +
  scale_linetype_manual("Dataset", values = c(3, 1), labels = c("Test", "Train"))

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
               control = list(adapt_delta = 0.9),
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
                  control = list(adapt_delta = 0.95),
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
plot(loo(m_gam_ar), diagnostic = "ESS") #' same as `diagnostic = "n_eff"`
layout(1)

#' Dynamic `{mvgam}` models contain draws for many quantities, all stored as MCMC
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

#' LV and LV_raw are latent variables from AR trend
#' ypred are predictions
#' mus are estimated means
#' trend

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

# generate forecasts for up to end of 2026
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
tail(preds_2026)
for(i in which(is.na(preds_2026$passengers))) {
  if(is.na(preds_2026$lag_12_passengers[i])) {
    preds_2026$lag_12_passengers[i] <- preds_2026$passengers[i - 12]
  }
  
  preds_2026$passengers[i] <-
    predict(m_gam_ar, preds_2026[i, ], type = "expected")[, "Estimate"]
}
tail(preds_2026)

# to stop from calculating score (since values are predicted)
preds_2026 <- rename(preds_2026, estimate = passengers)

# model predicts that passengers will exceed 8 million by the end of 1972
print(preds_2026, n = 13) #' `lag_12_passengers` for row 13 is `estimate` for 1
max(preds_2026$dec_date)
plot_mvgam_fc(m_gam_ar, newdata = preds_2026, ylim = c(0, 8e9 / 1e3),
              realisations = TRUE)
plot_mvgam_fc(m_gam_ar, newdata = preds_2026, ylim = c(0, 8e9 / 1e3))
abline(v = 516, lty = "dashed")
filter(preds_2026, time == 516)

#' `plot(forecast(m_gam_ar, newdata = preds_2026))` fails with error:
#' `arguments imply differing number of rows: 780, 840`
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

# the different scales with count data:
# - link: values are all real numbers (+ or -); scale is additive
# - expected: values are > 0; scale is multiplicative
# - response: values are > 0; scale is multiplicative

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
#' intercept = est. mean response, averaged across all smooths, on expected scale
draw(m_gam_ar$trend_mgcv_model, fun = exp) # relative change in passengers
draw(m_gam_ar$trend_mgcv_model, # partial change in passengers
     fun = \(y) exp(coef(m_gam_ar$trend_mgcv_model)["(Intercept)"]) * exp(y))
plot(forecast(m_gam_ar, type = "expected"))
head(as.vector(forecast(m_gam_ar, type = "expected")$forecasts$series1))

#' *NOTE:* in `{mgcv}` and `{gratia}`, predictions on the "response" scale are
#'         actually on the expected scale, since they only include uncertainty
#'         in the mean
#'         in `{mvgam}`, predictions on the response scale include uncertainty
#'         at the observation level, rather than just the mean
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

pp_check(m_gam_ar, type = "stat", ndraws = 1000, stat = "mean", binwidth = 1)
pp_check(m_gam_ar, type = "stat", ndraws = 1000, stat = "median", binwidth = 2.5)
pp_check(m_gam_ar, type = "stat", ndraws = 1000, stat = "sd", binwidth = 2.5)
pp_check(m_gam_ar, type = "stat", ndraws = 1000, stat = "var", binwidth = 250)

pp_check(m_gam_ar, type = "stat_2d", ndraws = 1000, stat = c("mean", "sd"))

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

#' expected log posterior density: `mean(log(lik(y|model)))`
#' higher ELPD is better, but ELPD can be affected substantially by outliers
#' `m_gam` is best
#' `m_gam_ar` is worse
#' `m_bad` is clearly much worse
#' but ELPD can give too much importance to data in the tails because of it uses
#' the logged likelihood
#' as we'll see later, the score for `m_gam_ar` is penalized too much due to a
#' single uncertain prediction with `k_psis > 0.7`
#' `m_gam_ar` is actually the best model of the three
loo_compare(m_gam_ar, m_bad, m_gam)

#' can calculate ELPD quickly and easily with `{mvgam}`
#' requires `type = "link"`
score(forecast(m_gam_ar, type = "link"), score = "elpd")$series1 %>% head()

tibble(
  model_name = c("m_bad", "m_gam", "m_gam_ar"),
  elpd_data = map(model_name, \(m_n) {
    m <- get(m_n)
    bind_cols(data_test,
              score(forecast(m, type = "link"), score = "elpd")$series1)
  })) %>%
  unnest(elpd_data) %>%
  rename(ELPD = score) %>%
  ggplot(aes(eval_horizon, ELPD, color = model_name)) +
  geom_line(lwd = 1) +
  xlab("Forecast horizon (months)") +
  scale_color_highcontrast() +
  theme(legend.position = "top")

# predicting from new data (can't predict from bad model: no terms to plot!)
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

#' **break**

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
p_forecasts <- plot_grid(plot(forecast(m_bad)) + ggtitle("m_bad"),
                         plot(forecast(m_gam)) + ggtitle("m_gam"),
                         plot(forecast(m_gam_ar)) + ggtitle("m_gam_ar"),
                         ncol = 1)
p_forecasts

# point-based forecast evaluation
#' forecast error: `e = obs - pred` at a sufficiently distant future time
#' *NOTE:* forecasting at the response scale because that is what we observe
forecasts <- tibble(
  model_name = c("m_bad", "m_gam", "m_gam_ar"),
  forecasts = map(model_name, \(m_n) {
    m <- get(m_n)
    mutate(data_test,
           Estimate = predict(m, newdata = data_test, type = "response")[, 1],
           e = passengers - Estimate)
  })) %>%
  unnest(forecasts)

ggplot(forecasts, aes(time, e, color = model_name)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_line(lwd = 1) +
  scale_color_highcontrast() +
  labs(x = "Forecast horizon (months)", y = "Observed - Predicted") +
  theme(legend.position = "top")

# averaged estimates: average across time series at different forecast horizons
sim_ts <-
  tibble(ts = 1:100,
         sim = map(ts, function(i) {
           data_test %>%
             mutate(sim = i,
                    passengers = rpois(n = n(), lambda = passengers))
         })) %>%
  unnest(sim) %>%
  mutate(forecast_horizon = time - max(data_train$time)) %>%
  nest(sims = everything()) %>%
  expand_grid(
    model_name = c("m_bad", "m_gam", "m_gam_ar")) %>%
  mutate(forecasts = map2(model_name, sims, \(m_n, .d) {
    m <- get(m_n)
    bind_cols(.d, predict(m, newdata = .d, type = "response")) %>%
      mutate(e = passengers - Estimate)
  })) %>%
  select(! sims) %>%
  unnest(forecasts)

ggplot(sim_ts) +
  geom_line(aes(time, passengers, group = sim), alpha = 0.1) +
  geom_line(aes(time, passengers), data_test, color = "darkorange", lwd = 0.5)

pointwise_forecast_scores <-
  sim_ts %>%
  group_by(time, model_name) %>%
  summarize(
    #' mean absolute error
    mae = mean(abs(e)),
    #' mean squared error (similar to variance)
    mse = mean(e^2),
    #' - root mean squared error; similar to SD
    rmse = sqrt(mse),
    #' mean abs % error: `100 * mean(abs(e / k))`, where  `k` can be the
    #' observation or another benchmark
    #' for more info: https://www.youtube.com/watch?v=ek5xLEoQN3E
    mape = 100 * mean(abs(e / passengers)),
    .groups = "drop") %>%
  pivot_longer(mae:mape, names_to = "metric") %>%
  mutate(metric = toupper(metric))

ggplot(pointwise_forecast_scores, aes(time, value, color = model_name)) +
  facet_grid(metric ~ ., scales = "free_y", switch = "y") +
  geom_line() +
  labs(x = "Forecast horizon (months)", y = NULL) +
  scale_color_highcontrast(name = "Model") +
  theme(legend.position = "top", strip.placement = "outside",
        strip.background = element_blank())

p_forecasts # compare to forecast plots

# interval-based forecast evaluation
# the scaled interval score (SIS; https://doi.org/10.1198/016214506000001437 )
# evaluates forecasts based on deviation from a (credible) interval of y
#' `l`, `u` = lower and upper bounds of an interval
#' `alpha` = 1 - interval coverage
#' `y` = observation without measurement error
calculate_sis <- function(l, u, alpha, y) {
  case_when(y >= l & y <= u ~ u - l,
            y < l           ~ u - l + (l - y) / (alpha / 2),
            y > u           ~ u - l + (y - u) / (alpha / 2))
}

sim_ts <- sim_ts %>%
  mutate(sis = calculate_sis(Q2.5, Q97.5, alpha = 0.05, y = passengers))

interval_forecast_scores <-
  sim_ts %>%
  group_by(time, model_name) %>%
  summarize(mean_sis = mean(sis), .groups = "drop")

ggplot(interval_forecast_scores, aes(time, mean_sis, color = model_name)) +
  geom_line() +
  labs(x = "Forecast horizon (months)", y = "Mean SIS") +
  scale_color_highcontrast(name = "Model") +
  theme(legend.position = "top")

p_forecasts # compare to forecast plots

#' can calculate SIS quickly and easily with `{mvgam}`
score(forecast(m_gam_ar), score = "sis")$series1 %>% head()

# Probabilistic (i.e., distribution-based) forecast evaluation
# 
# rather than focusing on specific points or intervals, it's better to look at
# the full forecast distribution of the forecast and account for the performance
# of all parameters at once (even latent ones).
# this is what we've been doing when comparing forecasts with DRPS/CRPS, the
# Discrete/Continuous Ranked Probability Score
#' for more info, see: https://doi.org/10.1371/journal.pcbi.1008618
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

plot_grid(
  # probability densities
  ggplot() +
    geom_area(aes(y, dens), d_probs, fill = "grey", color = "black") +
    geom_rug(aes(y), d_obs, lwd = 1, color = "red4") +
    geom_segment(aes(x = y, xend = y, y = 0, yend = dens),
                 d_obs, color = "red4") +
    geom_segment(aes(x = y, xend = -Inf, y = dens, yend = dens),
                 d_obs, color = "red4",
                 arrow = arrow(angle = 15, type = "closed")) +
    ylab("Probability density"),
  
  # log probability density
  ggplot() +
    geom_area(aes(y, log_dens), d_probs, fill = "grey", color = "black") +
    geom_rug(aes(y), d_obs, lwd = 1, color = "red4", sides = "t") +
    geom_segment(aes(x = y, xend = y, y = 0, yend = log_dens),
                 d_obs, color = "red4") +
    geom_segment(aes(x = y, xend = -Inf, y = log_dens, yend = log_dens),
                 d_obs, color = "red4",
                 arrow = arrow(angle = 15, type = "closed")) +
    ylab("Log(probability density)"),
  ncol = 1)

# taking log(density):
# - helps keep numbers more manageable
# - moves operations to the additive scale
# - is sensitive to outliers (i.e., very negative log(density), like with ELPD)
# - requires distributional/parametric assumptions

# Continuous Ranked Probability Score
# "-1" forces a flip in the pdf, making the observation the MLE value
if (FALSE) { # do not run: for illustration purposes only
  if_else(est < y,
          integrate(f = (pnorm(est)    )^2, lower = -Inf, upper = Inf),
          integrate(f = (pnorm(est) - 1)^2, lower = -Inf, upper = Inf))
}
# SIS converges to CRPS when evaluating many of equally spaced intervals for SIS

# Discrete Ranked Probability Score (CRPS for discrete random variables)
# "-1" forces a flip in the pdf, making the observation the MLE value
if (FALSE) { # do not run: for illustration purposes only
  if_else(est < y,
          sum(f = (ppois(est)    )^2, lower = 0, upper = Inf),
          sum(f = (ppois(est) - 1)^2, lower = 0, upper = Inf))
}

#' `{mvgam}` can produce scores quickly and easily
score(forecast(m_gam_ar))$series1 %>% head() # CRPS by default; = DRPS in output
score(forecast(m_gam_ar))$all_series %>% head()

# DRPS produces same values since P(Y = x) = 0 at non-integer values of x
score(forecast(m_gam_ar), score = "drps")$series1 %>% head()
score(forecast(m_gam_ar), score = "crps")$series1 %>% head()

dprs_scores <- tibble(
  model_name = c("m_bad", "m_gam", "m_gam_ar"),
  dprs_tib = map(model_name, \(m_n) {
    mutate(data_test,
           DRPS = score(forecast(get(m_n)), score = "drps")$series1$score)
  })) %>%
  unnest(dprs_tib)

ggplot(dprs_scores, aes(time, DRPS, color = model_name)) +
  geom_line() +
  labs(x = "Forecast horizon (months)") +
  scale_color_highcontrast(name = "Model") +
  theme(legend.position = "top")
