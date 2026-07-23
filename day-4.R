source("packages.R") # attach necessary packages

#' recap:
#' - ARMA and CAR models help deal with autocorrelation in data and residuals
#' - GAM's smooth terms are a continuous version of random effects
#' - to estimate change, data should have 3+ observations per period of interest
#' - state space models separate the observation process from the latent process
#' - forecasting models are best assessed using out-of-sample prediction, and
#'   with distributional scores such as DRPS and CRPS

#' today's topics:
#' - evaluating model forecasts
#' - multivariate ecological time series
#' - vector autoregressive processes
#' - dynamic factor models
#' - multivariate forecast evaluation

# assessing model fits with forecasts ----
air_passengers <-
  tibble(time = 1:length(AirPassengers),
         dec_date = time(AirPassengers) %>% as.numeric(),
         year = floor(dec_date),
         month = round((dec_date - year) * 12) + 1, # to ensure January is 1
         passengers = as.numeric(AirPassengers)) %>% # in thousands
  mutate(lag_12_passengers = lag(passengers, 12))
data_train <- filter(air_passengers, year <= 1955) %>%
  filter(! is.na(lag_12_passengers)) # lagged value before 1st observation is NA
data_test <- filter(air_passengers, year > 1955)

m_bad <- mvgam(formula = passengers ~ 0,
               trend_formula = ~ 1,
               trend_model = AR(p = 1), # AR(1) model
               noncentred = TRUE, # use a noncentered AR(1) model
               family = poisson(link = "log"),
               data = data_train, # calculate forecast while fitting
               newdata = data_test,
               chains = 4,
               burnin = 750,
               samples = 500,
               parallel = TRUE,
               silent = 2)
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
  geom_line(aes(time, passengers), data_test, color = "darkorange", lwd = 1)

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
xf
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

#' *break*

# modelling multiple time series ----
# multiple cores from the same lake: multiple time series for the same pigments
#' core locations: https://onlinelibrary.wiley.com/cms/asset/a02d5fe1-044e-4dd2-b50c-ff33f328d953/fwb14192-fig-0001-m.jpg
SAMPLING_DATE <- lubridate::decimal_date(as.POSIXlt("2014-04-01"))

pigments <-
  bind_rows(
    read.xlsx("https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%201%20April%202014.xlsx") %>%
      mutate(core = "Core 1"),
    read.xlsx("https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%202%20April%202014.xlsx") %>%
      mutate(core = "Core 2") %>%
      rename(CHLA = CHL_A),
    read.xlsx("https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%203%20April%202014.xlsx") %>%
      mutate(core = "Core 3"),
    read.xlsx("https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%204%20April%202014.xlsx") %>%
      mutate(core = "Core 4")) %>%
  select(core, YEAR, DIATOX) %>%
  rename_with(tolower) %>%
  filter(! is.na(year)) %>%
  mutate(year = round(year),
         interval = c(SAMPLING_DATE, year[-length(year)]) - year,
         weight = interval / mean(interval),
         .by = core) %>%
  summarize(diatox = mean(diatox),
            .by = c(core, year)) %>%
  as_tibble() %>%
  right_join(expand_grid(year = seq(min(.$year), max(.$year), by = 1),
                         core = unique(.$core)),
             by = c("year", "core")) %>%
  mutate(core = factor(core),
         series = core,
         time = year) %>%
  arrange(time, core) %>%
  mutate(diatox = if_else(diatox > 200 & core == "Core", NA_real_, diatox))

ggplot(pigments, aes(year, diatox)) +
  facet_wrap(. ~ core) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "gam", color = "black", formula = y ~ s(x, k = 15)) +
  ylim(c(0, NA))

# multivariate ecological time series ----
# we want to allow the model to learn about trends ...
ggplot(pigments, aes(year, diatox, color = core, fill = core)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 15)) +
  scale_fill_bright(name = "Core") +
  scale_color_bright(name = "Core") +
  ylim(c(0, NA))

# ... but limit the effects of outlier/odd data and periods ...
ggplot(pigments, aes(year, diatox, color = core)) +
  geom_point(alpha = 0.4) +
  scale_color_bright(name = "Core") +
  ylim(c(0, NA))

# ... and only fit the model to data prior to 1950
d_train <- filter(pigments, year < 1950)
d_test <- filter(pigments, year >= 1950) %>%
  # round to nearest 5 for easier comparisons
  # scoring across series need all values to be non-NA
  mutate(year = year - year %% 5,
         time = year) %>%
  summarize(diatox = mean(diatox, na.rm = TRUE),
            .by = c(core, year, series, time)) %>%
  right_join(expand_grid(core = unique(.$core),
                         year = 1950:2010) %>%
               mutate(series = core,
                      time = year),
             by = c("year", "time", "core", "series")) %>%
  arrange(time, series)

ggplot(mapping = aes(year, diatox, color = core)) +
  geom_point(data = d_train, alpha = 0.4) +
  geom_vline(xintercept = 1950, lty = "dashed") +
  geom_point(data = d_test, alpha = 0.4, pch = 4) +
  scale_color_bright(name = "Core") +
  ylim(c(0, NA))

# fit a null model
m_null <- mvgam(formula = diatox ~ s(core, bs = "re"),
                family = Gamma(link = "log"),
                data = d_train,
                newdata = d_test,
                chains = 4,
                burnin = 500,
                sample = 750,
                parallel = TRUE,
                silent = 2)
plot(loo(m_null))

# unsurprisingly, the fit is horrible
plot_grid(plot(forecast(m_null), series = 1),
          plot(forecast(m_null), series = 2),
          plot(forecast(m_null), series = 3),
          plot(forecast(m_null), series = 4))

# fit a hierarchical GAM
#' https://doi.org/10.7717/peerj.6876/fig-4
m_gam <- mvgam(formula = diatox ~ s(year, core, bs = "fs", k = 10),
               family = Gamma(link = "log"),
               data = d_train,
               newdata = d_test,
               chains = 4,
               burnin = 500,
               samples = 750,
               parallel = TRUE)
summary(m_gam) # no issues with diagnostics
plot(m_gam) # clear autocorrelation at lag 1
draw(m_gam$mgcv_model) # notice the small wiggles in the fs smooths
plot(loo(m_gam))

# accuracy of forecasts isn't bad, but it's also quite uncertain
# the forecasts assume the trends continue over time
#' NOTE: `plot(m_gam, type = "forecast")` uses base `graphics::plot()` 
plot_grid(plot(forecast(m_gam), series = 1),
          plot(forecast(m_gam), series = 2),
          plot(forecast(m_gam), series = 3),
          plot(forecast(m_gam), series = 4))

# add AR(1) process
m_gam_ar1 <- mvgam(formula = diatox ~ s(year, core, bs = "fs", k = 10),
                   trend_model = AR(1), # autoregressive 1 process
                   noncentred = TRUE,
                   family = Gamma(link = "log"),
                   data = d_train,
                   newdata = d_test,
                   chains = 4,
                   burnin = 500,
                   sample = 750,
                   parallel = TRUE)

summary(m_gam_ar1) # no issues with diagnostics
plot(m_gam_ar1) # no appreciable autocorrelation at lag 1

#' credible intervals are wider than with `m_gam`, esp. between data points
#' but trends in `m_gam_ar1` are more jagged
layout(matrix(1:8, ncol = 2))
for(i in 1:4) plot(m_gam, type = "forecast", series = i)
for(i in 1:4) plot(m_gam_ar1, type = "forecast", series = i)
layout(1)

# trend term peaks when data deviate most from the estimated mean
plot_grid(plot(m_gam_ar1, type = "trend", series = 1),
          plot(m_gam_ar1, type = "trend", series = 2),
          plot(m_gam_ar1, type = "trend", series = 3),
          plot(m_gam_ar1, type = "trend", series = 4))

# again, forecasts aren't great
plot_grid(plot(forecast(m_gam_ar1), series = 1),
          plot(forecast(m_gam_ar1), series = 2),
          plot(forecast(m_gam_ar1), series = 3),
          plot(forecast(m_gam_ar1), series = 4))

# vector autoregressive processes ----
#' y values are correlated in time; use correlated `AR(1)` processes: `VAR(1)`
#' `VAR(1)` processes allow us to:
#' - predict missing observations using AR processes from other series
#' - predict "pulses" that may not be estimable efficiently with scarce data
#' - model proxies efficiently if they are more convenient to monitor
#' 
#' we fit `VAR(1)` processes by setting `y_t ∼ Normal(A ∗ y_{t−1}, Σ)`, where:
#' - we have `s` series,
#' - `y_t` is the vector with length `s` of observations at time `t`,
#' - `y_{t-1}` is the vector of `s` observations at time `t-1`,
#' - `A` is the `s * s` matrix of autocorrelations for `y_t` and `y_{t-1}` values,
#' - `Σ` is the covariance matrix that determines the correlation across values
#' can safely ignore messages about rejections of initial values
#' for more info: nicholasjclark.github.io/mvgam/articles/trend_formulas.html

d_train <- mutate(d_train,
                  log_diatox = log(diatox),
                  sample_mean = mean(log_diatox, na.rm = TRUE),
                  sample_sd = sd(log_diatox, na.rm = TRUE),
                  z_log = (log_diatox - sample_mean) / sample_sd,
                  .by = core)
d_test <- left_join(d_test,
                    d_train %>% slice(1, .by = core) %>%
                      select(core, sample_mean, sample_sd),
                    by = "core") %>%
  mutate(log_diatox = log(diatox),
         z_log = (log_diatox - sample_mean) / sample_sd)

# plot centered data
ggplot(d_train, aes(year, z_log)) +
  facet_wrap(~ core) +
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_vline(xintercept = min(d_test$year), lty = "dashed") +
  geom_point() +
  geom_point(data = d_test, pch = 4)

#' *fits in ~ 5 minutes*
#' drop intercepts to avoid non-identifiability issues (`noncentered = FALSE`)
#' trying to fit the model with a Gamma distribution results in horrible chains
m_gam_var1 <- mvgam(formula = z_log ~ 0 + s(year, core, bs = "fs", k = 10),
                    trend_model = VAR(1, cor = TRUE), # vector AR(1) process
                    noncentred = FALSE, # cannot use if using for VAR() models
                    family = student_t(), # since we centered the data
                    data = d_train,
                    newdata = d_test,
                    chains = 4,
                    burnin = 200,
                    sample = 750,
                    control = list(adapt_delta = 0.95),
                    parallel = TRUE)

mcmc_plot(m_gam_var1, type = "trace", variable = ".", regex = TRUE)
summary(m_gam_var1)

matrix(summary(m_gam_var1)$parameters$var_coefficient_matrix[, "50%"],
       ncol = 4, byrow = TRUE)

plot(m_gam_var1) # no appreciable issues with ACF or pACF, but overdispersed
draw(m_gam_var1$mgcv_model)

plot_grid(plot(m_gam_var1, type = "trend", series = 1),
          plot(m_gam_var1, type = "trend", series = 2),
          plot(m_gam_var1, type = "trend", series = 3),
          plot(m_gam_var1, type = "trend", series = 4))

# forecasts are on the centered scale
plot_grid(plot(forecast(m_gam_var1), series = 1),
          plot(forecast(m_gam_var1), series = 2),
          plot(forecast(m_gam_var1), series = 3),
          plot(forecast(m_gam_var1), series = 4))

# predictions require back-transformation, and the CI for core 1 is unreliable
phi_tbl <-
  tibble(core = paste("Core", 1:4),
         phi = summary(m_gam_var1)$parameters$standard_deviation[, "50%"])

pigments %>%
  left_join(
    d_train %>%
      slice(1, .by = core) %>%
      select(core, sample_mean, sample_sd) %>%
      left_join(phi_tbl, by = "core")) %>%
  bind_cols(predict(m_gam_var1, newdata = ., type = "expected")) %>%
  mutate(Estimate = exp(Estimate * sample_sd + sample_mean),
         Q2.5 = exp(Q2.5 * sample_sd + sample_mean),
         Q97.5 = exp(Q97.5 * sample_sd + sample_mean)) %>%
  ggplot() +
  facet_wrap(~ core) +
  coord_cartesian(ylim = c(0, 400)) +
  geom_ribbon(aes(year, ymin = Q2.5, ymax = Q97.5), alpha = 0.3) +
  geom_point(aes(year, diatox), alpha = 0.5) +
  geom_line(aes(year, Estimate)) +
  geom_vline(xintercept = 1950, lty= "dashed")

# dynamic factor models ----
#' - a dynamic factor (DF) are latent time series
#' - DF loadings are coefficients that relate DFs to observed time series
#' - DFMs are similar to PCA, but they allow a vcov matrix with cov != 0
#' - this allows for correlations across time series
#' - allow observations to depend on latent variables: `y_t = θ * z_t`, where:
#'   - `y_t` is a vector of observations at time `t`,
#'   - `z_t` is a vector of dynamic factor estimates at time `t`,
#'   - `θ` is a matrix of loading coefficients that controls how each series in
#'     `y` depends on the `z` dynamic factors (i.e., dimension-reduced ts)
#' - factor: 
#' for more info, see:
#' - https://atsa-es.github.io/atsa2017/Labs/Week%206%20dynamic%20factor%20analysis/Intro_to_DFA.html
#' - https://www.youtube.com/watch?v=FMLh8c_Sa8s
#' - cran.r-project.org/web/packages/dfms/vignettes/dynamic_factor_models.pdf
#' note: can ignore warning about integer division. it's caused by forced
#' rounding from a decimal to an integer in the stan code
m_df <- mvgam(formula = diatox ~ s(core, bs ="re"), # assume stationarity
              use_lv = TRUE, n_lv = 2, # needs to be chosen manually
              trend_model = AR(1),
              noncentred = FALSE, # since trend formula is empty
              family = Gamma(link = "log"),
              data = d_train,
              newdata = d_test,
              chains = 4,
              burnin = 500,
              sample = 750,
              control = list(adapt_delta = 0.95),
              parallel = TRUE)

summary(m_df) #' convergence issues REs and `ar1` coefs
plot(m_df) # strong but uncertain autocorrelation at lag 1
m_df$mgcv_model # there is no mgcv model now (other than intercepts)

# trend term includes a large portion of the change over time
# forecasts are quite poor because they assume stationarity & mean reversion...
plot_grid(plot(m_df, type = "trend", series = 1),
          plot(m_df, type = "trend", series = 2),
          plot(m_df, type = "trend", series = 3),
          plot(m_df, type = "trend", series = 4))

# ... but the trends do follow the data closely, for being such a simple model
plot_grid(plot(forecast(m_df), series = 1),
          plot(forecast(m_df), series = 2),
          plot(forecast(m_df), series = 3),
          plot(forecast(m_df), series = 4))

# add a GAM term (fits faster because it converges better)
m_gam_df <- mvgam(formula = diatox ~ s(year, core, bs = "fs", k = 10),
                  use_lv = TRUE, n_lv = 2,
                  trend_model = AR(1),
                  noncentred = FALSE, #' since `trend_formula` is not empty
                  family = Gamma(link = "log"),
                  data = d_train,
                  newdata = d_test,
                  chains = 4,
                  burnin = 500,
                  sample = 750,
                  control = list(adapt_delta = 0.9),
                  parallel = TRUE)

summary(m_gam_df) # no issues with diagnostics
plot(m_gam_df) # no appreciable autocorrelations
draw(m_gam_df$mgcv_model) # terms are quite smooth
draw(m_gam$mgcv_model)

# trend terms now change less over time because smooth terms account for change
plot_grid(plot(m_gam_df, type = "trend", series = 1),
          plot(m_gam_df, type = "trend", series = 2),
          plot(m_gam_df, type = "trend", series = 3),
          plot(m_gam_df, type = "trend", series = 4))

plot_grid(plot(m_df, type = "trend", series = 1),
          plot(m_df, type = "trend", series = 2),
          plot(m_df, type = "trend", series = 3),
          plot(m_df, type = "trend", series = 4))

# forecasts aren't too bad
plot_grid(plot(forecast(m_gam_df), series = 1),
          plot(forecast(m_gam_df), series = 2),
          plot(forecast(m_gam_df), series = 3),
          plot(forecast(m_gam_df), series = 4))

# the other element of the list is a list of posterior correlations
lv_correlations(object = m_df)$mean_correlations # strong correlations if no GAM
lv_correlations(object = m_gam_df)$mean_correlations # moderate for 2, 3, 4

#' some scores may be unreliable because of high pareto shape values
loo_compare(m_null, m_gam, m_gam_ar1, m_gam_var1, m_df, m_gam_df)

# advantages of VAR processes:
# - can use for assessing whether one series can predict another (often called
#   "Granger causality", but it's not causality -- "precedence" is more correct)
# - can allow us to develop complex correlations across many different variables
#   across time, but can be very unstable and data-hungry
# advantages of DF models:
# - can explain complex models with very few terms
# - can be hard to interpret

#' **break**

# multivariate forecast evaluation ----

# the Energy Score generalizes the continuous ranked probability score to
# multivariate forecasts
#' `ES(F, y) = (mean(Enorm(F_i − y))) − (mean(Enorm(F_i − F_j)))`,
#' where:
#' - `F` is a vector of `m` forecasts with elements indicated by `F_i` and `F_j`,
#' - `y` is a vector of future observations for the corresponding `m` timeseries,
#' - `Enorm(k)` indicates the length or Euclidean norm of a vector, and
#' - `mean(F_i − F_j) = 1 / (2 * m^2) * sum(F_i − F_j)`.
#' A lower energy score is better: it indicates less deviation from the data.
#' 
#' the Energy Score is a weighted distance between the distribution of the
#' forecasts and the distribution of the observations
#' 
#' useful for generalizing CPRS but maintaining a simple structure while also
#' accounting for:
#' - sharpness (i.e., concentration of the forecasts, inverse of uncertainty), &
#' - calibration (i.e., consistency between forecasts and observed data)
#' https://sites.stat.washington.edu/raftery/Research/PDF/Gneiting2007jrssb.pdf
#' 
#' however, does not account for correlations in the test data that were absent
#' in the training data

#' the variogram score quantifies the correlation structure in observations
#' using a measure of the squared distances between pairs of observations:
#' `VG(F',y) = mean(w_{ij} * (sqrt(abs(y_i − y_j)) − sqrt(abs(F'_i − F'_j)))^2`,
#' where:
#' - `F'` is a vector of forecast summaries (e.g., mean or median),
#' - `y` is a vector of corresponding future observations, and
#' - `w_{ij}` are non-negative wights for observations `i` and `j` (often 1).
#' A lower variogram score is better: it indicates less deviation from the data.

#' useful for detecting whether the predictions match the correlation structures
#' in the future data, which may be unobserved in the training data
#' but it relies on a summary (`F'`) of the forecasts for a given time, so it
#' does not account for forecast sharpness or calibration

# extending penalties for better forecasts ---- 
# AR, VAR, and DFA models only account for stochastic component
# we can extend the penalty for the GAM using b-splines
# but we need separate smooths for each core because fs terms can't use bs bases
# b-splines can be useful, but extrapolating with GAMs is hard if predictors
# aren't truly deterministic (i.e., not a smooth of time)
#' for more info: fromthebottomoftheheap.net/2020/06/03/extrapolating-with-gams/
#' using `by` smooths because `fs` smooths don't support `bs` basis
m_gam_ar1_bs <- mvgam(formula = diatox ~
                        core + #' `by` smooths require explicit intercepts
                        s(year, by = core, bs = "bs", m = c(3, 1), k = 10),
                      knots = list(year = c(1800, 2010)),
                      trend_model = AR(1), # autoregressive 1 process
                      noncentred = TRUE,
                      family = Gamma(link = "log"),
                      data = d_train,
                      newdata = d_test,
                      chains = 4,
                      burnin = 500,
                      sample = 750,
                      parallel = TRUE)

summary(m_gam_ar1_bs) # no issues with diagnostics
plot(m_gam_ar1_bs) # no appreciable autocorrelations
draw(m_gam_ar1_bs$mgcv_model) # each term has its own smoothness parameter

# stable trend term after 1950 (1950 - 1801 = 149 on x axis)
plot_grid(plot(m_gam_ar1_bs, type = "trend", series = 1),
          plot(m_gam_ar1_bs, type = "trend", series = 2),
          plot(m_gam_ar1_bs, type = "trend", series = 3),
          plot(m_gam_ar1_bs, type = "trend", series = 4))

# predictions are better than other models when we have correct assumptions of
# how the penalty carries forward. we penalized deviations from the slope
plot_grid(plot(forecast(m_gam_ar1_bs), series = 1),
          plot(forecast(m_gam_ar1_bs), series = 2),
          plot(forecast(m_gam_ar1_bs), series = 3),
          plot(forecast(m_gam_ar1_bs), series = 4))

# compare models based on their forecast scores
models <-
  tibble(name = factor(c("null", "GAM only", "GAM with AR(1)", "GAM with VAR(1)",
                         "AR(1) and DF(2)", "GAM with AR(1) and DF(2)",
                         "GAM with AR(1) and b-splines"),
                       levels = c("null", "AR(1) and DF(2)", "GAM only",
                                  "GAM with AR(1)", "GAM with VAR(1)",
                                  "GAM with AR(1) and DF(2)",
                                  "GAM with AR(1) and b-splines")),
         model = list(m_null, m_gam, m_gam_ar1, m_gam_var1,
                      m_df, m_gam_df, m_gam_ar1_bs),
         forecast = map2(model, name, function(.m, .n) {
           out <- forecast(.m)
           if (.n == "GAM with VAR(1)") {
             for (.core in paste("Core", 1:4)) {
               sm <- filter(d_train, core == .core)$sample_mean[1]
               ssd <- filter(d_train, core == .core)$sample_sd[1]
               out$forecasts[[.core]] <- exp(out$forecasts[[.core]] * ssd + sm)
             }
           }
           }),
         scores = map(forecast, function(.f) {
           tibble(
             year = unique(d_test$year),
             time = unique(d_test$time),
             energy_score = score(.f, score = "energy")$all_series$score,
             variogram_score = score(.f, score = "variogram")$all_series$score)
         }))

#' lower energy and variogram scores are better: less deviation from data
#' *NOTE*: cannot compare the GAM with VAR(1) because it used transformed data
models %>%
  select(name, scores) %>%
  unnest(scores) %>%
  filter(! is.na(energy_score)) %>%
  pivot_longer(c(energy_score, variogram_score), names_to = "type",
               values_to = "score") %>%
  mutate(type = gsub("_score", "", type)) %>%
  ggplot(aes(year, score, color = name, fill = name)) +
  facet_grid(type ~ ., scales = "free_y") +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x), se = FALSE) +
  scale_color_bright(name = "Model") +
  scale_fill_bright(name = "Model")

# extra material if time permits:
# - https://ctmm-initiative.github.io/ctmm/articles/variogram.html
# - https://www.ab.mpg.de/518835/ctmm-1_IntroductionToContinousTimeMovementModelling.pdf
# - https://link.springer.com/article/10.1186/s40462-019-0177-1
