source("packages.R") # attach necessary packages

#' today's topics:
#' - multivariate ecological time series
#' - vector autoregressive processes
#' - dynamic factor models
#' - multivariate forecast evaluation

# multiple cores from the same lake: multiple time series for the same pigments
#' core locations: `https://onlinelibrary.wiley.com/cms/asset/a02d5fe1-044e-4dd2-b50c-ff33f328d953/fwb14192-fig-0001-m.jpg`
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

# Multivariate ecological time series ----
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
m_null <- mvgam(formula = diatox ~ core,
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
#' `https://doi.org/10.7717/peerj.6876/fig-4`
m_gam <- mvgam(formula = diatox ~ core + s(year, core, bs = "fs", k = 10),
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
plot_grid(draw(m_gam$mgcv_model), # notice the small wiggles
          draw(m_gam_ar1$mgcv_model), #' terms are smoother than for `m_gam`
          ncol = 1)

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

#' **break**

# Vector autoregressive processes ----
#' y values are correlated in time; use correlated `AR(1)` processes: `VAR(1)`
#' `VAR(1)` processes allow us to:
#' - predict missing observations d using AR processes from other series
#' - predict "pulses" that may not be estimable efficiently with scarce data
#' - model proxies efficiently if they are more convenient to monitor
#' 
#' we fit `VAR(1)` processes by setting `y_t ∼ Normal(A ∗ y_{t−1}, Σ)`, where:
#' - we have `s` series,
#' - `y_t` is the vector with length `s` of observations at time `t`,
#' - `y_{t-1}` is the vector of `s` observations at time `t-1`,
#' - `A` is the `s * s` matrix of correlations for `y_t` and `y_{t-1}` values,
#' - `Σ` is the covariance matrix that determines the correlation across errors
#' can safely ignore messages about rejections of initial values
#' for more info: `nicholasjclark.github.io/mvgam/articles/trend_formulas.html`

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

# fits in ~ 2-3 minutes
#' drop intercepts to avoid non-identifiability issues (`noncentered = FALSE`)
#' trying to fit the model with a Gamma distribution results in horrible chains
m_gam_var1 <- mvgam(formula = z_log ~ 0 + s(year, core, bs = "fs", k = 10),
                    trend_model = VAR(1, cor = TRUE), # vector AR(1) process
                    noncentred = FALSE, # cannot use if using for VAR() models
                    family = gaussian(), # since we centered the data
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

plot(m_gam_var1) # no appreciable issues with ACF or pACF
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
# also, estimated means are the geometric mean because of the log-transform
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

# Dynamic factor models ----
#' allow observations to depend on latent variables: `y_t = θ * z_t`, where:
#' - `y_t` is a vector of observations at time `t`,
#' - `z_t` is a vector of dynamic factor estimates at time `t`,
#' - `θ` is a matrix of loading coefficients that controls how each series in
#'   `y` depends on the factors in `z`
#' for more info, see:
#' - `https://www.youtube.com/watch?v=FMLh8c_Sa8s`
#' - `cran.r-project.org/web/packages/dfms/vignettes/dynamic_factor_models.pdf`
#' note: can ignore warning about integer division. it's caused by forced
#' rounding from a decimal to an integer in the stan code
m_df <- mvgam(formula = diatox ~ s(core, bs ="re"), # assume stationarity
              use_lv = TRUE, n_lv = 4, # needs to be chosen manually
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
                  use_lv = TRUE, n_lv = 4,
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
draw(m_gam_df$mgcv_model) # terms are again quite smooth

# trend terms now change less over time because smooth terms account for change
plot_grid(plot(m_gam_df, type = "trend", series = 1),
          plot(m_gam_df, type = "trend", series = 2),
          plot(m_gam_df, type = "trend", series = 3),
          plot(m_gam_df, type = "trend", series = 4))

# predictions aren't too bad
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

# Multivariate forecast evaluation ----

# the Energy Score generalizes the continuous ranked probability score to
# multivariate forecasts
#' `ES(F, y) = mean(norm(F_i − y)) − mean(norm(F_i − F_j))`,
#' where:
#' - `F` is a vector of `m` forecasts with elements indicated by `F_i` and `F_j`,
#' - `y` is a vector of future observations for the corresponding `m` times,
#' - `norm(k)` indicates the length or Euclidean norm of a vector, and
#' - `mean(F_i − F_j) = 1 / (2 * m^2) * sum(F_i − F_j)`.
#' A lower energy score is better: it indicates less deviation from the data.
#' 
#' the Energy Score is a weighted distance between the distribution of the
#' forecasts and the distribution of the observations
#' 
#' useful for generalizing CPRS but maintaining a simple structure while also
#' accounting for sharpness and calibration
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

# AR, VAR, and DFA models only account for stochastic component
# we can extend the penalty for the GAM using b-splines
# but we need separate smooths for each core because fs terms can't use bs bases
# b-splines can be useful, but extrapolating with GAMs is hard if predictors
# aren't truly deterministic (i.e., not a smooth of time)
#' for more info: `fromthebottomoftheheap.net/2020/06/03/extrapolating-with-gams/`
#' using `by` smooths because `fs` smooths don't support `bs` basis
m_gam_ar1_bs <- mvgam(formula = diatox ~
                        core + #' `by` smooths require explicit intercepts
                        s(year, by = core, bs = "bs", m = c(3, 2), k = 10),
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

# predictions are worse than other models because they depend on assumptions of
# how the penalty carries forward. we penalized deviations from the curvature
plot_grid(plot(forecast(m_gam_ar1_bs), series = 1),
          plot(forecast(m_gam_ar1_bs), series = 2),
          plot(forecast(m_gam_ar1_bs), series = 3),
          plot(forecast(m_gam_ar1_bs), series = 4))

# compare models based on their forecast scores
models <-
  tibble(name = factor(c("null", "GAM only", "GAM with AR(1)", "GAM with VAR(1)",
                         "AR(1) and DF(2)", "GAM with AR(1) and DF(2)",
                         "GAM with AR(1) and b-splines"),
                       levels = c("null", "AR(1) and DF(4)", "GAM only",
                                  "GAM with AR(1)", "GAM with VAR(1)",
                                  "GAM with AR(1) and DF(4)",
                                  "GAM with AR(1) and b-splines")),
         model = list(m_null, m_gam, m_gam_ar1, m_gam_var1,
                      m_df, m_gam_df, m_gam_ar1_bs),
         forecast = map(model, function(.m) forecast(.m)),
         scores = map(forecast, function(.f) {
           tibble(
             year = unique(d_test$year),
             time = unique(d_test$time),
             energy_score = score(.f, score = "energy")$all_series$score,
             variogram_score = score(.f, score = "variogram")$all_series$score)
         }))

# lower energy and variogram scores are better: less deviation from data
# cannot compare the GAM with VAR(1) because it used transformed data
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
