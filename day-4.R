source('packages.R') # attach necessary packages

# multiple cores from the same lake: multiple time series for the same pigments
#' core locations: `https://onlinelibrary.wiley.com/cms/asset/a02d5fe1-044e-4dd2-b50c-ff33f328d953/fwb14192-fig-0001-m.jpg`
SAMPLING_DATE <- lubridate::decimal_date(as.POSIXlt('2014-04-01'))

pigments <-
  bind_rows(
    read.xlsx('https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%201%20April%202014.xlsx') %>%
      mutate(core = 'Core 1'),
    read.xlsx('https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%202%20April%202014.xlsx') %>%
      mutate(core = 'Core 2') %>%
      rename(CHLA = CHL_A),
    read.xlsx('https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%203%20April%202014.xlsx') %>%
      mutate(core = 'Core 3'),
    read.xlsx('https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%204%20April%202014.xlsx') %>%
      mutate(core = 'Core 4')) %>%
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
             by = c('year', 'core')) %>%
  mutate(core = factor(core),
         series = core,
         time = year) %>%
  arrange(time, core)

View(pigments)

ggplot(pigments, aes(year, diatox)) +
  facet_wrap(. ~ core) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'gam', color = 'black', formula = y ~ s(x, k = 15)) +
  ylim(c(0, NA))

# Multivariate ecological time series ----
# we want to allow the model to learn about trends ...
ggplot(pigments, aes(year, diatox, color = core, fill = core)) +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 15)) +
  scale_fill_bright(name = 'Core') +
  scale_color_bright(name = 'Core') +
  ylim(c(0, NA))

# ... but limit the effects of outlier/odd data and periods ...
ggplot(pigments, aes(year, diatox, color = core)) +
  geom_point(alpha = 0.4) +
  scale_color_bright(name = 'Core') +
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
             by = c('year', 'time', 'core', 'series')) %>%
  arrange(time, series)

ggplot(mapping = aes(year, diatox, color = core)) +
  geom_point(data = d_train, alpha = 0.4) +
  geom_vline(xintercept = 1950, lty = 'dashed') +
  geom_point(data = d_test, alpha = 0.4, pch = 4) +
  scale_color_bright(name = 'Core') +
  ylim(c(0, NA))

# fit a null model
m_null <- mvgam(formula = diatox ~ s(core, bs = 're'),
                family = Gamma(link = 'log'),
                data = d_train,
                newdata = d_test,
                chains = 4,
                burnin = 500,
                sample = 750,
                parallel = TRUE)

# unsurprisingly, the fit is horrible
plot_grid(plot(forecast(m_null), series = 1),
          plot(forecast(m_null), series = 2),
          plot(forecast(m_null), series = 3),
          plot(forecast(m_null), series = 4))

# fit a hierarchical GAM
#' `https://doi.org/10.7717/peerj.6876/fig-4`
m_gam <- mvgam(formula = diatox ~ 0,
               trend_formula = ~
                 s(year, bs = 'tp', k = 15) + # common trend across cores
                 s(year, core, bs = 'fs', k = 15), # deviations from common trend
               noncentred = TRUE,
               family = Gamma(link = 'log'),
               data = d_train,
               newdata = d_test,
               chains = 4,
               burnin = 500,
               samples = 750,
               parallel = TRUE)
summary(m_gam) # no issues with diagnostics
plot(m_gam) # likely autocorrelation at lag 1
draw(m_gam$trend_mgcv_model)

plot_grid(plot(m_gam, type = 'trend', series = 1),
          plot(m_gam, type = 'trend', series = 2),
          plot(m_gam, type = 'trend', series = 3),
          plot(m_gam, type = 'trend', series = 4))

# accuracy of forecasts isn't vary inaccurate, but it's also quite uncertain
#' NOTE: `plot(m_gam, type = 'forecast')` uses base `graphics::plot()` 
plot_grid(plot(forecast(m_gam), series = 1),
          plot(forecast(m_gam), series = 2),
          plot(forecast(m_gam), series = 3),
          plot(forecast(m_gam), series = 4))

# add AR(1) process
m_gam_ar1 <- mvgam(formula = diatox ~ 0,
                   trend_formula = ~
                     s(year, bs = 'tp', k = 15) +
                     s(year, core, bs = 'fs', k = 15),
                   trend_model = AR(1), # autoregressive 1 process
                   noncentred = TRUE,
                   family = Gamma(link = 'log'),
                   data = d_train,
                   newdata = d_test,
                   chains = 4,
                   burnin = 500,
                   sample = 750,
                   parallel = TRUE)

summary(m_gam_ar1) # no issues with diagnostics
plot(m_gam_ar1) # no appreciable autocorrelation at lag 1
draw(m_gam_ar1$trend_mgcv_model) #' fs terms are smoother than for `m_gam`

plot_grid(plot(m_gam_ar1, type = 'trend', series = 1),
          plot(m_gam_ar1, type = 'trend', series = 2),
          plot(m_gam_ar1, type = 'trend', series = 3),
          plot(m_gam_ar1, type = 'trend', series = 4))

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
if(file.exists('models/day-4-m_gam_var1.rds')) {
  m_gam_var1 <- readRDS('models/day-4-m_gam_var1.rds')
} else {
  #' TODO: add better priors. see `https://doi.org/10.3390/e19100555`
  # fits in ~ 10 minutes
  m_gam_var1 <- mvgam(formula = diatox ~ 0,
                      trend_formula = ~ # global smooth term explained by VAR(1)
                        s(year, core, bs = 'fs', k = 15),
                      trend_model = VAR(1), # vector AR(1) process
                      noncentred = FALSE, # cannot use if using for VAR() models
                      # priors = ,
                      family = Gamma(link = 'log'),
                      data = d_train,
                      newdata = d_test,
                      chains = 4,
                      burnin = 500,
                      sample = 750,
                      control = list(max_treedepth = 20, adapt_delta = 0.9),
                      parallel = TRUE)
  saveRDS(m_gam_var1, 'models/day-4-m_gam_var1.rds')
}

summary(m_gam_var1)
plot(m_gam_var1) # no appreciable issues with ACF or pACF
draw(m_gam_var1$trend_mgcv_model)

plot_grid(plot(m_gam_var1, type = 'trend', series = 1),
          plot(m_gam_var1, type = 'trend', series = 2),
          plot(m_gam_var1, type = 'trend', series = 3),
          plot(m_gam_var1, type = 'trend', series = 4))

#
plot_grid(plot(forecast(m_gam_var1), series = 1),
          plot(forecast(m_gam_var1), series = 2),
          plot(forecast(m_gam_var1), series = 3),
          plot(forecast(m_gam_var1), series = 4))

#' *ADD B SPLINE TERMS?*
# trend_formula = ~
#   s(year, bs = 'bs', k = 15) + # common trend across cores
#   s(year, core, bs = 'fs', k = 15,# deviations from common trend
#     xt = list(bs = 'tp')), # cannot be b spline
# trend_knots = list(year = c(1800, 1801, 1950, 2010)),

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
m_df <- mvgam(formula = diatox ~ s(core, bs = 're'), # assume stationarity
              use_lv = TRUE, n_lv = 4, # needs to be chosan manually
              trend_model = AR(1),
              noncentred = FALSE, # since trend formula is empty
              family = Gamma(link = 'log'),
              data = d_train,
              newdata = d_test,
              chains = 4,
              burnin = 500,
              sample = 750,
              control = list(max_treedepth = 20, adapt_delta = 0.9),
              parallel = TRUE)

summary(m_df) # no issues with diagnostics
plot(m_df) # strong but uncertain autocorrelation at lag 1
m_df$trend_mgcv_model # there is no mgcv model now

# forecasts are quite poor because they assume stationarity...
plot_grid(plot(m_df, type = 'trend', series = 1),
          plot(m_df, type = 'trend', series = 2),
          plot(m_df, type = 'trend', series = 3),
          plot(m_df, type = 'trend', series = 4))

# ... but the trends do follow the data closely, for being such a simple model
plot_grid(plot(forecast(m_df), series = 1),
          plot(forecast(m_df), series = 2),
          plot(forecast(m_df), series = 3),
          plot(forecast(m_df), series = 4))

# add a GAM term
m_gam_df <- mvgam(formula = diatox ~ 1,
                  trend_formula = ~ # global smooth term explained by DFA terms
                    s(year, core, bs = 'fs', k = 15),
                  use_lv = TRUE, n_lv = 4,
                  trend_model = AR(1),
                  noncentred = TRUE, #' since `trend_formula` is not empty
                  family = Gamma(link = 'log'),
                  data = d_train,
                  newdata = d_test,
                  chains = 4,
                  burnin = 500,
                  sample = 750,
                  control = list(max_treedepth = 20, adapt_delta = 0.9),
                  parallel = TRUE)

summary(m_gam_df) # no issues with diagnostics
plot(m_gam_df) # no appreciable autocorrelations
draw(m_gam_df$trend_mgcv_model) # terms are again quite smooth

# predictions are now tighter and more informed because they depend on
# time-varying coefficients (i.e., the GAM portion)
plot_grid(plot(m_gam_df, type = 'trend', series = 1),
          plot(m_gam_df, type = 'trend', series = 2),
          plot(m_gam_df, type = 'trend', series = 3),
          plot(m_gam_df, type = 'trend', series = 4))

# predictions are similar to previous models and still not great
plot_grid(plot(forecast(m_gam_df), series = 1),
          plot(forecast(m_gam_df), series = 2),
          plot(forecast(m_gam_df), series = 3),
          plot(forecast(m_gam_df), series = 4))

# other element of the list is a list of posterior correlations
lv_correlations(object = m_df)$mean_correlations # strong correlations if no GAM
lv_correlations(object = m_gam_df)$mean_correlations # moderate for 2, 3, 4

# adding the GAM reduced the magnitude of most coefficients because it explains
# the trends over time
abs(lv_correlations(object = m_gam_df)$mean_correlations) <
  abs(lv_correlations(object = m_df)$mean_correlations)

# but correlations are still moderately strong for (2, 3), (2, 4), and (3, 4)
lv_correlations(object = m_gam_df)$mean_correlations > 0.3

#' FIXME: scores are unreliable because the null model is listed as the best one
#' need to see which data point(s) are causing the values to be unreliable
loo_compare(m_null, m_gam, m_gam_ar1, m_gam_var1, m_df, m_gam_df)

# advantages of VAR processes:
# - can use for assessing whether one series can predict another (often called
#   "Granger causality", but it's not causality -- "precedence" is more correct)
# - can allow us to develop complex correlations across many different variables
#   across time
# - can be very unstable and data-hungry
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
m_gam_ar1_bs <- mvgam(formula = diatox ~ 0,
                      trend_formula = ~
                        core + #' `by` smooths require explicit intercepts
                        s(year, by = core, bs = 'bs', m = c(3, 2), k = 15),
                      knots = list(year = c(1800, 2010)),
                      trend_model = AR(1), # autoregressive 1 process
                      noncentred = TRUE,
                      family = Gamma(link = 'log'),
                      data = d_train,
                      newdata = d_test,
                      chains = 4,
                      burnin = 500,
                      sample = 750,
                      parallel = TRUE)

summary(m_gam_ar1_bs) # no issues with diagnostics
plot(m_gam_ar1_bs) # no appreciable autocorrelations
draw(m_gam_ar1_bs$trend_mgcv_model) # each term has its own smoothness parameter

# predictions are now tighter and more informed because they depend on
# time-varying coefficients (i.e., the GAM portion)
plot_grid(plot(m_gam_ar1_bs, type = 'trend', series = 1),
          plot(m_gam_ar1_bs, type = 'trend', series = 2),
          plot(m_gam_ar1_bs, type = 'trend', series = 3),
          plot(m_gam_ar1_bs, type = 'trend', series = 4))

# predictions are similar to previous models and still not great
plot_grid(plot(forecast(m_gam_ar1_bs), series = 1),
          plot(forecast(m_gam_ar1_bs), series = 2),
          plot(forecast(m_gam_ar1_bs), series = 3),
          plot(forecast(m_gam_ar1_bs), series = 4))

#' `m_gam_ar1` has much lower uncertainty
plot_grid(plot(forecast(m_gam_ar1), series = 1),
          plot(forecast(m_gam_ar1), series = 2),
          plot(forecast(m_gam_ar1), series = 3),
          plot(forecast(m_gam_ar1), series = 4))

# compare models based on their forecast scores
models <-
  tibble(name = factor(c('null', 'GAM only', 'GAM with AR(1)', 'GAM with VAR(1)',
                         'AR(1) and DF(2)', 'GAM with AR(1) and DF(2)',
                         'GAM with AR(1) and b-splines'),
                       levels = c('null', 'AR(1) and DF(4)', 'GAM only',
                                  'GAM with AR(1)', 'GAM with VAR(1)',
                                  'GAM with AR(1) and DF(4)',
                                  'GAM with AR(1) and b-splines')),
         model = list(m_null, m_gam, m_gam_ar1, m_gam_var1,
                      m_df, m_gam_df, m_gam_ar1_bs),
         forecast = map(model, function(.m) forecast(.m)),
         scores = map(forecast, function(.f) {
           tibble(
             year = unique(d_test$year),
             time = unique(d_test$time),
             energy_score = score(.f, score = 'energy')$all_series$score,
             variogram_score = score(.f, score = 'variogram')$all_series$score)
         }))

# lower energy and variogram scores are better: less deviation from data
models %>%
  select(name, scores) %>%
  unnest(scores) %>%
  filter(! is.na(energy_score)) %>%
  pivot_longer(c(energy_score, variogram_score), names_to = 'type',
               values_to = 'score') %>%
  mutate(type = gsub('_score', '', type)) %>%
  ggplot(aes(year, score, color = name, fill = name)) +
  facet_grid(type ~ ., scales = 'free_y') +
  geom_point() +
  geom_smooth(method = 'gam', formula = y ~ s(x), se = FALSE) +
  scale_color_bright(name = 'Model') +
  scale_fill_bright(name = 'Model')
