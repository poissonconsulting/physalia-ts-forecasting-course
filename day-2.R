source("packages.R") # attach necessary packages
source("gaussian-process-functions.R") # for plotting GP covariance function

#' yesterday's topics:
#' - ARIMA models:
#'   - AR if *data* are correlated through time
#'   - I if data need to be detrended by taking differences
#'   - MA if *errors* are correlated through time
#'  - what counts as a change depends on the scale of reference

#' today's topics:
#' - modelling nonlinear trends
#' - state space models
#' - fitting SSMs in `{mvgam}`
#' - dealing with continuous correlations over time (CAR processes)
#' - gaussian processes
#' - dynamic coefficient models (if time permits)

# modeling nonlinear trends ----
#' fitting a GLM with a polynomial term (note: still technically linear terms)
#' The terms can't be functions of each other, so we need to add columns of the
#' polynomial that are independent of each other (i.e., orthogonal) to avoid
#' complete collinearity and non-identifiability issues when fitting.
d_temp <- airquality %>%
  rename_with(stringr::str_to_snake, everything()) %>% # convert to snake_case
  mutate(date = as_date(paste0("1973-", month, "-", day)),
         doy = yday(date),
         week_re = factor(week(date)),
         time = order(date)) %>% #' required for `trend_model` in `mvgam()`
  select(temp, month, date, doy, week_re, time) %>%
  as_tibble() %>%
  bind_cols(.,
            poly(.$doy, degree = 3) %>%
              as.data.frame() %>%
              rename(doy_1 = 1, doy_2 = 2, doy_3 = 3))

#' `mu = b[1] + b[2] * doy + b[3] * doy^2 + b[4] * doy^3`
m_temp_poly_prior <- mvgam(temp ~ doy_1 + doy_2 + doy_3,
                           family = gaussian(),
                           data = d_temp, prior_simulation = TRUE,
                           samples = 2000, silent = 2)

#' prior for `(Intercept)` is for `temp` for `mean(doy)`, not when `doy = 0`
m_temp_poly_prior$model_output@sim$samples[[1]] %>%
  select(matches("b\\[[1-4]\\]")) %>% # select coefficients only
  pivot_longer(everything(), values_to = "sample", names_to = "coef") %>%
  ggplot() +
  facet_wrap(~ coef, scales = "free") +
  geom_histogram(aes(sample), fill = "grey", color = "black", bins = 100) +
  labs(x = "Prior", y = "Count")

# prior predictive samples for the polynomial model
m_temp_poly_prior$model_output@sim$samples[[1]] %>%
  select(matches("b\\[[1-4]\\]")) %>%
  as_tibble() %>%
  mutate(preds = map(1:n(), function(i) {
    mutate(d_temp,
           sample_id = i,
           y = `b[1]`[i] + `b[2]`[i] * doy_1 + `b[3]`[i] * doy_2 +
             `b[4]`[i] * doy_3)
  })) %>%
  unnest(preds) %>%
  ggplot() +
  geom_line(aes(x = doy - mean(doy), y = y, group = sample_id),
            color = "red3", alpha = 0.1) +
  geom_point(aes(doy - mean(doy), temp), d_temp)

# prior sample curves are too flat because x covariates have been rescaled to
# be quite close to 0:
range(d_temp$doy_1); range(d_temp$doy_2); range(d_temp$doy_3)

# the model still fits quickly because it's relatively simple
# but more complex models may need more careful choices of samples
m_temp_poly <- mvgam(temp ~ doy_1 + doy_2 + doy_3,
                     family = gaussian(),
                     data = d_temp, samples = 1000, burnin = 1000)

#' since `{mvgam}` fits Bayesian models with `Stan`, we should check that all
#' chains converged properly: check `Rhat`, `n_eff`, and Stan MCMC diagnostics
#' each chain is one of the "paths" the model took to estimate the parameters
#' the model will spend more time near the better parameter estimates and less
#' near unlikely values, which estimates the posterior distribution.
#' ideally, all chains should be similar and indistinguishable.
#' for more info on HMC sampling: https://www.youtube.com/watch?v=a-wydhEuAm0s
summary(m_temp_poly)
mcmc_plot(m_temp_poly, type = "trace", variable = rownames(coef(m_temp_poly)))
plot(m_temp_poly, type = "residuals") #' can also use `plot_mvgam_resids()`

# can add predictions to the data...
d_temp %>%
  bind_cols(predict(m_temp_poly, type = "response") %>%
              as.data.frame() %>%
              rename(est_poly = Estimate,
                     se_poly = Est.Error,
                     q2.5_poly = Q2.5,
                     q97.5_poly = Q97.5))

# ... but there are also many useful built-in functions
plot(hindcast(m_temp_poly, type = "response")) # uncertainty in Y
plot(hindcast(m_temp_poly, type = "link")) # uncertainty in mu on link scale
plot(hindcast(m_temp_poly, type = "expected")) # uncertainty in mu
#' *NOTE:* link and response scale are the same for identity link functions
plot(m_temp_poly, type = "smooths") # only works with GAMs (see below)

# how can we model data that have irregular sampling over time? ----
# with irregular sampling, it may help to focus on rates of change and trends
# over time rather than changes over steps in discrete time
#' *NOTE:* many of the `{mvgam}` plots assume discrete-time sampling, so the
#' missing observations should be `NA` rather than missing the full row, as
#' long as none of the predictors have `NA` values.
d_temp_missing <- d_temp %>%
  mutate(time = 1:n()) %>%
  mutate(temp = if_else(month(date) == 6, NA_real_, temp),
         temp = if_else(time %in% sample(time, size = n() / 2),
                        NA_real_, temp)) %>%
  arrange(date)
d_temp_missing #' note the `NA`s in the `temp` column

ggplot(d_temp_missing, aes(date, temp)) +
  geom_point(alpha = 0.75) +
  labs(x = NULL, y = expression(bold(paste("Temperature (\U00B0", "F)"))))

#' fitting a *GAM* with a smooth term of `doy`
#' the smooth term is created using the `s()` function
#' greater model complexity requires a bit more sampling and burnin
m_temp_gam <- mvgam(temp ~ s(doy, k = 10, bs = "cr"),
                    family = gaussian(), data = d_temp_missing,
                    parallel = TRUE, burnin = 1e3, samples = 1e3,
                    #' tune parameters to improve sampling:
                    #' - `max_treedepth`: max n of binary choices when sampling
                    #' - `adapt_delta`: target average proposal acceptance prob. 
                    control = list(max_treedepth = 10, adapt_delta = 0.99))

#' `summary()` looks a bit different from the one for the GLM:
#' - `s(doy)` has `k - 1` coefficients
#' - each coefficient is multiplied by the respective basis function
summary(m_temp_gam)
coef(m_temp_gam$mgcv_model) # model coefficients
coef(m_temp_gam$mgcv_model)[-1] #' check `s(doy)` terms only
mcmc_plot(m_temp_gam, type = "trace", variable = "doy", regex = TRUE)

#' `mcmc_plot()` only plots intercept and rho terms by default
mcmc_plot(m_temp_gam, type = "trace", variable = ".", regex = TRUE)

#' visualize the cubic basis
draw(basis(s(doy, bs = "cr"), data = d_temp_missing)) # default cubic basis
draw(basis(s(doy, bs = "cr"), data = d_temp_missing, coefficients = 1:9,
           constraints = TRUE))
draw(basis(s(doy, bs = "cr"), data = d_temp_missing, # default basis * coefs
           coefficients = coef(m_temp_gam$mgcv_model)[-1], 
           constraints = TRUE))
draw(basis(m_temp_gam$mgcv_model)) & # fitted cubic basis
  geom_line(aes(doy, .fitted - coef(m_temp_gam$mgcv_model)[1]),
            fitted_values(m_temp_gam$mgcv_model), lwd = 1,
            inherit.aes = FALSE)

plot(hindcast(m_temp_gam, type = "response")) # uncertainty in Y
plot(hindcast(m_temp_gam, type = "link")) # uncertainty in mu on link scale
plot(hindcast(m_temp_gam, type = "expected")) # identity link: same as above

#' smooths are centered at 0
#' this is why the prior for the intercepts is for the average covariate values!
plot(m_temp_gam, type = "smooths") #' model terms; == `plot_mvgam_smooth()`

plot(m_temp_gam, type = "residuals") #' model diagnostics; `plot_mvgam_resids()`

# temporal random effects and temporal residual correlation structures ----
# fit a GAM with a random effect of week
ggplot(d_temp_missing) + geom_point(aes(week_re, temp), alpha = 0.3)

m_temp_re <- mvgam(temp ~ s(week_re, bs = "re"), family = gaussian(),
                   data = d_temp_missing, parallel = TRUE, silent = 2)
summary(m_temp_re)
draw(m_temp_re$mgcv_model) # each week has a coefficient

plot(hindcast(m_temp_re, type = "response")) # response scale; uncertainty in Y

#' *Q:* how do we choose the window width? is 2 weeks or 10 days better than 7?

#' smooth terms in GAMs can be thought of as a continuous version of these
#' discrete-time random effects. The random effects in the GAMs are the basis
#' coefficients
plot_grid(
  draw(basis(s(doy, bs = "cr"), data = d_temp_missing)), # default cubic basis
  draw(basis(m_temp_gam$mgcv_model), residuals = TRUE) + # fitted cubic basis
    geom_line(aes(doy, .fitted - coef(m_temp_gam$mgcv_model)[1]),
              fitted_values(m_temp_gam$mgcv_model), lwd = 1,
              inherit.aes = FALSE),
  ncol = 1)

#' smooth terms are better at dealing with gaps and irregular sampling. they
#' also don't require choosing a window size, but you do need to choose `k`.
ggplot() +
  geom_line(aes(doy, Estimate, color = "s(doy)"),
            bind_cols(d_temp, predict(m_temp_gam, d_temp)),
            lwd = 1, inherit.aes = FALSE) +
  geom_line(aes(doy, Estimate, color = "s(week, bs = \"re\")"),
            bind_cols(d_temp, predict(m_temp_re, d_temp)),
            lwd = 1, inherit.aes = FALSE) +
  geom_point(aes(doy, temp), d_temp_missing, alpha = 0.3, inherit.aes = FALSE) +
  scale_color_highcontrast(name = "Model") +
  theme(legend.position = "top")

# clear autocorrelation at lag 1
plot(m_temp_gam)

#' add an `AR(1)` component
m_temp_ar <- mvgam(formula = temp ~ s(doy, k = 10, bs = "cr"),
                   trend_model = AR(1),
                   family = gaussian(), data = d_temp_missing, noncentred = TRUE,
                   parallel = TRUE, burnin = 500, samples = 1000,
                   control = list(adapt_delta = 0.99))

# chains are occasionally not well-mixed even if Rhat is near 1
# there seem to be two alternative fits the model is trying to decide between
mcmc_plot(m_temp_ar, type = "trace")
# chains are not well mixed: conflicting coefficients
mcmc_plot(m_temp_ar, type = "trace_highlight", highlight = 3)
summary(m_temp_ar) # Rhat can be deceiving
plot(m_temp_ar) # no appreciable correlations at any lags
# uncertainty at greater lags is because of missing data

# compare the autocorrelation to that of the previous model
plot(m_temp_gam)

# compare fits
plot_mvgam_fc(m_temp_gam)

# can't place plots in grids because they use a mix of ggplot and base plot
# additionally, some plots functions don't return the plots (can't be assigned)
plot_mvgam_fc(m_temp_ar)
plot(m_temp_ar, type = "smooth")
plot(m_temp_ar, type = "trend")

#' *break*

# SSMs ----
#' *State space models*
#' - process model:                    `mu_proc = b0 + b1 * x1 + ...`
#' - process output (states):          `Y_proc ~ MVN(mu_proc, s_proc)`
#' - process observations at time `t`: `O_t ~ MVN(Y_proc, s_obs)`
#' 
#' but estimating the process model requires us to work backwards:
#' - model the space of possible states (i.e., outcomes, responses)
#' - process observations at time `t`: `O_t ~ MVN(Y_proc, s_obs)`
#' - process output (states):          `Y_proc ~ MVN(mu_proc, s_proc)`
#' - process model:                    `mu_proc = b0 + b1 * x1 + ...`
#' uncertainty needs to be propagated accordingly across each step

# example with count data: global number of international air passengers ----
data("AirPassengers")

AirPassengers
class(AirPassengers)

air_passengers <-
  tibble(time = 1:length(AirPassengers),
         dec_date = time(AirPassengers) %>% as.numeric(),
         year = floor(dec_date),
         month = round((dec_date - year) * 12) + 1, # to ensure January is 1
         passengers = as.numeric(AirPassengers)) # in thousands
air_passengers

# visualize time as a single line
ggplot(air_passengers, aes(dec_date, passengers)) +
  geom_line() +
  labs(x = "Year CE", y = "International airline passengers (thousands)")

# visualize time as a repeating cycle
ggplot(air_passengers, aes(month, passengers, group = year, color = year)) +
  facet_wrap(~ year) +
  geom_line() +
  scale_x_continuous("Month", expand = c(0, 0)) +
  scale_y_continuous("International airline passengers (thousands)")

# visualize time as a surface
ggplot(air_passengers, aes(year, month, fill = passengers)) +
  geom_raster() +
  scale_x_continuous("Year CE", expand = c(0, 0)) +
  scale_y_continuous("Month", expand = c(0, 0), breaks = 1:12,
                     labels = month.name) +
  scale_fill_bam(name = "International airline passengers (thousands)") +
  theme(legend.position = "top")

# modeling ----
# split data into training and testing sets
ggplot(air_passengers, aes(dec_date, passengers, lty = year < 1958)) +
  geom_line() +
  geom_vline(xintercept = 1958, lty = "dashed") +
  labs(x = "Year CE", y = "International airline passengers (thousands)") +
  scale_linetype_manual("Dataset", values = c(3, 1), labels = c("Test", "Train"))

data_train <- filter(air_passengers, year < 1958)
data_test <- filter(air_passengers, year >= 1958)

#' fit a simple gam with `dec_date`
#' *NOTE:* fitting times are a lot slower if `passengers` is multiplied by 1000
m_gam <- mvgam(formula = passengers ~ s(dec_date, k = 30),
               family = poisson(link = "log"),
               data = data_train,
               newdata = data_test,
               chains = 4,
               burnin = 500,
               samples = 500,
               control = list(max_treedepth = 12, adapt_delta = 0.9),
               parallel = TRUE)

# diagnostics look ok
mcmc_plot(m_gam, type = "trace", variable = ".", regex = TRUE)
summary(m_gam)
plot(m_gam) # but there's clear unaccounted autocorrelation!

# predictions for the test dataset are quite good
plot(hindcast(m_gam))

# predictions for the test dataset are quite poor
plot(m_gam, type = "forecast")

# basis information only extends to the range of the data, so predictions past
# the last time only continue the trend at the last timestamp
draw(basis(m_gam$mgcv_model))

#' reduce model complexity to improve prediction accuracy
m_gam_smooth <- mvgam(formula = passengers ~ s(dec_date, k = 10),
                      family = poisson(),
                      data = data_train,
                      newdata = data_test,
                      chains = 4,
                      burnin = 500,
                      samples = 500,
                      control = list(adapt_delta = 0.9),
                      parallel = TRUE,
                      silent = 2)

summary(m_gam_smooth)

# predictions are better, but the model is missing the seasonality
plot(hindcast(m_gam_smooth))
plot(m_gam_smooth, type = "forecast")
draw(basis(m_gam_smooth$mgcv_model)) # simpler bases
plot(m_gam_smooth) # autocorrelation is even worse

#' we can improve the predictions somewhat by extending the basis, but this
#' still depends on the model and the complexity of the trends...
#' for more info, see:
#' https://fromthebottomoftheheap.net/2020/06/03/extrapolating-with-gams/
#' fit a GAM with a cubic B spline whose curvature and slope are penalized
#' note that `k` is quite high now
m_gam_bs <- mvgam(formula = passengers ~ s(dec_date, k = 50, bs = "bs",
                                           m = c(3, 2, 1)),
                  knots = list(dec_date = c(min(data_train$dec_date),
                                            max(data_test$dec_date))),
                  family = poisson(),
                  data = data_train,
                  newdata = data_test,
                  chains = 4,
                  burnin = 500,
                  samples = 500,
                  control = list(adapt_delta = 0.9),
                  parallel = TRUE)

#' warnings indicate there are no data for a portion of the spline, as expected
draw(basis(m_gam_bs$mgcv_model, data = air_passengers)) +  # model bases
  geom_rug(aes(dec_date), data_train, inherit.aes = FALSE) # rug plot of data

summary(m_gam_bs)

# forecasts for the test data are better, but still not good
# the model doesn't know there are seasonal cycles. it only understands that
# the values go up and down, but not *why* they do.
layout(1:3)
plot(m_gam, type = "forecast")
plot(m_gam_smooth, type = "forecast") # lower DRPS implies better forecasts
plot(m_gam_bs, type = "forecast")
layout(1)

#' more rigid models tend to extrapolate better, so we can decompose the trend
#' into `doy` and `year` trends to allow the model to learn the seasonal cycles
m_gam_month <- mvgam(passengers ~
                       s(year, k = 9, bs = "bs", # k must be <= unique(year)
                         m = c(3, 2, 1)) +
                       s(month, k = 10, bs = "cc"),
                     knots = list(year = c(range(air_passengers$year)),
                                  month = c(0.5, 12.5)), # ensures smooth cycles
                     family = poisson(),
                     data = data_train,
                     newdata = data_test,
                     chains = 4,
                     burnin = 500,
                     samples = 500,
                     control = list(adapt_delta = 0.9),
                     parallel = TRUE,
                     silent = 2)

plot(hindcast(m_gam_month)) # predicts similar oscillations over the years
plot(m_gam_month, type = "smooths") #' trends decomposed into `year` and `month`

#' value, slope, and curvature match at `month = 12.5 = 0.5`
draw(m_gam_month, select = "s(month)",
     data = tibble(month = seq(0, 24, by = 0.001), year = 0))

layout(matrix(1:4, ncol = 2))
plot(m_gam, type = "forecast")
plot(m_gam_smooth, type = "forecast")
plot(m_gam_bs, type = "forecast")
plot(m_gam_month, type = "forecast") # best model, but under-predicts a bit
layout(1)

# assuming the seasonal cycle repeats across years allows us to reduce the
# extrapolation to only be an extrapolation across years but not across months

# autocorrelation present at lags 1 and 12, but much better than previous models
plot(m_gam_month)

# state space models (SSMs) separate the true trend from observational trend by
# assuming that:
# (1) true values of the response are autocorrelated (e.g., y_t will likely be
# similar to y_{t-1}), and
# (2) residuals are independent once we account for the model (including
# previous values of y)
# for more info, see:
# - Auger-Méthé et al. (2021; https://doi.org/10.1002/ecm.1470 )
# - https://nicholasjclark.github.io/mvgam/articles/trend_formulas.html

# fit a GAM with an AR(1) process
#' the `m_gam_ar` model is:
#' `passengers ∼ Poisson(λ_t)`       # observation
#' `log(λ_t) = l_t + z_t`            # trend of mean obs on log scale
#' `l_t = b_0 + s(year) + s(month)`  # latent process trend varies w time
#' `z_t ∼ Normal(z_{t−1} * a, σ)`    # latent stochastic component on log scale
#' where `a` is the coefficient of the `AR(1)` process
m_gam_ar <- mvgam(formula = passengers ~ 0,
                  trend_formula = ~
                    s(year, k = 9, bs = "tp") +
                    s(month, k = 10, bs = "cc"),
                  trend_model = AR(p = 1), # AR(1) model
                  noncentred = TRUE, # use a noncentered AR(1) model
                  knots = list(month = c(0.5, 12.5)),
                  family = poisson(),
                  data = data_train,
                  newdata = data_test,
                  chains = 4,
                  burnin = 750,
                  samples = 500,
                  parallel = TRUE,
                  silent = 2)

summary(m_gam_ar) # diagnostics are ok

# there's still some autocorrelation at lag 12, but coef is quite uncertain
# intervals are 50%, 60%, and 90%, estimated from samples
plot(m_gam_ar)
plot(forecast(m_gam_ar))
plot(hindcast(m_gam_ar))

# add a term of lag-12 passengers
#' *NOTE:* we need to drop rows with `NA` lagged y values because the model
#' can't use them. specifying the AR process in the `trends_model` argument
#' avoids this issue. This is why we shouldn't just add the lagged values as
#' predictors in a GAM.
data_train_12 <- mutate(data_train, lag_12_passengers = lag(passengers, 12)) %>%
  filter(! is.na(lag_12_passengers))
data_test_12 <- air_passengers %>%
  mutate(lag_12_passengers = lag(passengers, 12)) %>%
  filter(year >= 1958, ! is.na(lag_12_passengers))

# add a term to crudely account for a 12-month autocorrelation
# fits in ~20 seconds
m_gam_ar_12 <- mvgam(formula = passengers ~ 0,
                     trend_formula = ~
                       s(year, k = 8, bs = "tp") +
                       s(month, k = 10, bs = "cc") +
                       log(lag_12_passengers),
                     trend_model = AR(p = 1),
                     noncentred = TRUE,
                     knots = list(month = c(0.5, 12.5)),
                     family = poisson(),
                     data = data_train_12,
                     newdata = data_test_12,
                     chains = 4,
                     burnin = 750,
                     samples = 500,
                     parallel = TRUE)

summary(m_gam_ar_12)
mcmc_plot(m_gam_ar_12, type = "trace", variable = ".", regex = TRUE)
coef(m_gam_ar_12$trend_mgcv_model)["log(lag_12_passengers)"]
plot_grid(plot(m_gam_ar), # see values at lag 12 for ACF and pACF
          plot(m_gam_ar_12)) # values at lag 12 are smaller, but lag 1 is larger

#' how does `{mvgam}` handle many missing data?
#' `{mvgam}` does not drop rows with `NA` response values, so it keeps track of
#' which values are temporally adjacent. for a comparison with `{brms}` see:
#' https://github.com/nicholasjclark/physalia-forecasting-course/blob/main/day2/tutorial_2_physalia.html
#' other advantages of using `{mvgam}` over `{brms}` include:
#' - `{mvgam}` allows each time series to have different AR1 parameters
#' - `{mvgam}` can model the correlations among errors of each time series
#' - `{mvgam}` can fit the dynamic processes using State-Space Models (SSMs)
data_train_missing <- data_train %>%
  mutate(passengers = if_else(1:n() %in% sample(1:n(), n() * 0.9), NA_real_,
                              passengers))
data_train_missing
ggplot(data_train_missing, aes(dec_date, passengers)) +
  geom_line() +
  geom_point() +
  labs(x = "Year CE", y = "International airline passengers (thousands)") +
  scale_linetype_manual("Dataset", values = c(3, 1), labels = c("Test", "Train"))

m_gam_ar_missing <-
  mvgam(formula = passengers ~ 0,
        trend_formula = ~
          s(year, k = 9, bs = "tp") +
          s(month, k = 10, bs = "cc"),
        trend_model = AR(p = 1),
        knots = list(month = c(0.5, 12.5)),
        family = poisson(),
        data = data_train_missing,
        newdata = data_test,
        chains = 4,
        burnin = 750,
        samples = 500,
        parallel = TRUE,
        silent = 2)

plot(m_gam_ar_missing, type = "forecast")
#' `{brms}` would assume that each observation follows the previous row!

# compare to a simple GAM (i.e., without AR(1) term)
m_gam_missing <-
  mvgam(formula = passengers ~ 0,
        trend_formula = ~
          s(year, k = 9, bs = "tp") +
          s(month, k = 10, bs = "cc"),
        knots = list(month = c(0.5, 12.5)),
        family = poisson(),
        data = data_train_missing,
        newdata = data_test,
        chains = 4,
        burnin = 500,
        samples = 500,
        parallel = TRUE,
        silent = 2)

# AR GAM has much more uncertainty
plot_grid(
  plot(forecast(m_gam_missing)) +
    geom_point(aes(time, passengers), air_passengers, shape = 4, size = 0.75) +
    geom_point(aes(time, passengers), data_train_missing, na.rm = TRUE)+
    ylim(c(0, 1e3)) +
    ggtitle("SSM GAM with uncorrelated error"),
  plot(forecast(m_gam_ar_missing)) +
    geom_point(aes(time, passengers), air_passengers, shape = 4, size = 0.75) +
    geom_point(aes(time, passengers), data_train_missing, na.rm = TRUE)+
    ylim(c(0, 1e3)) +
    ggtitle("SSM GAM with AR(1) trend"),
  ncol = 1)

#' **break**

# smooth correlations over time ----
# continuous auto-regressive (CAR) processes
#' data from: Gushulak et al. (2023); https://doi.org/10.1111/fwb.14192
#' study site: https://onlinelibrary.wiley.com/cms/asset/a02d5fe1-044e-4dd2-b50c-ff33f328d953/fwb14192-fig-0001-m.jpg
pigments <- read.xlsx("https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%201%20April%202014.xlsx") %>%
  as_tibble() %>%
  rename_with(stringr::str_to_snake, everything()) %>%
  select(mid_depth_cm, year, diatox, percentn) %>%
  rename(percent_n = percentn) %>%
  # collapse two duplicate years into the same row
  summarize(mid_depth_cm = mean(mid_depth_cm, na.rm = TRUE),
            year = mean(year, na.rm = TRUE),
            diatox = mean(diatox, na.rm = TRUE),
            percent_n = mean(percent_n),
            .by = year) %>%
  arrange(year) %>%
  mutate(interval = year - lag(year),
         series = factor("core 1")) %>%
  #' add a column of time for the `CAR(1)` process
  #' can use `CAR(1)` if times are not integers (equivalent otherwise)
  #' requires a consecutive series of `time` values
  right_join(tibble(year = seq(min(.$year), max(.$year), by = 1)),
             by = "year") %>%
  mutate(time = year,
         series = factor("core 1")) %>%
  arrange(time)

# sediment age decreases nonlinearly with sample depth
ggplot(pigments, aes(year, mid_depth_cm)) +
  geom_point(alpha = 0.75) +
  geom_path() +
  xlab("Year CE") +
  scale_y_reverse("Sample depth (cm)") +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, lwd = 2,
           arrow = arrow(length = unit(0.5, "cm"), ends = "last",
                         type = "closed"), color = "darkorange") +
  # Prevents clipping so left side of arrow is visible
  coord_cartesian(clip = "off") +
  theme(axis.line.y = element_line(color = "darkorange"))

# time intervals vary substantially across years
ggplot(pigments, aes(interval, mid_depth_cm)) +
  geom_path() +
  geom_point(alpha = 0.75) +
  xlab("Time interval between samples (years)") +
  scale_y_reverse("Sample depth (cm)") +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf,
           arrow = arrow(length = unit(0.5, "cm"), ends = "last",
                         type = "closed")) +
  coord_cartesian(clip = "off")

# plot an example time series (diatoxantin is a pigment produced by diatoms)
# diatoms are glass-like algae: https://en.wikipedia.org/wiki/Diatom
lab_diatox <- expression(bold(Diatoxanthin~concentration~"(nmol"~g^{"-1"}~"C)"))

ggplot(pigments, aes(year, diatox)) +
  geom_path() +
  geom_point(alpha = 0.75) +
  labs(x = "Year CE", y = lab_diatox)

#' fit a basic GAM
#' using observation formula because we can't separate it from the latent trend
#' technically, we could add `s(year)` to both formulas, since the values
#' changed over the years (trend process) and the pigment possibly decayed over
#' time (observation process), but we can't disentangle the two with only one
#' time series
m_diatox_0 <- mvgam(formula = diatox ~ s(year, k = 20),
                    family = Gamma(link = "log"),
                    data = pigments,
                    chains = 4,
                    burnin = 500,
                    samples = 500,
                    parallel = TRUE,
                    silent = 2)

plot(m_diatox_0, type = "series") # time series not accounting for the model
plot(m_diatox_0) # ACF is much lower once the model accounts for time
summary(m_diatox_0)
mcmc_plot(m_diatox_0, type = "trace", variable = ".", regex = TRUE)
draw(m_diatox_0$mgcv_model, n = 200)

#' add a `CAR(1)` term to account for continuous-time autocorrelation
#' the model is:
#' `diatox ∼ Gamma(mu_t, theta)`        # observation
#' `log(mu_t) = s(year) + l_t + z_t`    # trend of mean obs on log scale
#' `l_t = 0`                            # latent process trend varies w time
#' `z_t ∼ Normal(z_{t−dt} * a^{dt}, σ)` # latent stochastic component
#' where `0 < a < 1` is the coef of the `CAR(1)` process for time difference `dt`
#' 
#' note: when `dt = 1`, the model becomes a simple `AR(1)` process:
#' `E(z_t) = z_{t−dt} * a^{dt} = z_{t−1} * a^1 = z_{t−1} * a`
#' 
#' note: when `dt = 0`, `z_t` becomes equal to itself:
#' `E(z_t) = z_{t−dt} * a^{dt} = z_{t-0} * a^0 = z_{t} * 1`
#' 
#' note: as `dt` becomes large, `z_t` becomes independent from `z_{t-dt}`
#' `E(z_t) = z_{t−dt} * a^{Inf} = z_{t-dt} * 0 = z_{t−Inf} * 0`
m_diatox_car <- mvgam(formula = diatox ~ s(year, k = 6),
                      trend_formula = ~ 0,
                      trend_model = CAR(1),
                      noncentred = TRUE,
                      family = Gamma(link = "log"),
                      data = pigments,
                      chains = 4,
                      burnin = 500,
                      samples = 500,
                      control = list(adapt_delta = 0.95),
                      parallel = TRUE)

plot(m_diatox_car, type = "residuals") # diagnostics look great
summary(m_diatox_car) # summary looks good
mcmc_plot(m_diatox_car, type = "trace", variable = ".", regex = TRUE) # good

plot_predictions(m_diatox_car, "year") #' `s(year)` is very smooth...
plot(hindcast(m_diatox_car)) # predictions with data points are not smooth...
