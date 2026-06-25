source("packages.R") # attach necessary packages

#' today's topics:
#' - introduction to time series and time series visualization
#' - ARIMA models: theory, assumptions, limitations, and applications
#' - GLMs and GAMs for ecological modelling
#' - fitting nonlinear models with `{mvgam}`
#' - dealing with irregular sampling over time

# introduction to time series and time series visualization ----
# a simple time series: daily temperature values
data("airquality")
?airquality
head(airquality) #' *NOTE:* temperature is in Fahrenheit degrees

# clean up the data format
d_temp <- airquality %>%
  rename_with(stringr::str_to_snake, everything()) %>% # convert to snake_case
  mutate(date = as_date(paste0("1973-", month, "-", day)),
         doy = yday(date),
         week_re = factor(week(date)),
         time = order(date)) %>% #' required for `trend_model` in `mvgam()`
  select(temp, month, date, doy, week_re, time) %>%
  as_tibble()
d_temp

# plot the data
p_temp <-
  ggplot(d_temp, aes(date, temp)) +
  geom_point(alpha = 0.75) +
  labs(x = NULL, y = expression(bold(paste("Temperature (\U00B0", "F)"))))
p_temp

# plot the data with a smooth model
p_temp + geom_smooth(color = "darkorange", fill = "darkorange", method = "gam",
                     formula = y ~ s(x, k = 5))

# plot the data with a wiggly model
p_temp + geom_smooth(color = "darkorange", fill = "darkorange", method = "gam",
                     formula = y ~ s(x, k = 50), n = 400)

# plot the data with an extremely wiggly model
p_temp + geom_smooth(color = "darkorange", fill = "darkorange", method = "gam",
                     formula = y ~ s(x, k = nrow(d_temp) - 1), n = 400,
                     method.args = list(gamma = 0.00001))

#' *questions:*
#' - which model is better?
#' - what counts as a change in the time series?
#' - what underlying process are we trying to estimate from the data?

# some traditional time series models and their assumptions ----
# traditional time series models often focus on discrete-time correlations
# across previous observations, e.g.: at times t -3, t - 2, and t - 1.
# example of a time series where values correlate with the previous value
d_ar <- tibble(t = 1:1e3, z = rnorm(length(t)),
               zero = NA_real_, pos = NA_real_, neg = NA_real_)
d_ar

for(i in 1:nrow(d_ar)) {
  d_ar$zero[i] <- d_ar$pos[i] <- d_ar$neg[i] <- d_ar$z[i]
  # if not first lag, add a portion of the previous value
  if(i > 1) {
    d_ar$zero[i] <- d_ar$zero[i] + d_ar$zero[i - 1] * 0.0001
    d_ar$pos[i] <- d_ar$pos[i] + d_ar$pos[i - 1] * 0.99
    d_ar$neg[i] <- d_ar$neg[i] + d_ar$neg[i - 1] * (-0.99)
  }
}
d_ar

# coefficient near 0: series is similar to random noise
ggplot(d_ar, aes(t, zero)) + geom_line(alpha = 0.5) + geom_point()

# coefficient just below 1: series struggles to stay near starting point
ggplot(d_ar, aes(t, pos)) + geom_line(alpha = 0.5) + geom_point()

# coefficient just above -1: series stays near starting point; oscillates a lot
ggplot(d_ar, aes(t, neg)) + geom_line(alpha = 0.5) + geom_point()

# diagnostics for these models are based on the Auto-Correlation Function (ACF)
# and the Partial Auto-Correlation Function (PACF), which show estimated
# correlations at different lags
layout(matrix(1:6, ncol = 2))
#' correlations between time points `t` and `t-l`
acf(d_ar$zero)  # lag-1 ACF will always be 1, the rest are negligible
acf(d_ar$pos)   # strong positive correlation with previous lags; decays w lag 
acf(d_ar$neg)   # as above, but sign flips every lag: strong lag-1 negative coef
#' correlations between time points after removing correlations of previous lags
pacf(d_ar$zero) # again, no appreciable correlations
pacf(d_ar$pos)  # strong, positive correlation at lag 1 only
pacf(d_ar$neg)  # strong, negative correlation at lag 1 only
layout(1)

#' the *ARIMA* model is a particularly common model. it is the combination of
#' three discrete-time models:
#' *AR*: AutoRegressive: observations are a function of previous values
#' *MA*: Moving Average: error terms are a linear combination of previous values
#' *I*: Integrated: AR & MA models apply to differences of consecutive values
#' we won't be covering the integrated model in detail because it only requires
#' applying the models to differences of the values. The number of differences
#' is determined by the second number in the ARIMA(a, d, m) model.
#' For example, an `ARIMA(2, 1, 3)` would have 2 AR coefficients, 1 difference,
#' and 3 MA coefficients. In this case, for values `y_t`, an `I(1)` model would
#' apply the AR and MA to `z_t = y_t - y_{t-1}`. nth-order differences remove
#' trends that are approximately nth-degree polynomials: a 1st-order difference
#' removes linear trends, while a 2nd-order difference removes parabolic trends
#' an example of a 1st-order difference (see `?AirPassengers` for more info):
layout(1:2)
plot(AirPassengers)
plot(diff(AirPassengers))
layout(1)

# some examples of AR, MA, and ARMA models
d_ts <- tibble(t = 1:1e3,
               w = rnorm(length(t)),
               ar_1 = if_else(t == 1, w, NA_real_), # AR(1) process
               ma_1 = if_else(t == 1, w, NA_real_), # MA(1) process
               arma_11 = if_else(t == 1, w, NA_real_))

for(i in 2:nrow(d_ts)) {
  d_ts$ar_1[i] <- d_ts$w[i] + d_ts$ar_1[i - 1] * 0.9
  d_ts$ma_1[i] <- d_ts$w[i] + d_ts$w[i - 1] * 0.9
  d_ts$arma_11[i] <- d_ts$w[i] + d_ts$ar_1[i - 1] * 0.9 + d_ts$w[i - 1] * (-0.5)
  rm(i)
}

#' *white noise*: `y_t` is IID Gaussian: `w_t`
ggplot(d_ts, aes(t, w)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_line()

layout(1:2)
acf(d_ts$w, ci = 0.99) # lag-0 AR will always be 1
pacf(d_ts$w, ci = 0.99) # no appreciable partial correlations
arima(d_ts$w, order = c(0, 0, 0)) %>% coef()

#' *AR(1) process*: `y_t` is a function of of `y_{t-i}` and `w_t`
ggplot(d_ts, aes(t, ar_1)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_line()

acf(d_ts$ar_1, ci = 0.99); pacf(d_ts$ar_1, ci = 0.99)
m_ar <- arima(d_ts$ar_1, order = c(1, 0, 0))
coef(m_ar)
acf(resid(m_ar)); pacf(resid(m_ar)) # residuals don't have appreciable signals

#' *MA(1) process*: errors of `y_t` are a function of `w_t` and `w_{t-1}`
ggplot(d_ts, aes(t, ma_1)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_line()

acf(d_ts$ma_1, ci = 0.99); pacf(d_ts$ma_1, ci = 0.99)
m_ma <- arima(d_ts$ma_1, order = c(0, 0, 1))
coef(m_ma)
acf(resid(m_ma)); pacf(resid(m_ma)) # residuals don't have appreciable signals

#' *ARMA(1, 1) process*: `AR(1) + MA(1)`
ggplot(d_ts, aes(t, arma_11)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_line()

acf(d_ts$arma_11, ci = 0.99); pacf(d_ts$arma_11, ci = 0.99) # AR + MA plots
m_arma <- arima(d_ts$arma_11, order = c(1, 0, 1))
coef(m_arma)
acf(resid(m_arma)); pacf(resid(m_arma)) # ARMA models can be hard to capture
layout(1)

# all three models predict that the data will be stationary and centered at the
# intercept with some random variation over time
n_preds <- 100

preds <-
  tibble(model_type = c("AR(1)", "MA(1)", "ARMA(1, 1)", "White noise"),
         model = list(m_ar, m_ma, m_arma, "White noise"),
         predictions = map(model, function(.m) {
           if(class(.m) == "character") {
             tibble(t = max(d_ts$t) + 1:n_preds,
                    pred = mean(d_ts$w),
                    se = sd(d_ts$w)) %>%
               return()
           } else {
             bind_cols(t = max(d_ts$t) + 1:n_preds,
                       predict(.m, n.ahead = n_preds) %>%
                         bind_cols()) %>%
               return()
           }
         })) %>%
  select(! model) %>%
  unnest(predictions) %>%
  mutate(lwr_90 = pred + qnorm(p = 0.05) * se,
         upr_90 = pred + qnorm(p = 0.95) * se)

d_ts %>%
  filter(t > 900) %>%
  pivot_longer(! t, names_to = "model_type") %>%
  mutate(model_type = case_when(model_type == "ar_1" ~ "AR(1)",
                                model_type == "ma_1" ~ "MA(1)",
                                model_type == "arma_11" ~ "ARMA(1, 1)",
                                model_type == "w" ~ "White noise")) %>%
  mutate(int = mean(value), .by = model_type) %>%
  ggplot() +
  facet_wrap(~ model_type) +
  geom_hline(aes(yintercept = int), lty = "dashed") +
  geom_point(aes(t, value), alpha = 0.5) +
  geom_ribbon(aes(t, ymin = lwr_90, ymax = upr_90), preds, alpha = 0.3,
              lwd = 0.5, , fill = "darkorange", color = "darkorange") +
  geom_line(aes(t, pred), preds, color = "darkorange", lwd = 1) +
  labs(x = "Timepoint", y = "Value")

#' see `?AirPassengers` for an example model that includes seasonal trends

#' *NOTE:* processes are stationary if `sum(coefs)` < 1 (not intercept)
plot_process <- function(ar = c(0), ma = c(0), n = 1e3, return_values = FALSE) {
  out <- rep(NA_real_, n) # vector of values to return
  w <- rnorm(n)
  
  for(i in 1:n) {
    #' take `length(coefs)` previous observations
    ar_indices <- i - (length(ar):1)
    ma_indices <- i - (length(ma):1)
    
    #' negative and `integer(0)` indices cause issues
    ar_indices[ar_indices <= 0] <- NA_integer_
    ma_indices[ma_indices <= 0] <- NA_integer_
    
    #' need to use `sum()` to drop the `NA`s
    out[i] <- sum(w[i], out[ar_indices] * ar, w[ma_indices] * ma,
                  na.rm = TRUE)
  }
  
  main_lab <- case_when(all(ar == 0, ma == 0) ~ "White noise",
                        all(ar != 0 & ma == 0) ~ "AR process",
                        all(ar == 0 & ma != 0) ~ "MA process",
                        all(ar != 0 & ma != 0) ~ "ARMA process",
                        .default = "Odd parameter input!")
  
  if(round(sum(ar, na.rm = TRUE), 10) >= 1) { # to rm floating point error
    main_lab <- paste(main_lab, "(non-stationary)")
  }
  
  layout(matrix(c(1, 1, 2, 3), ncol = 2, byrow = TRUE))
  plot(out, type = "l", ylab = "Value",
       main = main_lab)
  acf(out, main = "")
  pacf(out, main = "")
  layout(1)
  if(return_values) return(out)
}

plot_process(ar = 0, ma = 0) # white noise

# AR processes: stationary if coefficients sum to < 1
plot_process(ar = c(0.7, 0.2)) # stationary since 0.7 + 0.2 < 1
plot_process(ar = c(0.7, 0.6)) # non-stationary since 0.7 + 0.6 > 1
plot_process(ar = -c(0.7, 0.6)) # stationary since -0.7 - 0.6 < 1
plot_process(ar = c(0.5, 0.4, 0.3)) # non-stationary since sum > 1
plot_process(ar = c(0.4, 0.3, 0.2, 0.0999)) # stationary since sum < 1
plot_process(ar = c(0.4, 0.3, 0.2, 0.0999), n = 1e5) # stationary since sum < 1

# MA processes are always stationary: coefficients only scale the variance
plot_process(ma = c(0.9, 0.1))
plot_process(ma = c(7, 6))
plot_process(ma = c(700, 50))

# ARMA processes: stationary if AR is stationary
plot_process(ar = c(0.7, 0.2), ma = c(0.5, 0.3)) # stationary ARMA(2, 2)
plot_process(ar = c(0.7, 0.3), ma = c(0.5, 0.3)) # non-stationary ARMA(2, 2)

#' **break**

# applying the plots to the air temperature series ----
p_temp # what counts as a "change of interest"?

layout(t(1:2))
acf(d_temp$temp) #' correlation of data pairs at a times `t` and `t+Lag`
pacf(d_temp$temp) #' correl. between pairs, without previous lags' effects
# ACF decays smoothly; high pacf value at lag 1, other values are small
# AR(1) model is a good start

#' `y_t = 77.33 + 0.82 * (y_{t-1} - 77.33) + e_t`
m_ar_temp <- arima(d_temp$temp, order = c(1, 0, 0))
coef(m_ar_temp)
acf(resid(m_ar_temp)); pacf(resid(m_ar_temp)) # residuals are ok
layout(1)

# strengths and limitations of ARMA models:
# - because coefficient estimates depend on pairs of points, the models are
#   robust to missing data, whether random or even missing sections. however,
#   estimates are sensitive to sampling interval and data thinning (see below)
# - AR and MA models can help understand the properties of the data, but they
#   cannot be used to interpolate or estimate the mean trend in a time series
#   other than the intercept term
# - AR and MA models assume the models are stationary over time, so predictions
#   past the range of the data are often not useful (if not in the short term)
# - the models only work in discrete time with fixed sampling intervals

# ARMA models are sensitive to data thinning
d_thin <-
  tibble(thinning = c(1, 10, 20, 30, 40, 50, 60), # thinning intervals
         data = map(thinning, \(.t) d_ts %>% filter(1:n() %% .t == 0) %>%
                      rename(y = ar_1)),
         model = map(data, \(.d) arima(.d$y, c(1, 0, 0))),
         intercept = map_dbl(model, \(.m) coef(.m)["intercept"]),
         ar_1 = map_dbl(model, \(.m) coef(.m)["ar1"]))

d_thin %>%
  unnest(data) %>%
  ggplot() +
  facet_wrap(~ thinning) +
  geom_line(aes(t, y))

d_thin %>%
  select(thinning, ar_1, intercept) %>%
  pivot_longer(c(ar_1, intercept)) %>%
  ggplot(aes(thinning, value)) +
  facet_wrap(~ name, scales = "free_y", ncol = 1) +
  geom_point()

# data thinning is particularly problematic with animal movement data ----
#' simulate a movement path from a continuous-time stochastic process
mm <- ctmm(tau = c(180, 5) %#% "minutes", sigma = 10000, mu = c(0, 0))
track <- simulate(mm, nsim = 1, seed = 160, t = 1:1e3 %#% "minutes")
plot(track)

# calculate speed using a discrete-time approach: straight-line displacement
data.frame(track) %>%
  mutate(displacement = sqrt(x^2 + y^2),
         time_interval = t - lag(t),
         speed = displacement / time_interval) %>%
  ggplot(aes(t, speed)) +
  geom_line() +
  ylab("SLD divided by time interval (m/s)")

# calculate SLD for thinned time series
d_track <- tibble(thinning = 1:60,
                  sampling_interval = thinning * unique(diff(track$t)),
                  data = map(thinning, \(.thin) {
                    data.frame(track) %>%
                      slice(seq(1, n(), by = .thin)) %>%
                      mutate(displacement = sqrt(x^2 + y^2),
                             time_interval = t - lag(t),
                             speed = displacement / time_interval) %>%
                      return()
                  }),
                  mean_speed = map_dbl(data, \(.d) {
                    return(mean(.d$speed, na.rm = TRUE))
                  }))
d_track

# plot the estimated speed for each sampling interval
d_track %>%
  unnest(data) %>%
  filter(sampling_interval <= 540) %>%
  mutate(sampling_interval = paste(sampling_interval, "seconds") %>%
           factor(., unique(.))) %>%
  ggplot(aes(t, speed)) +
  facet_wrap(~ sampling_interval) +
  geom_line() +
  ylab("SLD divided by time interval (m/s)")

# estimated speed decreases substantially with sampling interval...
ggplot(d_track, aes(sampling_interval, mean_speed)) +
  geom_line() +
  labs(x = "Sampling interval (seconds)",
       y = "Straight-line displacement divided by time interval (m/s)")

# ... which is not surprising since we are losing data on the complexity of the
# tracks as we thin them. the issues compound if sampling intervals are
# irregular. for more info, see:
# - https://doi.org/10.1186/s40462-019-0177-1
# - https://doi.org/10.1086/675504
# - https://doi.org/10.1101/2025.07.17.665364
d_track %>%
  filter(thinning %in% c(1, 5, 10, 15, 30, 60)) %>%
  mutate(sampling_interval =
           paste(sampling_interval / 60,
                 if_else(sampling_interval == 60, "minute", "minutes")) %>%
           factor(., levels = unique(.))) %>%
  unnest(data) %>%
  ggplot(aes(x, y)) +
  coord_equal() +
  facet_wrap(~ sampling_interval) +
  geom_path(aes(x, y), d_track$data[[1]], alpha = 0.3) +
  geom_path() +
  geom_point() +
  scale_x_continuous(name = "x (meters)", breaks = (-3:3) * 100) +
  ylab("y (meters)")

# GLMs and GAMs for ecological modelling ----

#' three main parts to a GLM/GAM:
#' 1. *family* of distributions: distribution of response, Y
#' 2. *linear predictor*: sum of coefficients multiplied by predictor variables
#' 3. *link function*: connects linear predictor with parameter estimates

#' linear models are Gaussian GLMs:
#' 1. family is Gaussian with mean `mu` and variance `sigma^2`
#' 2. linear predictor is `beta_0 + x_1 * beta_1 + ...`
#' 3. link function is the identity function: `I(c) = c`: input = output

#' *fit the model to the data, not the data to the model!*
#' choose a family of distributions and link function based on:
#' 1. the possible values of the response variable
#' 2. the mean-variance relationship
#' 3. any additional considerations about the variance, such as overdispersion

?mvgam::mvgam_families #' families supported by `{mvgam}`
#' `gaussian()` for real-valued data
#' `student_t()` for heavy-tailed real-valued data
#' `lognormal()` for non-negative real-valued data
#' `Gamma()` for non-negative real-valued data
#' `betar()` for proportional data on (0,1)
#' `bernoulli()` for binary data
#' `poisson()` for count data
#' `nb()` for overdispersed count data
#' `tweedie()` for overdispersed count data (power parameter fixed at p = 1.5)
#' `binomial()` for count data with known number of trials
#' `beta_binomial()` for overdispersed count data with known number of trials
#' `nmix()` for count data with imperfect detection (unknown number of trials)
#' does not support any families that require a list of formulas
#' if you need more flexibility, see `?mgcv::family.mgcv` and
#' `?brms::family.brmsfit` for families supported by `{mgcv}` or `{brms}` but
#' not necessarily supported by `{mvgam}`.

#' *choose a link function based on the possible values for the distribution*
#' unbounded: identity; `I(-Inf, Inf) = (-Inf, Inf)`
#' `Y >= 0` or `Y > 0`: `log(0, Inf) = (-Inf, Inf)`
#' `0 <= Y <= 1`: `logit(0, 1) = log(odds(0, 1)) = log(0,Inf) = (-Inf,Inf)`
#' there are other options, but these are generally sufficient (esp. with GAMs)

#' *note:* link functions introduce two new terms:
#' - *response scale*: the original response values; e.g., (0, Inf), (0, 1)
#' - *link scale*: the transformed response values; generally (-Inf, Inf)

#' *note:* link function is applied to the *mean*, not to the data directly

#' `E(Y) = mu`
#' `g(mu) = eta = b_0 + b_1 * x_1 + b_2 * x_2`
#' `mu = g^{-1}(eta) = g^{-1}(b_0 + b_1 * x_1 + b_2 * x_2)`

# fitting models with {mvgam} ----
ggplot(d_temp, aes(date, temp)) +
  geom_point(alpha = 0.75) +
  labs(x = NULL, y = expression(bold(paste("Temperature (\U00B0", "F)"))))

#' how models are fit in `{mvgam}`: a quick intro to Bayesian modeling
#' - Bayesians view knowledge as a constantly updating set of information
#' - before fitting any models or viewing new data, we have some *prior beliefs*
#' - by observing some data, we can find the *most likely parameter values* for
#'   that specific dataset and model
#' - we can then update our prior knowledge state with the new data to produce a
#'   *posterior* knowledge state
#' - choose priors that include sensible parameters; prevent impossible
#'   parameter values, but don't force a specific hypothesis. for more info:
#' ..- `https://doi.org/10.3390/e19100555`
#' ..- `https://www.youtube.com/watch?v=ztbYkBPDOgU`
#' 
#' can get priors from `get_mvgam_priors()` without fitting a model
#' ensure to specify family and link scale so coefficients on the right scale
get_mvgam_priors(temp ~ doy, data = d_temp, family = gaussian())
#' uses `brms::get_prior()` to generate uninformative priors

#' to check whether priors are reasonable, we can use prior predictive checks:
#' - generate random coefficients by sampling the priors
#' - equivalent to checking whether our prior beliefs make sense
b1_prior_samples <- 79 + 8.9 * rt(n = 1e4, df = 3) # student_t(nu, mu, sigma)
b2_prior_samples <- 0 + 2 * rt(n = length(b1_prior_samples), df = 3)

#' prior for `(Intercept)` is for `temp` for `mean(doy)`, not when `doy = 0`
#' thus the prior for `(Intercept)` is close to `mean(temp)`
#' we will see why this is useful once we fit GAMs
mean(d_temp$temp)

# priors are quite uninformative (wide)
plot_grid(
  ggplot() +
    geom_histogram(aes(b1_prior_samples), fill = "grey", color = "black") +
    xlab("Prior for (Intercept)"),
  ggplot() +
    geom_histogram(aes(b2_prior_samples), fill = "grey", color = "black") +
    xlab("Prior for doy effect"),
  ncol = 1
)

# prior slopes are too steep
ggplot() +
  geom_abline(intercept = b1_prior_samples[1:1000],
              slope = b2_prior_samples[1:1000],
              color = "red3", alpha = 0.1) +
  geom_point(aes(doy - mean(doy), temp), d_temp)

#' can also simulate priors using `{mvgam}`
#' *NOTE:* `burnin` is ignored if `prior_simulation = TRUE`
m_temp_prior <- mvgam(temp ~ doy, #' model is `mu = b[1] + b[2] * doy`
                      family = gaussian(), #' *NOTE:* default family is Poisson
                      data = d_temp, prior_simulation = TRUE, samples = 2000,
                      silent = 2)

plot_grid(
  ggplot() +
    geom_histogram(aes(m_temp_prior$model_output@sim$samples[[1]]$`b[1]`),
                   fill = "grey", color = "black") +
    xlab("Prior for (Intercept)"),
  ggplot() +
    geom_histogram(aes(m_temp_prior$model_output@sim$samples[[1]]$`b[2]`),
                   fill = "grey", color = "black") +
    xlab("Prior for doy effect"),
  ncol = 1
)

ggplot() +
  geom_abline(intercept = m_temp_prior$model_output@sim$samples[[1]]$`b[1]`,
              slope = m_temp_prior$model_output@sim$samples[[1]]$`b[2]`,
              color = "red3", alpha = 0.1) +
  #' *NOTE:* `(Intercept)` is for `temp` for mean `doy`, not `temp` when `doy = 0`
  geom_point(aes(doy - mean(doy), temp), d_temp)

#' when fitting the model, `cmdstanr` uses Hamiltonian Monte Carlo (HMC),
#' specifically the No-U-Turn Sampler (NUTS). The method simulates a marble
#' rolling up and down surfaces of the *negative log likelihood* space by
#' leveraging Hamiltonian dynamics (i.e., physics laws of movement and energy).
#' 
#' warmup phase:
#' 0. choose a starting point based **FIX THIS** or the given initial values
#' 1. "push the marble": sample a value from the momentum distribution, `N(0, M)`
#' 2. take a series of steps based on the momentum, stopping if you make a
#'    U-turn, to avoid backtracking and sampling inefficiently
#' 3. calculate the likelihood at the starting point and the final point
#' 4. with probability alpha, accept the new candidate if the ball has lower
#'   -log(likelihood)
#' 5. repeat steps 1-5 while optimizing step length and covariance matrix `M`
#'    to achieve the target `P(acceptance)` for each warmup iteration
#' 
#' sampling phase:
#' - using the optimized step length and covariance matrix, repeat steps 1-5,
#'   but now store the values from each sample
#' - since worse choices can still be stored with some small probability,
#'   the sampling process approximates the posterior distribution efficiently
#' - using Hamiltonian dynamics and preventing U-turns results in more efficient
#'   sampling than traditional MCMC sampling
#'   
#' additionally:
#' - the target acceptance probability is the `adapt_delta` parameter
#' - the number of steps per sample is limited by the `max_treedepth` parameter,
#'   which is the max doubling of steps if no U-turns occur
#' - numbers of samples are controlled by the `burnin` and `samples` arguments
#' 
#' for more info:
#' - `https://mc-stan.org/docs/reference-manual/execution.html#random-initial-values`
#' - `https://arogozhnikov.github.io/2016/12/19/markov_chain_monte_carlo.html`
#' - `https://discourse.mc-stan.org/t/the-role-of-max-treedepth-in-no-u-turn/24155/2`
#' - `https://jwmi.github.io/BMB/18-Hamiltonian-Monte-Carlo-and-NUTS.pdf`
#' - `https://arxiv.org/abs/1701.02434`
#' - `https://mc-stan.org/learn-stan/diagnostics-warnings.html`

# fitting a simple GLM (fits in < 1 second)
m_temp <- mvgam(temp ~ doy,
                family = gaussian(),
                data = d_temp, samples = 1000, burnin = 1000)

# posterior distributions are much narrower than priors, so the data were more
# informative than the prior
layout(matrix(1:4, ncol = 2, byrow = TRUE))
hist(b1_prior_samples, main = "Prior for intercept")
hist(b2_prior_samples, main = "Prior for slope")
hist(m_temp$model_output@sim$samples[[1]]$`b[1]`, main = "Est. intercept",
     xlim = range(b1_prior_samples))
hist(m_temp$model_output@sim$samples[[1]]$`b[2]`, main = "Est. slope",
     xlim = range(b2_prior_samples))
layout(1)

layout(1:2)
hist(m_temp$model_output@sim$samples[[1]]$`b[1]`, main = "Est. intercept")
hist(m_temp$model_output@sim$samples[[1]]$`b[2]`, main = "Est. slope")
layout(1)

ggplot() +
  geom_abline(intercept = m_temp$model_output@sim$samples[[1]]$`b[1]`,
              slope = m_temp$model_output@sim$samples[[1]]$`b[2]`,
              color = "red3", alpha = 0.1) +
  #' *NOTE:* `(Intercept)` is now `temp` when `doy = 0` because no smooth terms
  geom_point(aes(doy, temp), d_temp)

summary(m_temp)

#' **break**

# modeling nonlinear trends ----
#' fitting a GLM with a polynomial term (note: still technically linear terms)
#' The terms can't be functions of each other, so we need to add columns of the
#' polynomial that are independent of each other (i.e., orthogonal) to avoid
#' complete collinearity and non-identifiability issues when fitting.
d_temp <- d_temp %>%
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
mcmc_plot(m_temp_ar, type = "trace_highlight", highlight = 2)
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
