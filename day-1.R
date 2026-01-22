source('packages.R') # attach necessary packages

## Introduction to time series and time series visualization ----
## a simple time series: daily temperature values
data('airquality')
?airquality
head(airquality) #' *NOTE:* temperature is in Fahrenheit degrees

## clean up the data format
d_temp <- airquality %>%
  janitor::clean_names() %>% # convert to snake_case
  mutate(date = as_date(paste0('1973-', month, '-', day)),
         doy = yday(date),
         week_re = factor(week(date)),
         time = order(date)) %>% #' required for `trend_model` in `mvgam()`
  select(temp, month, date, doy, week_re, time) %>%
  as_tibble()
d_temp

## plot the data
p_temp <-
  ggplot(d_temp, aes(date, temp)) +
  geom_point(alpha = 0.75) +
  labs(x = NULL, y = expression(bold(paste('Temperature (\U00B0', 'C)'))))
p_temp

## plot the data with a smooth model
p_temp + geom_smooth(color = 'darkorange', fill = 'darkorange', method = 'gam',
                     formula = y ~ s(x, k = 5))

## plot the data with a wiggly model
p_temp + geom_smooth(color = 'darkorange', fill = 'darkorange', method = 'gam',
                     formula = y ~ s(x, k = 50), n = 400)

## plot the data with an extremely wiggly model
p_temp + geom_smooth(color = 'darkorange', fill = 'darkorange', method = 'gam',
                     formula = y ~ s(x, k = nrow(d_temp) - 1), n = 400,
                     method.args = list(gamma = 0.00001))

#' *questions:*
#' - which model is better?
#' - what counts as a change in the time series?
#' - what underlying process are we trying to estimate from the data?

## some traditional time series models and their assumptions ----
## traditional time series models often focus on discrete-time correlations
## across previous observations, e.g.: at times t -3, t - 2, and t - 1.
## diagnostics for these models are based on the Auto-Correlation Function (ACF)
## and the Partial Auto-Correlation Function (PACF)
## the ARIMA model is a particularly common model. it is the combination of
## three discrete-time models:
## AR (AutoRegressive): observations are a linear combination of previous values
## I (Integrated): AR & MA models apply to differenced consecutive observations
## MA (Moving Average): error terms are a linear combination of previous values
## we won't be covering the integrated step in detail because it only requires
## applying the models to differences of the values. The order of the
## differences is determined by the second number in the ARIMA(a, d, m) model.
## For example, an ARIMA(2, 1, 3) would have 2 AR coefficients, 1 difference,
## and 3 MA coefficients. nth-order differences remove trends that are
## approximately nth-degree polynomials: a 1st-order difference removes linear
## trends, while a 2nd-order difference removes parabolic trends. an example of
## a 1st-order difference is:
air_passengers <- bind_cols(t = time(AirPassengers) %>% as.numeric(),
                            passengers = as.numeric(AirPassengers))
ggplot(air_passengers) +
  geom_line(aes(t, passengers)) +
  labs(x = 'Index', y = 'Observation')

ggplot(air_passengers) +
  geom_line(aes(t, c(NA_real_, diff(passengers)))) +
  labs(x = 'Index', y = '1st-order difference')

#' see `?AirPassengers` for more info

# some example AR, MA, and ARMA models
d_ts <- tibble(t = 1:1e3,
               w = rnorm(length(t)),
               ar_1 = if_else(t == 1, w, NA_real_), # AR(1) process
               ma_1 = if_else(t == 1, w, NA_real_), # MA(1) process
               arma_11 = if_else(t == 1, w, NA_real_))

for(i in 2:nrow(d_ts)) {
  d_ts$ar_1[i] <- d_ts$w[i] + d_ts$ar_1[i - 1] * 0.9
  d_ts$ma_1[i] <- d_ts$w[i] + d_ts$w[i - 1] * 0.9
  d_ts$arma_11[i] <- d_ts$w[i] + d_ts$ar_1[i - 1] * 0.9 + d_ts$w[i - 1] * 0.9
  rm(i)
}

#' *white noise*: `y_t` is IID Gaussian: `w_t`
ggplot(d_ts, aes(t, w)) +
  geom_hline(yintercept = 0, lty = 'dashed') +
  geom_line()

layout(t(1:2))
acf(d_ts$w, ci = 0.99) # first lag-0 AR will always be at 1
pacf(d_ts$w, ci = 0.99) # no appreciable partial correlations
arima(d_ts$w, order = c(0, 0, 0)) %>% coef()

#' *AR(1) process*: `y_t` is a weighted average of of `y_{t-i}` and `w_t`
ggplot(d_ts, aes(t, ar_1)) +
  geom_hline(yintercept = 0, lty = 'dashed') +
  geom_line()

acf(d_ts$ar_1, ci = 0.99); pacf(d_ts$ar_1, ci = 0.99)
m_ar <- arima(d_ts$ar_1, order = c(1, 0, 0))
coef(m_ar)
acf(resid(m_ar)); pacf(resid(m_ar)) # residuals don't have appreciable signals

#' *MA(1) process*: errors of `y_t` are a weighted average of `w_t` and `w_{t-1}`
ggplot(d_ts, aes(t, ma_1)) +
  geom_hline(yintercept = 0, lty = 'dashed') +
  geom_line()

acf(d_ts$ma_1, ci = 0.99); pacf(d_ts$ma_1, ci = 0.99)
m_ma <- arima(d_ts$ma_1, order = c(0, 0, 1))
coef(m_ma)
acf(resid(m_ma)); pacf(resid(m_ma))

#' *ARMA(1, 1) process*: `AR(1) + MA(1)`
ggplot(d_ts, aes(t, arma_11)) +
  geom_hline(yintercept = 0, lty = 'dashed') +
  geom_line()

acf(d_ts$arma_11, ci = 0.99); pacf(d_ts$arma_11, ci = 0.99)
m_arma <- arima(d_ts$ma_1, order = c(1, 0, 1))
coef(m_arma)
acf(resid(m_arma)); pacf(resid(m_arma))
layout(1)

## all three models predict that the data will be stationary and centered at the
## intercept with some random variation over time
n_preds <- 100

preds <-
  tibble(model_type = c('AR(1)', 'MA(1)', 'ARMA(1, 1)', 'White noise'),
         model = list(m_ar, m_ma, m_arma, 'White noise'),
         predictions = map(model, function(.m) {
           if(class(.m) == 'character') {
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
  pivot_longer(! t, names_to = 'model_type') %>%
  mutate(model_type = case_when(model_type == 'ar_1' ~ 'AR(1)',
                                model_type == 'ma_1' ~ 'MA(1)',
                                model_type == 'arma_11' ~ 'ARMA(1, 1)',
                                model_type == 'w' ~ 'White noise')) %>%
  ggplot() +
  facet_wrap(~ model_type) +
  geom_hline(aes(yintercept = mean(value)), lty = 'dashed') +
  geom_point(aes(t, value), alpha = 0.5) +
  geom_ribbon(aes(t, ymin = lwr_90, ymax = upr_90), preds, alpha = 0.3,
              lwd = 0.5, , fill = 'darkorange', color = 'darkorange') +
  geom_line(aes(t, pred), preds, color = 'darkorange', lwd = 1) +
  labs(x = 'Timepoint', y = 'Value')

#' see `?AirPassengers` for an example that includes seasonal trends

##' *NOTE:* processes are stationary if `sum(coefs)` < 1 (not intercept)
plot_process <- function(ar = c(0), ma = c(0), return_values = FALSE) {
  out <- rep(NA_real_, nrow(d_ts)) # vector of values to return
  w <- d_ts$w
  
  for(i in 1:nrow(d_ts)) {
    ## create indices to reference the previous observations
    ar_indices <- (i - length(ar)):(i - 1)
    ma_indices <- (i - length(ma)):(i - 1)
    
    ##' negative and `integer(0)` indices cause issues
    ar_indices[ar_indices <= 0] <- NA_integer_
    ma_indices[ma_indices <= 0] <- NA_integer_
    
    out[i] <- sum(w[i], sum(out[ar_indices] * ar), sum(w[ma_indices] * ma),
                  na.rm = TRUE)
  }
  
  layout(matrix(c(1, 1, 2, 3), ncol = 2, byrow = TRUE))
  plot(out, type = 'l')
  acf(out, main = NULL)
  pacf(out, main = NULL)
  layout(1)
  if(return_values) return(out)
}

plot_process(ar = 0, ma = 0) # white noise

# AR processes: stationary if coefficients sum to <= 1
plot_process(ar = c(0.7, 0.2)) # stationary since 0.7 + 0.2 < 1
plot_process(ar = c(0.7, 0.6)) # non-stationary since 0.7 + 0.6 > 1
plot_process(ar = c(0.5, 0.4, 0.3)) # non-stationary since sum > 1
plot_process(ar = c(0.4, 0.3, 0.2, 0.1)) # stationary since sum = 1

# MA processes are always stationary: coefficients only scale the variance
plot_process(ma = c(0.9, 0.1))
plot_process(ma = c(7, 6))
plot_process(ma = c(700, 50))

# ARMA processes: stationary if AR is stationary
plot_process(ar = c(0.7, 0.2), ma = c(0.5, 0.3)) # stationary ARMA(2, 2)

## applying the plots to the air temperature series ----
layout(t(1:2))
acf(d_temp$temp) #' correlation of data pairs at a times `t` and `t+Lag`
pacf(d_temp$temp) #' correl. between pairs, without previous lags' effects
## ACF decays smoothly; high pacf value at lag 1, other values are small
## AR(1) model is a good start

##' `y_t = 0.82 * y_{t-1} + 77.33`
m_ar_temp <- arima(d_temp$temp, order = c(1, 0, 0))
coef(m_ar_temp)
acf(resid(m_ar_temp)); pacf(resid(m_ar_temp)) # residuals are ok

## strengths and limitations of ARMA models:
## - because coefficient estimates depend on pairs of points, the models are
##   robust to missing data, whether random or even missing sections. however,
##   estimates are sensitive to sampling interval and data thinning (see below)
## - AR and MA models can help understand the properties of the data, but they
##   cannot be used to interpolate or estimate the mean trend in a time series
##   other than the intercept term
## - AR and MA models assume the models are stationary over time, so predictions
##   past the range of the data are often not useful (if not in the short term)
## - the models only work in discrete time with fixed sampling intervals

## ARMA models are sensitive to data thinning
d_thin <-
  tibble(thinning = c(1, 10, 20, 30, 40, 50, 60), # thinning intervals
         data = map(thinning, \(.t) d_ts %>% filter(1:n() %% .t == 0) %>%
                      rename(y = ar_1)),
         model = map(data, \(.d) arima(.d$y, c(1, 0, 0))),
         intercept = map_dbl(model, \(.m) coef(.m)['intercept']),
         ar_1 = map_dbl(model, \(.m) coef(.m)['ar1']))

d_thin %>%
  unnest(data) %>%
  ggplot() +
  facet_wrap(~ thinning) +
  geom_line(aes(t, y))

d_thin %>%
  select(thinning, ar_1, intercept) %>%
  pivot_longer(c(ar_1, intercept)) %>%
  ggplot(aes(thinning, value)) +
  facet_wrap(~ name, scales = 'free_y', ncol = 1) +
  geom_point()

# data thinning is particularly problematic with animal movement data ----
#' simulate a movement path from a continuous-time stochastic process
mm <- ctmm(tau = c(180, 5) %#% 'minutes', sigma = 10000, mu = c(0, 0))
track <- simulate(mm, nsim = 1, seed = 160, t = 1:1e3 %#% 'minutes')
plot(track)

# calculate speed using a discrete-time approach: straight-line displacement
data.frame(track) %>%
  mutate(displacement = sqrt(x^2 + y^2),
         time_interval = t - lag(t),
         speed = displacement / time_interval) %>%
  ggplot(aes(t, speed)) +
  geom_line()

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

## estimated speed decreases substantially with sampling interval...
ggplot(d_track, aes(sampling_interval, mean_speed)) +
  geom_line() +
  labs(x = 'Sampling interval (seconds)', y = 'Mean speed (m/s)')

## ... which is not surprising since we are losing data on the complexity of the
## tracks as we thin them. the issues compound if sampling intervals are
## irregular. for more info, see:
## - https://doi.org/10.1186/s40462-019-0177-1
## - https://doi.org/10.1086/675504
## - https://doi.org/10.1101/2025.07.17.665364
d_track %>%
  filter(thinning %in% c(1, 5, 10, 15, 30, 60)) %>%
  mutate(sampling_interval = paste0(sampling_interval / 60, ' (minutes)') %>%
           factor(., levels = unique(.))) %>%
  unnest(data) %>%
  ggplot(aes(x, y)) +
  coord_equal() +
  facet_wrap(~ sampling_interval) +
  geom_path(aes(x, y), d_track$data[[1]], alpha = 0.3) +
  geom_path() +
  geom_point() +
  scale_x_continuous(name = 'x (meters)', breaks = (-3:3) * 100) +
  ylab('y (meters)')

##' **break**

## GLMs and GAMs for ecological modelling ----

##' three main parts to a GLM/GAM:
##' 1. *family* of distributions: distribution of response, Y
##' 2. *linear predictor*: sum of coefficients multiplied by predictor variables
##' 3. *link function*: connects linear predictor with parameter estimates

##' linear models are Gaussian GLMs:
##' 1. family is Gaussian with mean `mu` and variance `sigma^2`
##' 2. linear predictor is `beta_0 + x_1 * beta_1 + ...`
##' 3. link function is the identity function: `I(c) = c`: input = output

##' *fit the model to the data, not the data to the model!*
##' choose a family of distributions and link function based on:
##' 1. the possible values of the response variable
##' 2. the mean-variance relationship
##' 3. any additional considerations about the variance, such as overdispersion

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

##' *choose a link function based on the support of the distribution*
##' unbounded: identity; `I(-Inf, Inf) = (-Inf, Inf)`
##' `Y >= 0` or `Y > 0`: `log(0, Inf) = (-Inf, Inf)`
##' `0 <= Y <= 1`: `logit(0, 1) = log(odds(0, 1)) = log(0,Inf) = (-Inf,Inf)`
##' there are other options, but these are generally sufficient (esp. with GAMs)

##' *note:* link functions introduce two new terms:
##' - *response scale*: the original response values; e.g., (0, Inf), (0, 1)
##' - *link scale*: the transformed response values; generally (-Inf, Inf)

##' *note:* link function is applied to the *mean*, not to the data directly

##' `E(Y) = mu`
##' `g(mu) = eta = b_0 + b_1 * x_1 + b_2 * x_2`

## how can we model data with irregular sampling over time?
##' *NOTE:* many of the `{mvgam}` plots assume discrete-time sampling, so the
##' missing observations should be `NA` rather than missing the full row, as
##' long as none of the predictors have `NA` values.
d_temp_missing <- d_temp %>%
  mutate(time = 1:n()) %>%
  mutate(temp = if_else(month(date) == 6, NA_real_, temp),
         temp = if_else(time %in% sample(time, size = n() / 2),
                        NA_real_, temp)) %>%
  arrange(date)

ggplot(d_temp_missing, aes(date, temp)) +
  geom_point(alpha = 0.75) +
  labs(x = NULL, y = expression(bold(paste('Temperature (\U00B0', 'C)'))))

## focus on rates of change and trends over time rather than changes over steps
## in discrete time

##' fitting a *GLM* with a polynomial term
##' The terms can't be functions of each other, so we need to add columns of the
##' polynomial that are independent of each other (i.e., orthogonal) to avoid
##' complete collinearity and non-identifiability issues when fitting.
d_temp_missing <- d_temp_missing %>%
  bind_cols(.,
            poly(.$doy, degree = 3) %>%
              as.data.frame() %>%
              rename(doy_1 = 1, doy_2 = 2, doy_3 = 3))
d_temp_missing ##' note the `NA`s in the `temp` column

m_temp_poly <- mvgam(temp ~ doy_1 + doy_2 + doy_3,
                     family = gaussian(), ##' *NOTE:* default family is Poisson
                     data = d_temp_missing)

##' since `{mvgam}` fits Bayesian models with `Stan`, we should check that all
##' chains converged properly: check `Rhat`, `n_eff`, and Stan MCMC diagnostics
summary(m_temp_poly)
plot(m_temp_poly, type = 'residuals') #' model diagnostics `plot_mvgam_resids()`

## can add predictions to the data...
d_temp_missing %>%
  bind_cols(predict(m_temp_poly, type = 'response') %>%
              as.data.frame() %>%
              rename(est_poly = Estimate,
                     se_poly = Est.Error,
                     q2.5_poly = Q2.5,
                     q97.5_poly = Q97.5))

## ... but there are also many useful built-in functions
plot(hindcast(m_temp_poly, type = 'response')) # uncertainty in Y
plot(hindcast(m_temp_poly, type = 'link')) # uncertainty in mu on link scale
plot(hindcast(m_temp_poly, type = 'expected')) # uncertainty in mu
plot(m_temp_poly, type = 'smooths') # only works with GAMs (see below)

##' fitting a *GAM* with a smooth term
##' the smooth term is created using the `s()` function
##' greater model complexity requires a bit more sampling and burnin
m_temp_gam <- mvgam(temp ~ s(doy, k = 10, bs = 'cr'),
                    family = gaussian(), data = d_temp_missing,
                    parallel = TRUE, burnin = 1e3, samples = 750,
                    #' increase `adapt_delta` to improve sampling:
                    #' - `max_treedepth`: max n of binary choices when sampling
                    #' - `adapt_delta`: target average proposal acceptance prob. 
                    control = list(max_treedepth = 10, adapt_delta = 0.9))

##' `summary()` looks a bit different from the one for the GLM:
##' - `s(doy)` has `k - 1` coefficients
##' - each coefficient is multiplied by the respective basis function
summary(m_temp_gam)
coef(m_temp_gam$mgcv_model) # model coeffients
coef(m_temp_gam$mgcv_model)[-1] # drop intercept term

#' visualize the cubic basis
draw(basis(s(doy, bs = 'cr'), data = d_temp_missing)) # default cubic basis
draw(basis(s(doy, bs = 'cr'), data = d_temp_missing, coefficients = 1:9,
           constraints = TRUE))
draw(basis(s(doy, bs = 'cr'), data = d_temp_missing, # default basis * coefs
           coefficients = coef(m_temp_gam$mgcv_model)[-1], 
           constraints = TRUE))
draw(basis(m_temp_gam$mgcv_model)) & # fitted cubic basis
  geom_line(aes(doy, .fitted - coef(m_temp_gam$mgcv_model)[1]),
            fitted_values(m_temp_gam$mgcv_model), lwd = 1,
            inherit.aes = FALSE)
plot(hindcast(m_temp_gam, type = 'response')) # uncertainty in Y
plot(hindcast(m_temp_gam, type = 'link')) # uncertainty in mu on link scale
plot(hindcast(m_temp_gam, type = 'expected')) # uncertainty in mu

##' smooths are centered at 0
plot(m_temp_gam, type = 'smooths') #' model terms; == `plot_mvgam_smooth()`

plot(m_temp_gam, type = 'residuals') #' model diagnostics; `plot_mvgam_resids()`

##' **break**

# Temporal random effects and temporal residual correlation structures ----
## fit a GAM with a random effect of week
ggplot(d_temp_missing) + geom_point(aes(week_re, temp), alpha = 0.3)

m_temp_re <- mvgam(temp ~ s(week_re, bs = 're'), family = gaussian(),
                   data = d_temp_missing, parallel = TRUE, silent = 2)
summary(m_temp_re)
draw(m_temp_re$mgcv_model)

plot(hindcast(m_temp_re, type = 'response')) # response scale; uncertainty in Y

##' *Q:* how do we choose the window width? is 2 weeks or 10 days better than 7?

##' smooth terms in GAMs can be thought of as a continuous version of these
##' discrete-time random effects. The random effects in the GAMs are the basis
##' coefficients
plot_grid(
  draw(basis(s(doy, bs = 'cr'), data = d_temp_missing)), # default cubic basis
  draw(basis(m_temp_gam$mgcv_model), residuals = TRUE) + # fitted cubic basis
    geom_line(aes(doy, .fitted - coef(m_temp_gam$mgcv_model)[1]),
              fitted_values(m_temp_gam$mgcv_model), lwd = 1,
              inherit.aes = FALSE),
  ncol = 1)

##' smooth terms are better at dealing with gaps and irregular sampling. they
##' also don't require choosing a window size, but you do need to choose `k`.
ggplot() +
  geom_line(aes(doy, Estimate, color = 's(doy)'),
            bind_cols(d_temp, predict(m_temp_gam, d_temp)),
            lwd = 1, inherit.aes = FALSE) +
  geom_line(aes(doy, Estimate, color = 's(week, bs = \'re\')'),
            bind_cols(d_temp, predict(m_temp_re, d_temp)),
            lwd = 1, inherit.aes = FALSE) +
  geom_point(aes(doy, temp), d_temp_missing, alpha = 0.3, inherit.aes = FALSE) +
  scale_color_highcontrast(name = 'Model') +
  theme(legend.position = 'top')
