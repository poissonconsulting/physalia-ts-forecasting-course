source('packages.R') # attach necessary packages

## Introduction to time series and time series visualization ----
## a simple time series: daily temperature values
data('airquality')
?aq
head(aq) #' *NOTE:* temperature is in Fahrenheit degrees

## clean up the data format
aq %<>%
  janitor::clean_names() %>% # convert to snake_case
  mutate(date = as_date(paste0('1973-', month, '-', day)),
         doy = yday(date)) %>%
  select(temp, month,day, date, doy) %>%
  as_tibble()
aq

## plot the data
p_aq <-
  ggplot(aq, aes(date, temp)) +
  geom_point(alpha = 0.75) +
  labs(x = NULL, y = expression(bold(paste('Temperature (\U00B0', 'C)'))))
p_aq

## plot the data with a smooth model
p_aq + geom_smooth(color = 'darkorange', fill = 'darkorange', method = 'gam',
                   formula = y ~ s(x, k = 5))

## plot the data with a wiggly model
p_aq + geom_smooth(color = 'darkorange', fill = 'darkorange', method = 'gam',
                   formula = y ~ s(x, k = 50), n = 400)

## plot the data with an extremely wiggly model
p_aq + geom_smooth(color = 'darkorange', fill = 'darkorange', method = 'gam',
                   formula = y ~ s(x, k = nrow(aq) - 1), n = 400,
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
acf(aq$temp) #' correlation of data pairs at a times `t` and `t+Lag`
pacf(aq$temp) #' correl. between pairs, without previous lags' effects
## ACF decays smoothly; high pacf value at lag 1, other values are small
## AR(1) model is a good start

##' `y_t = 0.82 * y_{t-1} + 77.33`
m_ar_temp <- arima(aq$temp, order = c(1, 0, 0))
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
## tracks as we thin them. for more info, see:
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

# GLMs and GAMs for ecological modelling ----

# TODO: continue here

p_aq
airquality

m_aq <- mvgam(temp ~ s(doy), data = airquality)

plot(m_aq, type = 'series') #' time series diagnostics; == `plot_mvgam_series()`
plot(m_aq, type = 'residuals') #' model diagnostics; == `plot_mvgam_resids()`
plot(m_aq, type = 'smooths') #' model terms; == `plot_mvgam_smooth()`

plot(hindcast(m_aq, type = 'response'))
plot(hindcast(m_aq, type = 'link'))
plot(hindcast(m_aq, type = 'expected'))


##' **break**

# Temporal random effects and temporal residual correlation structures ----


