source('packages.R') # attach necessary packages

# example with count data: global number of international air passengers ----
data('AirPassengers')

AirPassengers
class(AirPassengers)

air_passengers <-
  tibble(time = 1:length(AirPassengers),
         dec_date = time(AirPassengers) %>% as.numeric(),
         year = floor(dec_date),
         month = round((dec_date - year) * 12) + 1, # to ensure January is 1
         passengers = as.numeric(AirPassengers)) # in thousands
air_passengers

## visualize time as a single line
ggplot(air_passengers, aes(dec_date, passengers)) +
  geom_line() +
  labs(x = 'Year CE', y = 'International airline passengers (thousands)')

## visualize time as a repeating cycle
ggplot(air_passengers, aes(month, passengers, group = year, color = year)) +
  facet_wrap(~ year) +
  geom_line() +
  scale_x_continuous('Month', expand = c(0, 0)) +
  scale_y_continuous('International airline passengers (thousands)')

## visualize time as a surface
ggplot(air_passengers, aes(year, month, fill = passengers)) +
  geom_raster() +
  scale_x_continuous('Year CE', expand = c(0, 0)) +
  scale_y_continuous('Month', expand = c(0, 0), breaks = 1:12,
                     labels = month.name) +
  scale_fill_bam(name = 'International airline passengers (thousands)') +
  theme(legend.position = 'top')

## modeling ----
## split data into training and testing sets
ggplot(air_passengers, aes(dec_date, passengers, lty = year < 1958)) +
  geom_line() +
  geom_vline(xintercept = 1958, lty = 'dashed') +
  labs(x = 'Year CE', y = 'International airline passengers (thousands)') +
  scale_linetype_manual('Dataset', values = c(3, 1), labels = c('Test', 'Train'))

data_train <- filter(air_passengers, year <= 1958)
data_test <- filter(air_passengers, year > 1958)

##' fit a simple gam with `dec_date`
##' *NOTE:* fitting times are a lot slower if passengers is multiplied by 1000
m_gam <- mvgam(formula = passengers ~ s(dec_date, k = 30),
               family = poisson(),
               data = data_train,
               newdata = data_test,
               chains = 4,
               burnin = 500,
               samples = 500,
               control = list(max_treedepth = 20, adapt_delta = 0.9),
               parallel = TRUE)

# diagnostics look ok
mcmc_plot(m_gam, type = 'trace', variable = '.', regex = TRUE)
summary(m_gam)

## predictions for the test dataset are quite good
plot(hindcast(m_gam))

## predictions for the test dataset are quite poor
plot(m_gam, type = 'forecast')

## basis information only extends to the range of the data, so predictions past
## the last time only continue the trend at the last timestamp
draw(basis(m_gam$mgcv_model))

#' reduce model complexity to improve prediction accuracy
m_gam_smooth <- mvgam(formula = passengers ~ s(dec_date, k = 10),
                      family = poisson(),
                      data = data_train,
                      newdata = data_test,
                      chains = 4,
                      burnin = 500,
                      samples = 500,
                      control = list(max_treedepth = 20, adapt_delta = 0.9),
                      parallel = TRUE, silent = 2)

# predictions are better, but the model is missing the seasonality
plot(hindcast(m_gam_smooth))
plot(m_gam_smooth, type = 'forecast')

##' we can improve the predictions somewhat by extending the basis, but this
##' still depends on the model and the complexity of the trends...
##' for more info, see `https://fromthebottomoftheheap.net/2020/06/03/extrapolating-with-gams/`
##' fit a GAM with a cubic B spline whose curvature and slope are penalized
m_gam_bs <- mvgam(formula = passengers ~ s(dec_date, k = 30, bs = 'bs',
                                           m = c(3, 2, 1)),
                  knots = list(dec_date = c(min(data_train$dec_date),
                                            max(data_test$dec_date))),
                  family = poisson(),
                  data = data_train,
                  newdata = data_test,
                  chains = 4,
                  burnin = 500,
                  samples = 500,
                  control = list(max_treedepth = 20, adapt_delta = 0.9),
                  parallel = TRUE)

##' warnings indicate there is no data for a portion of the spline (as expected)
draw(basis(m_gam_bs$mgcv_model, data = air_passengers)) +  # model bases
  geom_rug(aes(dec_date), data_train, inherit.aes = FALSE) # rug plot of data

# predictions for training data aren't as good now
plot(hindcast(m_gam_bs))

## predictions for the test data are better, but still not good
## the model doesn't know there are seasonal cycles. it only understands that
## the values go up and down, but not why they do.
layout(1:2)
plot(m_gam, type = 'forecast')
plot(m_gam_bs, type = 'forecast') ## lower DRPS implies better forecasts
layout(1)

##' more rigid models tend to extrapolate better, so we can decompose the trend
##' into `doy` and `year` trends to allow the model to learn the seasonal cycles
m_gam_month <- mvgam(passengers ~
                       s(year, k = 9, bs = 'bs', # k must be <= unique(year)
                         m = c(3, 2, 1)) +
                       s(month, k = 10, bs = 'cc'),
                     knots = list(year = c(range(air_passengers$year)),
                                  month = c(0.5, 12.5)), # ensures smooth cycles
                     family = poisson(),
                     data = data_train,
                     newdata = data_test,
                     chains = 4,
                     burnin = 500,
                     samples = 500,
                     control = list(max_treedepth = 20, adapt_delta = 0.9),
                     parallel = TRUE, silent = 2)

plot(hindcast(m_gam_month)) # predicts similar oscillations over the years
plot(m_gam_month, type = 'smooths') #' trends decomposed into `year` and `month`

#' value, slope, and curvature match at `month = 12.5 = 0.5`
draw(m_gam_month, select = 's(month)',
     data = tibble(month = seq(0, 24, by = 0.001), year = 0))

layout(matrix(1:4, ncol = 2))
plot(m_gam, type = 'forecast')
plot(m_gam_smooth, type = 'forecast')
plot(m_gam_bs, type = 'forecast')
plot(m_gam_month, type = 'forecast') # best model, but over-predicts a bit
layout(1)

## assuming the seasonal cycle repeats across years allows us to reduce the
## extrapolation to only be an extrapolation across years but not across months

# autocorrelation present at lags 1 and 12
plot(m_gam_month)

## Auger-Méthé et al. (2021; https://doi.org/10.1002/ecm.1470):
## The assumptions that the hidden states are autocorrelated (e.g., that a large
## population in year t will likely lead to a large population in year t + 1),
## and that observations are independent once we account for their dependence on
## the states (Fig. 1a), allow SSMs to separate these two levels of
## stochasticity.

# fit a GAM with an AR(1) process
m_gam_ar <- mvgam(formula = passengers ~ 0, # assuming no error in response
                  trend_formula = ~
                    s(year, k = 9, bs = 'tp') +
                    s(month, k = 10, bs = 'cc'),
                  trend_model = AR(p = 1),
                  knots = list(month = c(0.5, 12.5)),
                  family = poisson(),
                  data = data_train,
                  newdata = data_test,
                  chains = 4,
                  burnin = 750,
                  samples = 500,
                  control = list(max_treedepth = 20, adapt_delta = 0.95),
                  parallel = TRUE, silent = 2)

##' the `m_gam_ar` model is:
##' `passengers ∼ Poisson(λ_t)`
##' `log(λ_t) = b_0 + s(year) + s(month) + z_t`
##' `z_t ∼ Normal(z_{t−1} * a, σ)`
##' `σ ∼ Exponential(2)`
##' where `a` is the coefficient of the `AR(1)` process

summary(m_gam_ar) # diagnostics are ok

# there's still some autocorrelation at lag 12, but coef is quite uncertain
# intervals are 50%, 60%, and 90%, estimated from samples
plot(m_gam_ar)

#' *NOTE:* we need to drop rows with `NA` predictors because the model can't use
#' them. specifying the AR process in the `trends_model` argument avoids this
#' issue. This is why we shouldn't just add the lagged values as predictors in a
#' GAM.
data_train_12 <- mutate(data_train, lag_12_passengers = lag(passengers, 12)) %>%
  filter(! is.na(lag_12_passengers))
data_test_12 <- air_passengers %>%
  mutate(lag_12_passengers = lag(passengers, 12)) %>%
  filter(year > 1958, ! is.na(lag_12_passengers))

# add a term to crudely account for a 12-month autocorrelation
m_gam_ar_12 <- mvgam(formula = passengers ~ 0,
                     trend_formula = ~
                       lag_12_passengers +
                       s(year, k = 9, bs = 'tp') +
                       s(month, k = 10, bs = 'cc'),
                     trend_model = AR(p = 1),
                     knots = list(month = c(0.5, 12.5)),
                     family = poisson(),
                     data = data_train_12,
                     newdata = data_test_12,
                     chains = 4,
                     burnin = 750,
                     samples = 500,
                     control = list(max_treedepth = 20, adapt_delta = 0.95),
                     parallel = TRUE)

coef(m_gam_ar_12$trend_mgcv_model)['lag_12_passengers']
plot(m_gam_ar) # see values at lag-12 for ACF and pACF
plot(m_gam_ar_12) # values at lag-12 are smaller, but lag-1 values are larger

##' how does `{mvgam}` handle many missing data?
data_train_missing <- data_train %>%
  mutate(passengers = if_else(1:n() %in% sample(1:n(), n() * 0.9), NA_real_,
                              passengers))

m_gam_ar_missing <-
  mvgam(formula = passengers ~ 0,
        trend_formula = ~
          s(year, k = 9, bs = 'tp') +
          s(month, k = 10, bs = 'cc'),
        trend_model = AR(p = 1),
        knots = list(month = c(0.5, 12.5)),
        family = poisson(),
        data = data_train_missing,
        newdata = data_test,
        chains = 4,
        burnin = 750,
        samples = 500,
        control = list(max_treedepth = 20, adapt_delta = 0.95),
        parallel = TRUE, silent = 2)

plot(m_gam_ar_missing, type = 'forecast')

## compare to a simple GAM
## GAM with AR process has 
m_gam_missing <-
  mvgam(formula = passengers ~
          s(year, k = 9, bs = 'tp') +
          s(month, k = 10, bs = 'cc'),
        knots = list(month = c(0.5, 12.5)),
        family = poisson(),
        data = data_train_missing,
        newdata = data_test,
        chains = 4,
        burnin = 750,
        samples = 500,
        control = list(max_treedepth = 20, adapt_delta = 0.95),
        parallel = TRUE, silent = 2)

## AR GAM has much more uncertainty
plot(forecast(m_gam_missing)) +
  geom_point(aes(time, passengers), air_passengers, color = 'white', size = 1.5) +
  geom_point(aes(time, passengers), air_passengers, color = 'black', size = 1) +
  ylim(c(0, 1e3)) +
plot(forecast(m_gam_ar_missing)) +
  geom_point(aes(time, passengers), air_passengers, color = 'white', size = 1.5) +
  geom_point(aes(time, passengers), air_passengers, color = 'black', size = 1) +
  ylim(c(0, 1e3))
  
#' *HERE*

## Gaussian Processes

## Dynamic coefficient models

## to get methods quickly with appropriate citations (need to add model terms)
how_to_cite(m_gam_ar)
