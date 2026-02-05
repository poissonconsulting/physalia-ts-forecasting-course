source('packages.R') # attach necessary packages
source('gaussian-process-functions.R') # for plotting GP covariance function

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
##' `{mvgam}` does not drop rows with `NA` response values, so it keeps track of
##' which values are temporally adjacent. for a comparison with `{brms}` see:
##' `https://github.com/nicholasjclark/physalia-forecasting-course/blob/main/day2/tutorial_2_physalia.html`
##' other advantages of using `{mvgam}` over `{brms}` include:
##' - `{mvgam}` allows each time series to have different AR1 parameters
##' - `{mvgam}` can model the correlations among errors of each time series
##' - `{mvgam}` can fit the dynamic processes using a State-Space approach
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
##' `{brms}` would assume that each observation follows the previous row!

## compare to a simple GAM
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
plot_grid(
  plot(forecast(m_gam_missing)) +
    geom_point(aes(time, passengers), air_passengers, color = 'white', size = 1.5) +
    geom_point(aes(time, passengers), air_passengers, color = 'black', size = 1) +
    ylim(c(0, 1e3)),
  plot(forecast(m_gam_ar_missing)) +
    geom_point(aes(time, passengers), air_passengers, color = 'white', size = 1.5) +
    geom_point(aes(time, passengers), air_passengers, color = 'black', size = 1) +
    ylim(c(0, 1e3)),
  ncol = 1)

#' **break**

## smooth correlations over time ----
## continuous auto-regressive (CAR) processes
pigments <- openxlsx::read.xlsx('https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%201%20April%202014.xlsx') %>%
  as_tibble() %>%
  rename_with(stringr::str_to_snake, everything()) %>%
  select(sample, mid_depth_cm, year, diatox) %>%
  arrange(desc(mid_depth_cm)) %>% #' sorting by `year` gives odd trends in plots
  mutate(interval = year - lag(year))

## sediment age decreases nonlinearly with sample depth
ggplot(pigments, aes(year, mid_depth_cm)) +
  geom_point(alpha = 0.75) +
  geom_path() +
  xlab('Year CE') +
  scale_y_reverse('Sample depth (cm)')

## time intervals vary substantially across years
ggplot(pigments, aes(interval, mid_depth_cm)) +
  geom_point(alpha = 0.75) +
  geom_path() +
  xlab('Time interval between samples (years)') +
  scale_y_reverse('Sample depth (cm)')

## plot an example time series (diatoxantin is a pigment produced by diatoms)
## diatoms are glass-like algae: https://en.wikipedia.org/wiki/Diatom
lab_diatox <- expression(bold(Diatoxanthin~concentration~'(nmol'~g^{'-1'}~'C)'))

ggplot(pigments, aes(year, diatox)) +
  geom_line() +
  geom_point(alpha = 0.75) +
  labs(x = 'Year CE', y = lab_diatox)

#' *Q:* how to distinguish the true trend from the error in the data?
## Auger-Méthé et al. (2021; https://doi.org/10.1002/ecm.1470):
## The assumptions that the hidden states are autocorrelated (e.g., that a large
## population in year t will likely lead to a large population in year t + 1),
## and that observations are independent once we account for their dependence on
## the states (Fig. 1a), allow SSMs to separate these two levels of
## stochasticity.
m_diatox_0 <- mvgam(formula = diatox ~ s(year, k = 20),
                    family = Gamma(link = 'log'),
                    data = pigments,
                    chains = 4,
                    burnin = 1000,
                    samples = 500,
                    control = list(max_treedepth = 20, adapt_delta = 0.9),
                    parallel = TRUE,
                    silent = 2)

summary(m_diatox_0)
mcmc_plot(m_diatox_0, type = 'trace', variable = '.', regex = TRUE)
draw(m_diatox_0$mgcv_model, n = 200)

plot(m_diatox_0, type = 'series')    # time series not accounting for the model
plot(m_diatox_0, type = 'residuals') # residuals from the model
plot(m_diatox_0, type = 'forecast')  # estimated mean over time with data points

## add a CAR(1) term to account for continuous-time autocorrelation
## add a column of time for the CAR(1) process
pigments_car <- pigments %>%
  mutate(time = year)

#' `AR(1)` fails because sampling is irregular
if(FALSE) {
  m_diatox_ar <- mvgam(formula = diatox ~ 0,
                       trend_formula = ~ s(year, k = 30),
                       trend_model = AR(),
                       family = Gamma(link = 'log'),
                       data = pigments_car,
                       chains = 4,
                       burnin = 750,
                       samples = 500,
                       control = list(max_treedepth = 20, adapt_delta = 0.9),
                       parallel = TRUE,
                       silent = 2)
}

#' fit a model with a continuous `AR(1)`
m_diatox_car <- mvgam(formula = diatox ~ 0,
                      trend_formula = ~ s(year, k = 30),
                      trend_model = CAR(),
                      family = Gamma(link = 'log'),
                      data = pigments_car,
                      chains = 4,
                      burnin = 750,
                      samples = 500,
                      control = list(max_treedepth = 20, adapt_delta = 0.9),
                      parallel = TRUE,
                      silent = 2)

summary(m_diatox_car)
mcmc_plot(m_diatox_car, type = 'trace', variable = '.', regex = TRUE)

plot(m_diatox_car, type = 'residuals') # residuals from the model
plot_predictions(m_diatox_car, 'year') # smooth term of year
plot(hindcast(m_diatox_car)) # predictions with data points
#' *x axis is wrong*
plot(m_diatox_car, type = 'forecast')  # predictions with data points

## why is the term so smooth? we can get a clue by looking at the posterior for
## the CAR(1) coefficient:
mcmc_plot(m_diatox_car, type = 'intervals', variable = 'ar1[1]')

## the trend is so smooth because the model has attributed the changes to the
## error process rather than to the biological process. this causes the model to
## interpret the pulse in diatom abundance as an autocorrelation process rather
## than part of the true trend.

## refit the model, but allow the term's smoothness to vary across the years
## the adaptive spline allows the wiggliness to vary over the years
?mgcv::smooth.construct.ad.smooth.spec

m_diatox_car_ad <- mvgam(formula = diatox ~ 0,
                         trend_formula = ~ s(year, bs = 'ad', k = 30),
                         trend_model = CAR(),
                         family = Gamma(link = 'log'),
                         data = pigments_car,
                         chains = 4,
                         burnin = 750,
                         samples = 1000,
                         control = list(max_treedepth = 30, adapt_delta = 0.95),
                         parallel = TRUE,
                         silent = 2)

summary(m_diatox_car_ad)
mcmc_plot(m_diatox_car_ad, type = 'trace', variable = '.', regex = TRUE)

plot(m_diatox_car_ad, type = 'residuals') # residuals from the model
plot_predictions(m_diatox_car_ad, 'year') # smooth term of year
plot(hindcast(m_diatox_car_ad))  # predictions with data points
#' *x axis is wrong*
plot(m_diatox_car_ad, type = 'forecast')  # predictions with data points

## CAR(1) coefficient estimate is about the same, but the posterior's much wider
plot_grid(mcmc_plot(m_diatox_car, type = 'intervals', variable = 'ar1[1]') +
            xlim(c(0.4, 1)),
          mcmc_plot(m_diatox_car_ad, type = 'intervals', variable = 'ar1[1]') +
            xlim(c(0.4, 1)),
          ncol = 1)

##' `s(year)` coefficients clearly show how the coefficients affect the basis
##' the `rho` coefficients are the smoothness coefficients
mcmc_plot(m_diatox_car_ad, type = 'intervals', variable = '.', regex = TRUE)

## get methods quickly with citations (you need to add the model terms)
how_to_cite(m_diatox_car_ad)

#' **break**

##' Gaussian Processes
##' for more info, see `katbailey.github.io/post/gaussian-processes-for-dummies`
##' 
##' rather than assuming the y values are independent, leverage the properties
##' of time series data by recognizing the autocorrelation across time. To do
##' this, we need to account for the correlation between pairs of observations.
## closer points are more similar (higher covariance)
## the second half of the plot is generally not useful
pigments %>%
  select(year, diatox) %>%
  mutate(row = 1:n()) %>%
  mutate(data_2 = list(rename_with(., \(x) paste0(x, '_2'), everything()))) %>%
  unnest(data_2) %>%
  filter(row_2 >= row) %>% # drop duplicate pairs (e.g., (2, 1) but not (1, 2))
  mutate(distance = abs(year - year_2)) %>%
  summarise(cov = var(diatox, diatox_2),
            #cov = cov(diatox, diatox_2),
            n = n(),
            .by = distance) %>%
  filter(! is.na(cov)) %>% # drop distances with only one value
  ggplot(aes(distance, cov)) +
  geom_point(size = 2.5) +
  geom_point(aes(color = sqrt(n))) +
  labs(x = 'Distance (years)',
       y = expression(bold(Covariance~ '(nmol'^'2'~g^{'-2'}~'C)'))) +
  scale_color_viridis_c(name = 'n', labels = \(x) x^2)

## closer points are less different (lower variance)
tibble(distance = 1:10,
       dat = pigments %>% select(year, diatox) %>% list()) %>%
  unnest(dat) %>%
  mutate(diatox_2 = map2_dbl(distance, year, \(.d, .y) {
    value <- filter(pigments, year == .y + .d)$diatox
    if(length(value) > 0) return(mean(value)) else return(NA_real_)
  })) %>%
  summarize(var = mean((diatox - diatox_2)^2, na.rm = TRUE),
            n = sum(! is.na(diatox_2)),
            .by = distance) %>%
  ggplot(aes(distance, var)) +
  geom_point(size = 2.5) +
  geom_point(aes(color = sqrt(n))) +
  labs(x = 'Distance (years)',
       y = expression(bold(Variance~ '(nmol'^'2'~g^{'-2'}~'C)'))) +
  scale_color_viridis_c(name = 'n', labels = \(x) x^2)

##' just like the previous GAMs were formed by spline bases multiplied by
##' coefficients, GPs can be interpreted as truly continuous, smooth, functions
##' of random effects
ggplot(pigments, aes(round(year / 50) * 50, diatox)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = 'cs', k = 5)) +
  labs(x = 'Year CE', y = lab_diatox)

ggplot(pigments, aes(round(year / 20) * 20, diatox)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = 'gam') +
  labs(x = 'Year CE', y = lab_diatox)

ggplot(pigments, aes(round(year, -1), diatox)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = 'gam') +
  labs(x = 'Year CE', y = lab_diatox)

ggplot(pigments, aes(round(year / 5) * 5, diatox)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = 'gam') +
  labs(x = 'Year CE', y = lab_diatox)

##' set up the response data as multivariate Gaussian rather than IID. the MVN
##' distribution has a mean matrix and a variance-covariance matrix. If the
##' response values are IID with variance 3.2, the vcov matrix will be
diag(10) * 3.2

##' note that the diagonals are the variances, which are all 3.2, while the
##' covariances are all 0, which implies independence (conditional on the model)
##' 
##' if we allow the covariances to be nonzero, we imply correlation (even after
##' accounting for the model), and we can thus learn the values of the response
##' by leveraging such correlation.
##' 
##' however, the number of cells in the vcov matrix is much larger than the
##' number of observations, so we cannot estimate them all individually.
##' 
##' instead, we need to set a prior over the smoothness of the function by
##' choosing a covariance function (the kernel) as a function of the distances
##' between pairs of observations (e.g, distance in time).
##' 
##' A common kernel function is the squared exponential kernel, also known as
##' the Gaussian kernel or radial basis function (RBF) kernel:
##' `K(x_i, x_j) = alpha^2 exp(- (x_i - x_j)^2 / (2 * phi^2))`,
##' where `alpha^2` is a scalar for the variance, `(x_i - x_j)` is the Euclidean
##' distance between the two `x` values, and `phi^2` determines the correlation
##' between observations with a distance of `x_i - x_j`. 
##' 
##' NOTE: many materials on GPs are from a machine learning perspective, which
##' uses fairly different terminology, so learning about GPs can be confusing.
##' `https://www.youtube.com/watch?v=Y2ZLt4iOrXU` is a good resource, but you
##' may need to watch earlier lectures to understand it fully
m_diatox_gp <- mvgam(formula = diatox ~
                       gp(time, # variable for calculating distances
                          c = 5/4, # scalar for range of predictions
                          k = 30, # n of basis functions for approx GPs
                          gr = FALSE, # grouping not supported by mvgam
                          scale = FALSE), # do not divide distances by max
                     trend_model = CAR(),
                     family = Gamma(link = 'log'),
                     data = pigments_car,
                     chains = 4,
                     burnin = 500,
                     samples = 500,
                     control = list(max_treedepth = 20, adapt_delta = 0.95),
                     parallel = TRUE,
                     silent = 2)

summary(m_diatox_gp)

plot(m_diatox_gp, type = 'residuals') # residuals from the model
plot_predictions(m_diatox_gp, 'year') # smooth term of year
plot(hindcast(m_diatox_gp)) # predictions with data points

## GPs allow users to evaluate the continuous-time correlation as a function of
## the distance between observations. In our model, observations are
## conditionally approximately independent after ~10 years.
## intervals are 60% and 90% CIs
as.data.frame(m_diatox_gp, variable = 'gp_', regex = TRUE) %>%
  plot_kernels(max_time = 20)

##' `rho`:
##' - length scale parameter; similar to SD in horizontal direction
##' - may be scaled so that the maximum euclidean distance between points is 1     
##' `alpha`: marginal variability; similar to variance in vertical direction

## Dynamic coefficient models
