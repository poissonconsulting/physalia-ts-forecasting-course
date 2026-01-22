source('packages.R') # attach necessary packages

# example with count data: global number of international air passengers ----
data('AirPassengers')

AirPassengers

air_passengers <-
  tibble(time = 1:length(AirPassengers),
         dec_date = time(AirPassengers) %>% as.numeric(),
         year = floor(dec_date),
         month = round((dec_date - year) * 12) + 1, # to ensure January is 1
         passengers = as.numeric(AirPassengers))
air_passengers

## visualize time as a single line
ggplot(air_passengers, aes(dec_date, passengers)) +
  geom_line() +
  labs(x = 'Year CE', y = 'International airline passengers')

## visualize time as a repeating cycle
ggplot(air_passengers, aes(month, passengers, group = year, color = year)) +
  facet_wrap(~ year) +
  geom_line() +
  scale_x_continuous('Month', expand = c(0, 0)) +
  scale_y_continuous('International airline passengers')

## visualize time as a surface
ggplot(air_passengers, aes(year, month, fill = passengers)) +
  geom_raster() +
  scale_x_continuous('Year CE', expand = c(0, 0)) +
  scale_y_continuous('Month', expand = c(0, 0), breaks = 1:12,
                     labels = month.name) +
  scale_fill_bam(name = 'International airline passengers') +
  theme(legend.position = 'top')

## modeling ----
## split data into training and testing sets
ggplot(air_passengers, aes(dec_date, passengers, lty = year > 1957)) +
  geom_line() +
  labs(x = 'Year CE', y = 'International airline passengers') +
  scale_linetype_manual('Dataset', values = c(1, 3), labels = c('Train', 'Test'))

data_train <- filter(air_passengers, year <= 1957)
data_test <- filter(air_passengers, year > 1957)

##' fit a simple gam with `dec_date`
m_gam <- mvgam(formula = passengers ~ s(dec_date, k = 30),
               family = poisson(),
               data = data_train,
               newdata = data_test,
               chains = 4,
               burnin = 500,
               samples = 500,
               control = list(max_treedepth = 20, adapt_delta = 0.9),
               parallel = TRUE)

summary(m_gam)

## predictions for the test dataset are quite good
plot(hindcast(m_gam))

## predictions for the test dataset are quite poor
plot(m_gam, type = 'forecast')

## basis information only extends to the range of the data, so predictions past
## the last time only continue the trend at the last timestamp
draw(basis(m_gam$mgcv_model))

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
draw(basis(m_gam_bs$mgcv_model, data = air_passengers)) +
  geom_rug(aes(dec_date), data_train, inherit.aes = FALSE)

# predictions for training data aren't as good now
plot(hindcast(m_gam_bs))

## predictions for the test data are better, but still not good
layout(1:2)
plot(m_gam, type = 'forecast')
plot(m_gam_bs, type = 'forecast') ## lower DRPS implies better forecasts
layout(1)

##' more rigid models tend to extrapolate better, so we can decompose the trend
##' into `doy` and `year` trends
m_gam_month <- mvgam(passengers ~
                       s(year, k = 9, bs = 'tp') + # k must be <= unique(year)
                       s(month, k = 10, bs = 'cc'),
                     knots = list(month = c(0.5, 12.5)), # ensures smooth cycles
                     family = poisson(),
                     data = data_train,
                     newdata = data_test,
                     chains = 4,
                     burnin = 500,
                     samples = 500,
                     control = list(max_treedepth = 20, adapt_delta = 0.9),
                     parallel = TRUE)

plot(m_gam_month, type = 'smooths')
plot(hindcast(m_gam_month))

layout(1:3)
plot(m_gam, type = 'forecast')
plot(m_gam_bs, type = 'forecast')
plot(m_gam_month, type = 'forecast') # best model
layout(1)

## assuming the seasonal cycle repeats across years allows us to reduce the
## extrapolation to only be an extrapolation across years but not across months

## TODO: add AR and maybe CAR terms

## Dynamic GLMs and Dynamic GAMs

## Autoregressive dynamic processes

## Gaussian Processes

## Dynamic coefficient models



