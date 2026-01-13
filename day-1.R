# attach necessary packages
source('packages.R')

# Introduction to time series and time series visualization ----
# a simple time series: daily temperature values
data('airquality')
?airquality
head(airquality)

airquality %<>%
  janitor::clean_names() %>% # convert to snake_case
  mutate(date = as_date(paste0('1973-', month, '-', day)),
         doy = yday(date)) %>%
  as_tibble()
airquality

# plot the data
p_0 <-
  ggplot(airquality, aes(date, temp)) +
  geom_point(alpha = 0.75) +
  labs(x = NULL, y = expression(bold(paste('Temperature (\U00B0', 'C)'))))
p_0

# plot the data with a smooth trend
p_0 + geom_smooth(color = 'darkorange', fill = 'darkorange', method = 'gam',
                  formula = y ~ s(x, k = 5))

# plot the data with a wiggly trend
p_0 + geom_smooth(color = 'darkorange', fill = 'darkorange', method = 'gam',
                  formula = y ~ s(x, k = 50), n = 400)

#' *questions:*
#' - which model is better?
#' - what counts as a change in the time series?
#' - what underlying process are we trying to estimate from the data?

#' *HERE*

# Some traditional time series models and their assumptions ----


# GLMs and GAMs for ecological modelling ----


# Temporal random effects and temporal residual correlation structures ----


