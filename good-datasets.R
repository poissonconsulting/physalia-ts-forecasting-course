source('packages.R')

# some useful datasets for examples or testing
plot(Nile)
layout(t(1:2))
acf(Nile)
pacf(Nile)
layout(1)

library('gamair')
data(chicago)
head(chicago)
ggplot(chicago, aes(time, death, color = log(pm10median))) +
  geom_line() +
  scale_color_batlowK()

data(chl)
head(chl)
ggplot(chl) +
  geom_point(aes(lon, lat, color = chl)) +
  scale_color_bamako(reverse = TRUE)

# but mvgam does not support the gev family
data(swer)
head(swer)
ggplot(swer, aes(year, exra, group = location)) +
  geom_line()

library('timeSeriesDataSets')
library('tsibble') #' to read in quarters correctly for `tourism_tbl_ts` below

# this one is excellent
pedestrian <- pedestrian_tbl_ts %>%
  rename_with(tolower) %>%
  mutate(sensor = factor(sensor))
n_distinct(pedestrian$sensor)
range(pedestrian$date)
range(pedestrian$time)

pacf(filter(pedestrian, sensor == sensor[1])$count)
acf(filter(pedestrian, sensor == sensor[1])$count)

pedestrian %>%
  select(count, sensor, date_time) %>%
  pivot_wider(values_from = count, names_from = sensor) %>%
  as_tibble() %>% # convert from tsibble (won't drop date otherwise)
  select(! date_time) %>%
  plot()

pedestrian %>%
  as_tibble() %>% # convert from tsibble (forcibly grouped by date otherwise)
  summarize(count = mean(count), .by = c(date, sensor)) %>%
  ggplot(aes(date, count, group = sensor)) +
  facet_wrap(~ sensor) +
  geom_line()

# this one is also very good
tourism_tbl_ts

acf(filter(tourism_tbl_ts, Region == Region[1])$Trips)
pacf(filter(tourism_tbl_ts, Region == Region[1])$Trips)

ggplot(tourism_tbl_ts, aes(as.Date(Quarter), Trips, group = Region, color = Purpose)) +
  facet_grid(State ~ Purpose) +
  geom_line() +
  scale_color_bright()

