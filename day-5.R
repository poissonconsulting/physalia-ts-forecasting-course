source('packages.R') # attach necessary packages
library('timeSeriesDataSets') # for datasets
library('tsibble') #' to read in quarters correctly for `tourism_tbl_ts` below

# Extended practical examples using {mvgam}
#' more good datasets included in `good-datasets.R`

# hourly pedestrian counts ----
ped <- timeSeriesDataSets::pedestrian_tbl_ts %>%
  rename_with(tolower) %>%
  mutate(date = as_date(date_time),
         year = year(date),
         doy = yday(date),
         dec_time = hour(date_time) + minute(date_time) / 60)
ped
table(ped$sensor) # four sensors at four crosswalks in Melbourne, Australia

pacf(filter(ped, sensor == sensor[1])$count)
acf(filter(ped, sensor == sensor[1])$count)

# relationships across four locations
ped %>%
  select(count, sensor, date_time) %>%
  pivot_wider(values_from = count, names_from = sensor) %>%
  as_tibble() %>% # convert from tsibble (won't drop date otherwise)
  select(! date_time) %>%
  plot()

# time series plots
ped %>%
  as_tibble() %>% # convert from tsibble (forcibly grouped by date otherwise)
  summarize(count = mean(count), .by = c(date, sensor)) %>%
  ggplot(aes(date, count)) +
  facet_wrap(~ sensor) +
  geom_line()

# plot by time of day: different trends on weekends?
ped %>%
  as_tibble() %>% # convert from tsibble (forcibly grouped by date otherwise)
  mutate(h = lubridate::hour(date_time)) %>%
  ggplot(aes(h, count, group = date)) +
  facet_wrap(~ sensor) +
  geom_line(alpha = 0.05)

# daily Australian domestic overnight trips ----
tourism_tbl_ts

layout(1:2)
acf(filter(tourism_tbl_ts, Region == Region[1])$Trips)
pacf(filter(tourism_tbl_ts, Region == Region[1])$Trips)
layout(1)

ggplot(tourism_tbl_ts, aes(as.Date(Quarter), Trips, group = Region, color = Purpose)) +
  facet_grid(State ~ Purpose) +
  geom_line() +
  scale_color_bright()
