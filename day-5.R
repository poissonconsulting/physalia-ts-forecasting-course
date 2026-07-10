source('packages.R') # attach necessary packages
library('timeSeriesDataSets') # for datasets
library('tsibble') #' to read in quarters correctly for `tourism_tbl_ts` below

#' Extended practical examples using `{mvgam}`
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

pacf(filter(ped, sensor == "Birrarung Marr")$count)
acf(filter(ped, sensor == "Birrarung Marr")$count)

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
  ggplot(aes(h, count, group = date, color = doy)) +
  facet_wrap(~ sensor) +
  geom_line(alpha = 0.2) +
  scale_color_bamO()

# daily Australian domestic overnight trips ----
tourism_tbl <- mutate(tourism_tbl_ts,
                      q_date = as.Date(Quarter),
                      year = year(q_date),
                      q = quarter(q_date))
tourism_tbl

layout(1:2)
acf(filter(tourism_tbl_ts, Region == "Adelaide")$Trips)
pacf(filter(tourism_tbl_ts, Region == "Adelaide")$Trips)
layout(1)

ggplot(tourism_tbl, aes(q_date, Trips, group = Region, color = Purpose)) +
  facet_grid(State ~ Purpose) +
  geom_line() +
  scale_color_bright() +
  xlab("Year") +
  theme(legend.position = "none")

ggplot(tourism_tbl, aes(q, Trips, group = paste(year, Region), color = year)) +
  facet_grid(State ~ Purpose) +
  geom_line() +
  scale_color_lipari()
