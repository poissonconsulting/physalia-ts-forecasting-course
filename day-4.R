source('packages.R') # attach necessary packages

# multiple cores from the same lake: multiple time series for the same pigments
#' core locations: `https://onlinelibrary.wiley.com/cms/asset/a02d5fe1-044e-4dd2-b50c-ff33f328d953/fwb14192-fig-0001-m.jpg`
SAMPLING_DATE <- lubridate::decimal_date(as.POSIXlt('2014-04-01'))

pigments <-
  bind_rows(
    read.xlsx('https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%201%20April%202014.xlsx') %>%
      mutate(core = 'Core 1'),
    read.xlsx('https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%202%20April%202014.xlsx') %>%
      mutate(core = 'Core 2') %>%
      rename(CHLA = CHL_A),
    read.xlsx('https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%203%20April%202014.xlsx') %>%
      mutate(core = 'Core 3'),
    read.xlsx('https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%204%20April%202014.xlsx') %>%
      mutate(core = 'Core 4')) %>%
  select(core, YEAR, DIATOX) %>%
  rename_with(tolower) %>%
  filter(! is.na(year)) %>%
  mutate(year = round(year),
         interval = c(SAMPLING_DATE, year[-length(year)]) - year,
         weight = interval / mean(interval),
         .by = core) %>%
  summarize(diatox = mean(diatox),
            .by = c(core, year)) %>%
  as_tibble() %>%
  right_join(expand_grid(year = seq(min(.$year), max(.$year), by = 1),
                         core = unique(.$core)),
             by = c('year', 'core')) %>%
  mutate(core = factor(core),
         series = core,
         time = year) %>%
  arrange(time, core)

View(pigments)

ggplot(pigments, aes(year, diatox)) +
  facet_wrap(. ~ core) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'gam', color = 'black', formula = y ~ s(x, k = 15)) +
  ylim(c(0, NA))


## Multivariate ecological time series

## Vector autoregressive processes
#' *model the same pigment across multiple cores with VARPs*

## Dynamic factor models

## Multivariate forecast evaluation
