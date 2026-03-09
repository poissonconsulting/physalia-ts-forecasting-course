# Extended practical examples using {mvgam}

# Review and open discussion 

# fitting and comparing models for analyzing a dataset of choice

# HGAM of pigments
SAMPLING_DATE <- lubridate::decimal_date(as.POSIXlt('2014-04-01'))

pigments <-
  bind_rows(
    read.xlsx('https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%201%20April%202014.xlsx') %>%
      mutate(core = 'Core 1',
             interval = c(SAMPLING_DATE, YEAR[-length(YEAR)]) - YEAR,
             weight = interval / mean(interval)),
    read.xlsx('https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%202%20April%202014.xlsx') %>%
      mutate(core = 'Core 2',
             CHLA = CHL_A,
             interval = c(SAMPLING_DATE, YEAR[-length(YEAR)]) - YEAR,
             weight = interval / mean(interval)),
    read.xlsx('https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%203%20April%202014.xlsx') %>%
      mutate(core = 'Core 3',
             interval = c(SAMPLING_DATE, YEAR[-length(YEAR)]) - YEAR,
             weight = interval / mean(interval)),
    read.xlsx('https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%204%20April%202014.xlsx') %>%
      mutate(core = 'Core 4',
             interval = c(SAMPLING_DATE, YEAR[-length(YEAR)]) - YEAR,
             weight = interval / mean(interval))) %>%
  select(core, YEAR, ALLOX, DIATOX, CANTH, PHEO_B, BCAROT, CHL_PHEO,
         interval, weight) %>%
  rename(allo = ALLOX,
         b_car = BCAROT) %>%
  rename_with(tolower)

pigments %>%
  pivot_longer(allo:chl_pheo, names_to = 'pigment', values_to = 'conc') %>%
  ggplot(aes(year, conc)) +
  facet_grid(pigment ~ core, scales = 'free') +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'gam', color = 'black', formula = y ~ s(x, k = 15)) +
  scale_fill_brewer(name = 'Core', aesthetics = c('fill', 'color'))
