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
