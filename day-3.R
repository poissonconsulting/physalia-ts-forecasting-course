source("packages.R") # attach necessary packages
source("gaussian-process-functions.R")

#' recap:
#' - ARIMA models:
#'   - AR if *data* are correlated through time
#'   - I if data need to be detrended by taking differences
#'   - MA if *errors* are correlated through time
#' - CAR models are a continuous-time version of AR models
#' - GAM's smooth terms are a continuous version of random effects
#' - to estimate change, data should have 3+ observations per period of interest
#' - state space models separate the observation process from the latent process

#' *State space models*:
#' - process model:                    `mu_proc = b0 + b1 * x1 + ...`
#' - process output (states):          `Y_proc ~ MVN(mu_proc, s_proc)`
#' - process observations at time `t`: `O_t ~ MVN(Y_proc, s_obs)`
#' 
#' but estimating the process model requires us to work backwards:
#' - model the space of possible states (i.e., true outcomes, true responses)
#' - process observations at time `t`: `O_t ~ MVN(Y_proc, s_obs)`
#' - process output (states):          `Y_proc ~ MVN(mu_proc, s_proc)`
#' - process model:                    `mu_proc = b0 + b1 * x1 + ...`
#' uncertainty needs to be propagated accordingly across each step

#' today's topics:
#' - finishing CAR(1) processes
#' - Gaussian processses
#' - dynamic processes
#' - forecasting from dynamic models
#' - interpreting the different types of predictions
#' - comparing models and assessing them with forecasts:
#'   - Point-based forecast evaluation
#'   - Probabilistic forecast evaluation
#'   - Bayesian posterior predictive checks

pigments <- read.xlsx("https://github.com/simpson-lab/wpg-mb-lakes/raw/refs/heads/main/data/mb/Manitoba%20pigs%20isotope%20Core%201%20April%202014.xlsx") %>%
  as_tibble() %>%
  rename_with(stringr::str_to_snake, everything()) %>%
  select(mid_depth_cm, year, diatox, percentn) %>%
  rename(percent_n = percentn) %>%
  # collapse two duplicate years into the same row
  summarize(mid_depth_cm = mean(mid_depth_cm, na.rm = TRUE),
            year = mean(year, na.rm = TRUE),
            diatox = mean(diatox, na.rm = TRUE),
            percent_n = mean(percent_n),
            .by = year) %>%
  arrange(year) %>%
  mutate(interval = year - lag(year),
         series = factor("core 1")) %>%
  #' add a column of time for the `CAR(1)` process
  #' can use `CAR(1)` if times are not integers (equivalent otherwise)
  #' requires a consecutive series of `time` values
  right_join(tibble(year = seq(min(.$year), max(.$year), by = 1)),
             by = "year") %>%
  mutate(time = year,
         series = factor("core 1")) %>%
  arrange(time)

# sediment age decreases nonlinearly with sample depth
ggplot(pigments, aes(year, mid_depth_cm)) +
  geom_point(alpha = 0.75) +
  geom_path() +
  xlab("Year CE") +
  scale_y_reverse("Sample depth (cm)") +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, lwd = 2,
           arrow = arrow(length = unit(0.5, "cm"), ends = "last",
                         type = "closed"), color = "darkorange") +
  # Prevents clipping so left side of arrow is visible
  coord_cartesian(clip = "off") +
  theme(axis.line.y = element_line(color = "darkorange"))

# time intervals vary substantially across years
ggplot(pigments, aes(interval, mid_depth_cm)) +
  geom_path() +
  geom_point(alpha = 0.75) +
  xlab("Time interval between samples (years)") +
  scale_y_reverse("Sample depth (cm)") +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf,
           arrow = arrow(length = unit(0.5, "cm"), ends = "last",
                         type = "closed")) +
  coord_cartesian(clip = "off")

# plot an example time series (diatoxantin is a pigment produced by diatoms)
# diatoms are glass-like algae: https://en.wikipedia.org/wiki/Diatom
lab_diatox <- expression(bold(Diatoxanthin~concentration~"(nmol"~g^{"-1"}~"C)"))

ggplot(pigments, aes(year, diatox)) +
  geom_path() +
  geom_point(alpha = 0.75) +
  labs(x = "Year CE", y = lab_diatox)

#' fit a basic GAM
#' using observation formula because we can't separate it from the latent trend
#' technically, we could add `s(year)` to both formulas, since the values
#' changed over the years (trend process) and the pigment possibly decayed over
#' time (observation process), but we can't disentangle the two with only one
#' time series
m_diatox_0 <- mvgam(formula = diatox ~ s(year, k = 20),
                    family = Gamma(link = "log"),
                    data = pigments,
                    chains = 4,
                    burnin = 500,
                    samples = 500,
                    parallel = TRUE,
                    silent = 2)

#' add a `CAR(1)` term to account for continuous-time autocorrelation
#' the model is:
#' `diatox ∼ Gamma(mu_t, theta)`        # observation
#' `log(mu_t) = s(year) + l_t + z_t`    # trend of mean obs on log scale
#' `l_t = 0`                            # latent process trend varies w time
#' `z_t ∼ Normal(z_{t−dt} * a^{dt}, σ)` # latent stochastic component
#' where `0 < a < 1` is the coef of the `CAR(1)` process for time difference `dt`
#' 
#' note: when `dt = 1`, the model becomes a simple `AR(1)` process:
#' `E(z_t) = z_{t−dt} * a^{dt} = z_{t−1} * a^1 = z_{t−1} * a`
#' 
#' note: when `dt = 0`, `z_t` becomes equal to itself:
#' `E(z_t) = z_{t−dt} * a^{dt} = z_{t-0} * a^0 = z_{t} * 1`
#' 
#' note: as `dt` becomes large, `z_t` becomes independent from `z_{t-dt}`
#' `E(z_t) = z_{t−dt} * a^{Inf} = z_{t-dt} * 0 = z_{t−Inf} * 0`
m_diatox_car <- mvgam(formula = diatox ~ s(year, k = 6),
                      trend_formula = ~ 0,
                      trend_model = CAR(1),
                      noncentred = TRUE,
                      family = Gamma(link = "log"),
                      data = pigments,
                      chains = 4,
                      burnin = 500,
                      samples = 500,
                      control = list(adapt_delta = 0.95),
                      parallel = TRUE)

plot(m_diatox_car, type = "residuals") # diagnostics look great
summary(m_diatox_car) # summary looks good
mcmc_plot(m_diatox_car, type = "trace", variable = ".", regex = TRUE) # good

plot_predictions(m_diatox_car, "year") #' `s(year)` is very smooth...
plot(hindcast(m_diatox_car)) # predictions with data points are not smooth...

# posterior predictive checks: very high uncertainty
pp_check(m_diatox_car, "dens_overlay", ndraws = 100)
pp_check(m_diatox_car, "ecdf_overlay", ndraws = 100)
pp_check(m_diatox_car, "intervals") # misses the spike entirely
pp_check(m_diatox_car, "ribbon")
pp_check(m_diatox_car, "error_scatter_avg") # error is proportional to y
pp_check(m_diatox_car, "scatter_avg") # spike is clearly visible
plot_mvgam_trend(m_diatox_car) # CAR(1) process
pp_check(m_diatox_car, "resid_ribbon") # residuals are uncorrelated after CAR(1)

# why is the year term so smooth and uncertain?
# we can get a clue by looking at the posterior for the CAR(1) coefficient:
#' the model attributes the spike to the `CAR(1)` process instead of `s(year)`
plot_grid(mcmc_plot(m_diatox_car, type = "intervals", variable = "ar1[1]"),
          plot_predictions(m_diatox_car, "year")) # smooth term of year

# the trend is so smooth because the model has attributed the changes to the
# error process rather than to the biological process. this causes the model to
# interpret the pulse in diatom abundance as an autocorrelation process rather
# than part of the true trend.

# refit the model, but allow the term's smoothness to vary across the years
# the adaptive spline allows the wiggliness to vary over the years
?mgcv::smooth.construct.ad.smooth.spec

m_diatox_car_ad <- mvgam(formula = diatox ~ s(year, bs = "ad", k = 15),
                         trend_model = CAR(),
                         noncentred = TRUE,
                         family = Gamma(link = "log"),
                         data = pigments,
                         chains = 4,
                         burnin = 500,
                         samples = 1000,
                         control = list(adapt_delta = 0.95),
                         parallel = TRUE,
                         silent = 2)

summary(m_diatox_car_ad)
mcmc_plot(m_diatox_car_ad, type = "trace", variable = ".", regex = TRUE)

plot(m_diatox_car_ad, type = "residuals") # residuals from the model
plot_predictions(m_diatox_car_ad, "year") # smooth term of year
plot(hindcast(m_diatox_car_ad))  # predictions with data points

# adaptive smooth allows to attribute more change to the smooth and less to CAR
plot_grid(mcmc_plot(m_diatox_car, type = "intervals", variable = "ar1[1]") +
            xlim(c(0, 1)) +
            ggtitle("Thin plate regression spline"),
          mcmc_plot(m_diatox_car_ad, type = "intervals", variable = "ar1[1]") +
            xlim(c(0, 1)) +
            ggtitle("Adaptive spline"),
          ncol = 1)

#' `s(year)` coefficients clearly show how the coefficients affect the basis
#' the `rho` coefficients are the smoothness coefficients
mcmc_plot(m_diatox_car, type = "intervals", variable = ".", regex = TRUE)
mcmc_plot(m_diatox_car_ad, type = "intervals", variable = ".", regex = TRUE)

# get methods quickly with citations (you need to add the model terms)
how_to_cite(m_diatox_car_ad)

# Gaussian Processes ----
#' for more info, see katbailey.github.io/post/gaussian-processes-for-dummies
#' applications of state-space models:
#' Kalman filter and Apollo missions: https://doi.org/10.1109/MCS.2010.936465
#' 
#' rather than assuming the y values are independent, leverage the properties
#' of time series data by recognizing the autocorrelation across time. To do
#' this, we need to account for the correlation between pairs of observations.
# closer points are more similar (higher covariance)
# the second half of the plot is generally not useful
pigments %>%
  select(year, diatox) %>%
  mutate(row = 1:n()) %>%
  mutate(data_2 = list(rename_with(., \(x) paste0(x, "_2"), everything()))) %>%
  unnest(data_2) %>%
  filter(row_2 >= row) %>% # drop duplicate pairs (e.g., (2, 1) but not (1, 2))
  mutate(distance = abs(year - year_2)) %>%
  summarise(cov = var(diatox, diatox_2, na.rm = TRUE),
            #cov = cov(diatox, diatox_2),
            n = n(),
            .by = distance) %>%
  filter(! is.na(cov)) %>% # drop distances with only one value
  ggplot(aes(distance, cov)) +
  geom_point(size = 2.5) +
  geom_point(aes(color = sqrt(n))) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "ad")) +
  labs(x = "Distance (years)",
       y = expression(bold(Covariance~ "(nmol"^"2"~g^{"-2"}~"C)"))) +
  scale_color_viridis_c(name = "n", labels = \(x) x^2)

# closer points are less different (lower variance)
tibble(distance = 1:75,
       dat = pigments %>% select(year, diatox) %>% list()) %>%
  unnest(dat) %>%
  mutate(diatox_2 = map2_dbl(distance, year, \(.d, .y) {
    value <- filter(pigments, year == .y + .d)$diatox
    if(length(value) > 0) {
      return(mean(value, na.rm = TRUE))
    } else {
      return(NA_real_)
    }
  })) %>%
  summarize(var = mean((diatox - diatox_2)^2, na.rm = TRUE),
            n = sum(! (is.na(diatox_2) | is.na(diatox))),
            .by = distance) %>%
  ggplot(aes(distance, var)) +
  geom_point(size = 2.5) +
  geom_point(aes(color = sqrt(n))) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5)) +
  labs(x = "Distance (years)",
       y = expression(bold(Variance~ "(nmol"^"2"~g^{"-2"}~"C)"))) +
  scale_color_viridis_c(name = "n", labels = \(x) x^2, limits = c(0, NA))

#' just like the previous GAMs were formed by spline bases multiplied by
#' coefficients, GPs can be interpreted as truly continuous, smooth, functions
#' of random effects
ggplot(pigments, aes(round(year / 50) * 50, diatox)) +
  geom_point(alpha = 0.75) +
  stat_summary(fun = "mean", geom = "point", col = "darkorange", size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5)) +
  labs(x = "Year CE", y = lab_diatox)

ggplot(pigments, aes(round(year / 20) * 20, diatox)) +
  geom_point(alpha = 0.75) +
  stat_summary(fun = "mean", geom = "point", col = "darkorange", size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 10)) +
  labs(x = "Year CE", y = lab_diatox)

ggplot(pigments, aes(round(year, -1), diatox)) +
  geom_point(alpha = 0.75) +
  stat_summary(fun = "mean", geom = "point", col = "darkorange", size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 20)) +
  labs(x = "Year CE", y = lab_diatox)

ggplot(pigments, aes(round(year / 5) * 5, diatox)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 20)) +
  labs(x = "Year CE", y = lab_diatox)

#' set up the response data as multivariate Gaussian rather than IID. the MVN
#' distribution has a mean matrix and a variance-covariance matrix. If the
#' response values are IID with variance 3.2, the vcov matrix will be
diag(10) # I_{10} identity matrix
3.2 * diag(10)

#' note that the diagonals are the variances, which are all 3.2, while the
#' covariances are all 0, which implies independence (conditional on the model)
#' 
#' if we allow the covariances to be nonzero, we imply correlation (even after
#' accounting for the model), and we can thus learn the values of the response
#' by leveraging such correlation.
#' 
#' however, the number of cells in the vcov matrix is much larger than the
#' number of observations, so we cannot estimate them all individually.
#' 
#' instead, we need to set a prior over the smoothness of the function by
#' choosing a covariance function (the kernel) as a function of the distances
#' between pairs of observations (e.g, distance in time).
#' 
#' A common kernel function is the squared exponential kernel, also known as
#' the Gaussian kernel or radial basis function (RBF) kernel:
#' `K(x_i, x_j) = alpha^2 exp(- (x_i - x_j)^2 / (2 * phi^2))`,
#' where `alpha^2` is a scalar for the variance, `(x_i - x_j)` is the Euclidean
#' distance between the two `x` values, and `phi^2` determines the correlation
#' between observations with a distance of `x_i - x_j`. 
#' 
#' NOTE: many materials on GPs are from a machine learning perspective, which
#' uses fairly different terminology, so learning about GPs can be confusing.
#' https://youtu.be/MtXg7fxQgeA?si=FfPjVZuJHm94XktP is a good resource, but you
#' may need to watch earlier lectures to understand it fully
m_diatox_gp <- mvgam(formula = diatox ~
                       gp(year, # variable for calculating distances
                          c = 5/4, # scalar for range of predictions
                          k = 30, # n of basis functions for approx GPs
                          gr = FALSE, # grouping not supported by mvgam
                          scale = FALSE), # do not divide distances by max
                     trend_model = CAR(),
                     family = Gamma(link = "log"),
                     data = pigments,
                     chains = 4,
                     burnin = 500,
                     samples = 500,
                     control = list(adapt_delta = 0.95),
                     parallel = TRUE,
                     silent = 2)

summary(m_diatox_gp)

plot(m_diatox_gp, type = "residuals") # residuals from the model
plot_predictions(m_diatox_gp, "year") # smooth term of year
plot(hindcast(m_diatox_gp)) # predictions with data points

# GPs allow users to evaluate the continuous-time correlation as a function of
# the distance between observations. In our model, observations are
# conditionally approximately independent after ~40 years.
# intervals are 60% and 90% CIs
as.data.frame(m_diatox_gp, variable = "gp_", regex = TRUE) %>%
  plot_kernels(max_time = 50)

#' `rho`:
#' - length scale parameter; similar to SD in horizontal direction
#' - may be scaled so that the maximum euclidean distance between points is 1     
#' `alpha`: marginal variability; similar to variance in vertical direction

#' *break*

# Dynamic coefficient models ----
#' the GP term above can be seen as the change in the intercept term over time,
#' i.e., an interaction between time and the intercept. We can also use GPs to
#' create interactions of a slope over time, which gives us dynamic coefficient
#' models.
#' add a time-varying effect of percent nitrogen in the (dry) soil to estimate
#' the effects of nutrient input and eutrophication on diatom abundance
#' this model will only work if `time == year`
ggplot(pigments) +
  geom_point(aes(percent_n, diatox))

m_diatox_pn <- mvgam(formula = diatox ~
                       dynamic(percent_n, # variable varying over time
                               k = 20, # n of basis functions for the GP
                               rho = 40, # need to specify manually
                               scale = FALSE), # do not divide distances by max
                     trend_model = CAR(),
                     noncentred = TRUE, # helps avoid conflations between terms
                     family = Gamma(link = "log"),
                     #' fails if there are missing `percent_n` values
                     data = filter(pigments, ! is.na(percent_n)),
                     chains = 4,
                     burnin = 500,
                     samples = 500,
                     parallel = TRUE,
                     silent = 2)

summary(m_diatox_pn)

plot_predictions(m_diatox_pn, "percent_n", type = "expected")

expand_grid(percent_n = c(0.4, 0.6, 0.8),
            time = gratia:::seq_min_max(pigments$year, n = 400),
            series = unique(pigments$series)) %>%
  bind_cols(predict(m_diatox_pn, newdata = ., type = "expected")) %>%
  ggplot(aes(time, Estimate, group = percent_n)) +
  coord_cartesian(ylim = c(0, 200)) +
  geom_ribbon(aes(time, ymin = Q2.5, ymax = Q97.5, fill = percent_n),
              alpha = 0.2) +
  geom_line(lwd = 2) +
  geom_line(aes(color = percent_n), lwd = 1) +
  geom_point(aes(year, diatox), pigments, size = 2) +
  geom_point(aes(year, diatox, color = percent_n), pigments, size = 1) +
  labs(x = NULL, y = lab_diatox) +
  scale_fill_acton(name = "% N (dry weight)", breaks = c(0.4, 0.6, 0.8, 1),
                   reverse = TRUE) +
  scale_color_acton(name = "% N (dry weight)", breaks = c(0.4, 0.6, 0.8, 1),
                    reverse = TRUE) +
  theme(legend.position = "top")

# diatox vs % N
expand_grid(
  percent_n = gratia:::seq_min_max(pigments$percent_n, n = 100),
  year = gratia:::seq_min_max(pigments$year, n = 5),
  series = unique(pigments$series)) %>%
  mutate(time = year) %>%
  bind_cols(., predict(m_diatox_pn, newdata = ., type = "expected")) %>%
  ggplot(aes(percent_n, Estimate, group = year)) +
  coord_cartesian(ylim = c(0, 200)) +
  geom_ribbon(aes(percent_n, ymin = Q2.5, ymax = Q97.5, fill = year),
              alpha = 0.2) +
  geom_line(lwd = 2) +
  geom_line(aes(color = year), lwd = 1) +
  geom_point(aes(percent_n, diatox), pigments, size = 2) +
  geom_point(aes(percent_n, diatox, color = year), pigments, size = 1) +
  labs(x = "% N (dry weight)", y = lab_diatox) +
  scale_color_iridescent(name = "Year") +
  scale_fill_iridescent(name = "Year")

# surface plot
expand_grid(
  percent_n = gratia:::seq_min_max(pigments$percent_n, n = 100),
  year = gratia:::seq_min_max(pigments$year, n = 100),
  series = unique(pigments$series)) %>%
  mutate(time = year) %>%
  bind_cols(., predict(m_diatox_pn, newdata = ., type = "expected")) %>%
  ggplot(aes(year, percent_n, fill = Estimate)) +
  geom_raster() +
  scale_x_continuous("Year CE", expand = c(0, 0)) +
  scale_y_continuous("% N (dry weight)", expand = c(0, 0)) +
  scale_fill_bamako(name = lab_diatox, limits = c(0, 200), na.value = "white") +
  theme(legend.position = "top")

#' allow the effect of `percent_n` to vary smoothly
m_diatox_pn_ti <- mvgam(formula = diatox ~
                          gp(year, c = 5/4, k = 30, gr = FALSE, scale = FALSE) +
                          s(percent_n, k = 5, bs = "tp") +
                          ti(year, percent_n, k = c(10, 5), bs = c("gp", "tp")),
                        trend_model = CAR(),
                        noncentred = TRUE, # helps avoid conflations between terms
                        family = Gamma(link = "log"),
                        data = filter(pigments, ! is.na(percent_n)),
                        chains = 4,
                        burnin = 500,
                        samples = 500,
                        parallel = TRUE,
                        silent = 2)

# diatox vs % N
expand_grid(
  percent_n = gratia:::seq_min_max(pigments$percent_n, n = 100),
  year = gratia:::seq_min_max(pigments$year, n = 5),
  series = unique(pigments$series)) %>%
  mutate(time = year) %>%
  bind_cols(., predict(m_diatox_pn_ti, newdata = ., type = "expected")) %>%
  ggplot(aes(percent_n, Estimate, group = year)) +
  coord_cartesian(ylim = c(0, 200)) +
  geom_ribbon(aes(percent_n, ymin = Q2.5, ymax = Q97.5, fill = year,
                  color = year), alpha = 0.2) +
  geom_line(lwd = 2) +
  geom_line(aes(color = year), lwd = 1) +
  geom_point(aes(percent_n, diatox), pigments, size = 2.5) +
  geom_point(aes(percent_n, diatox, color = year), pigments, size = 1) +
  labs(x = "% N (dry weight)", y = lab_diatox) +
  scale_color_iridescent(name = "Year") +
  scale_fill_iridescent(name = "Year")

# diatox vs year
expand_grid(
  percent_n = gratia:::seq_min_max(pigments$percent_n, n = 3),
  year = gratia:::seq_min_max(pigments$year, n = 400),
  series = unique(pigments$series)) %>%
  mutate(time = year) %>%
  bind_cols(., predict(m_diatox_pn_ti, newdata = ., type = "expected")) %>%
  ggplot(aes(year, Estimate, group = percent_n)) +
  coord_cartesian(ylim = c(0, 200)) +
  geom_ribbon(aes(year, ymin = Q2.5, ymax = Q97.5, fill = percent_n,
                  color = percent_n),
              alpha = 0.2) +
  geom_line(lwd = 2) +
  geom_line(aes(color = percent_n), lwd = 1) +
  geom_point(aes(year, diatox), pigments, size = 2) +
  geom_point(aes(year, diatox, color = percent_n), pigments, size = 1) +
  labs(x = NULL, y = lab_diatox) +
  scale_fill_acton(name = "% N (dry weight)", breaks = c(0.4, 0.6, 0.8, 1),
                   reverse = TRUE) +
  scale_color_acton(name = "% N (dry weight)", breaks = c(0.4, 0.6, 0.8, 1),
                    reverse = TRUE) +
  theme(legend.position = "top")

# surface plot
expand_grid(
  percent_n = gratia:::seq_min_max(pigments$percent_n, n = 100),
  year = gratia:::seq_min_max(pigments$year, n = 100),
  series = unique(pigments$series)) %>%
  mutate(time = year) %>%
  bind_cols(., predict(m_diatox_pn_ti, newdata = ., type = "expected")) %>%
  ggplot(aes(year, percent_n, fill = Estimate)) +
  geom_raster() +
  scale_x_continuous("Year CE", expand = c(0, 0)) +
  scale_y_continuous("% N (dry weight)", expand = c(0, 0)) +
  scale_fill_bamako(name = lab_diatox, limits = c(0, 200), na.value = "white") +
  theme(legend.position = "top")

#' dynamic coefficient models in `{mgcv}`
#' in `{mgcv}`, you can fit the term using `s(year, by = percent_n)`, or, more
#' specifically `s(year, by = percent_n, bs = "gp")`
m_diatox_pn$mgcv_model
#' in `{brms}`, you can also fit the term using `gp(year, by = percent_n, ...)`

# forecasting from dynamic models ----

air_passengers <-
  tibble(time = 1:length(AirPassengers),
         dec_date = time(AirPassengers) %>% as.numeric(),
         year = floor(dec_date),
         month = round((dec_date - year) * 12) + 1, # to ensure January is 1
         passengers = as.numeric(AirPassengers)) # in thousands

air_passengers <- mutate(air_passengers, lag_12_passengers = lag(passengers, 12))
data_train <- filter(air_passengers, year <= 1955) %>%
  filter(! is.na(lag_12_passengers)) # lagged value before 1st observation is NA
data_test <- filter(air_passengers, year > 1955)

ggplot(air_passengers, aes(dec_date, passengers, lty = year < 1955)) +
  geom_line() +
  geom_vline(xintercept = 1955, lty = "dashed") +
  labs(x = "Year CE", y = "International airline passengers (thousands)") +
  scale_linetype_manual("Dataset", values = c(3, 1), labels = c("Test", "Train"))

m_gam <- mvgam(formula = passengers ~ 0, # no error in observation process
               trend_formula = ~
                 log(lag_12_passengers) + # since we are on the log link scale
                 s(year, k = 5, bs = "tp") +
                 s(month, k = 10, bs = "cc"),
               trend_model = "None",
               noncentred = TRUE,
               knots = list(month = c(0.5, 12.5)),
               family = poisson(link = "log"),
               data = data_train,
               newdata = data_test, # calculate forecast while fitting
               chains = 4,
               burnin = 750,
               samples = 500,
               control = list(adapt_delta = 0.9),
               parallel = TRUE,
               silent = 2)

summary(m_gam) # check diagnostics
plot(m_gam)

# fit a GAM with an AR(1) term
# not fixing the AR(1) coefficient causes issues with PSIS diagnostics
get_mvgam_priors(formula = passengers ~ 0, # no error in observation process
                 trend_formula = ~
                   log(lag_12_passengers) +
                   s(year, k = 5, bs = "tp") +
                   s(month, k = 10, bs = "cc"),
                 trend_model = AR(p = 1),
                 data = data_train)

m_gam_ar <- mvgam(formula = passengers ~ 0, # no error in observation process
                  trend_formula = ~
                    log(lag_12_passengers) +
                    s(year, k = 5, bs = "tp") +
                    s(month, k = 10, bs = "cc"),
                  trend_model = AR(p = 1), # AR(1) model
                  # adding bounds for AR(1) range to improve diagnostics
                  priors = prior(normal(0.4, 0.01), class = ar1,
                                 lb = 0.35, ub = 0.45),
                  noncentred = TRUE, # use a noncentered AR(1) model
                  knots = list(month = c(0.5, 12.5)),
                  family = poisson(link = "log"),
                  data = data_train,
                  newdata = data_test, # calculate forecast while fitting
                  chains = 4,
                  burnin = 750,
                  samples = 1000,
                  control = list(adapt_delta = 0.95),
                  parallel = TRUE, silent = 2)

plot(m_gam_ar)
summary(m_gam_ar) # check diagnostics
mcmc_plot(m_gam_ar, type = "trace", variable = ".", regex = TRUE)

# plot diagnostics
# one point is slightly problematic
# not restricting the range of the AR(1) coef results in Pareto shapes > 5 due
# to very low effective sample sizes as the AR(1) fights with the smooth terms
layout(matrix(c(1, 1:3), ncol = 2, byrow = TRUE))
plot(m_gam_ar, type = "forecast")
plot(loo(m_gam_ar), diagnostic = "k")
plot(loo(m_gam_ar), diagnostic = "ESS") #' same as `diagnostic = "n_eff"`
layout(1)

#' Dynamic `{mvgam}` models contain draws for many quantities, all stored as MCMC
#' draws in an object of class `stanfit` in the `model_output` slot:
#' - `β` coefficients for linear predictor terms (called `b`)
#' - Family-specific shape/scale parameters:
#'    - `ϕ` for Negative Binomial,
#'    - `σ_obs` for Normal / LogNormal
#' - Trend-specific parameters:
#'    - `α` and `ρ` for GP trends,
#'    - `σ` and `ar1` for AR trends
#' - In-sample posterior predictions: `ypred`
#' - In-sample posterior trend estimates: `trend`

#' LV and LV_raw are latent variables from AR trend
#' ypred are predictions
#' mus are estimated means
#' trend

class(m_gam_ar$model_output)
m_gam_ar$model_output@model_pars # names of model parameters
m_gam_ar$model_output@par_dims # dimensions of parameter vectors/matrices
m_gam_ar$model_output # summary table of parameter samples
m_gam_ar$model_output@sim$samples[[1]][1:10, 1:8] # df of samples for chain 1

# view posterior draws of the trend
plot(m_gam_ar, type = "forecast") # with base plot
plot(forecast(m_gam_ar)) # with ggplot2
plot(forecast(m_gam_ar), realisations = TRUE) # CIs = summaries of realizations

# random draws from the posterior (NOTE: x axis is time since first observation)
plot(m_gam_ar, type = "forecast", realisations = TRUE, n_realisations = 10)
plot(m_gam_ar, type = "trend", realisations = TRUE, n_realisations = 10) +
  geom_vline(xintercept = nrow(data_train), lty = "dashed")
plot(m_gam_ar, type = "smooths", realisations = TRUE, n_realisations = 10,
     trend_effects = TRUE)

# draws summarized to credible intervals (NOTE: x axis is time since first obs)
plot(m_gam_ar, type = "forecast")
plot(m_gam_ar, type = "trend") +
  geom_vline(xintercept = nrow(data_train), lty = "dashed")
plot(m_gam_ar, type = "smooths", trend_effects = TRUE, )

# generate forecasts for up to end of 2026
# predicting later is useful if data are not available or too large to add
data_test
preds_2026 <- tibble(time = max(data_train$time) + 1:(12 * (2027 - 1963)),
                     dec_date = max(data_train$dec_date) + time/12,
                     year = floor(dec_date),
                     month = round((dec_date - year) * 12) + 1) %>%
  left_join(air_passengers %>% select(time, passengers, lag_12_passengers),
            by = "time")

#' need to predict with `for` loop since lag-12 values are not always avaiable
#' not the best way to predict: it does not include uncertainty in lagged values
tail(preds_2026)
for(i in which(is.na(preds_2026$passengers))) {
  if(is.na(preds_2026$lag_12_passengers[i])) {
    preds_2026$lag_12_passengers[i] <- preds_2026$passengers[i - 12]
  }
  
  preds_2026$passengers[i] <-
    predict(m_gam_ar, preds_2026[i, ], type = "expected")[, "Estimate"]
}
tail(preds_2026)

# to stop from calculating score (since values are predicted)
preds_2026 <- rename(preds_2026, estimate = passengers)

# model predicts that passengers will exceed 8 million by the end of 1972
print(preds_2026, n = 13) #' `lag_12_passengers` for row 13 is `estimate` for 1
max(preds_2026$dec_date)
plot_mvgam_fc(m_gam_ar, newdata = preds_2026, ylim = c(0, 8e9 / 1e3),
              realisations = TRUE)
plot_mvgam_fc(m_gam_ar, newdata = preds_2026, ylim = c(0, 8e9 / 1e3))
abline(v = 516, lty = "dashed")
filter(preds_2026, time == 516)

#' `plot(forecast(m_gam_ar, newdata = preds_2026))` fails with error:
#' `arguments imply differing number of rows: 780, 840`
#' function is adding the test data twice: `nrow(data_test)` is 60

# interpreting predictions ----
summary(m_gam_ar)

# - coefficients are often hard to interpret for GAMs, especially if non-gaussian
# - p-values are often too small for smooth terms bc they ignore uncertainty in
#   the smoothness parameter
# - predictions are more interpretable than coefficients

plot_preds <- function(.d, scale){
  .d %>%
    as.data.frame() %>%
    mutate(time = 1:n()) %>%
    ggplot() +
    geom_ribbon(aes(time, ymin = Q2.5, ymax = Q97.5), alpha = 0.3) +
    geom_line(aes(time, Estimate)) +
    labs(x = "Time", y = paste0("Estiamted effect (", scale, " scale)"))
}

# on link scale: for understanding coefficients of the linear predictor
predict(m_gam_ar, type = "link") %>% #' = `brms::posterior_linpred()`
  plot_preds(scale = "link") +
  ylim(c(0, 410))

# on expected scale: for understanding effects on the mean response
predict(m_gam_ar, type = "expected") %>% #' = `brms::posterior_epred()`
  plot_preds(scale = "expected") +
  ylim(c(0, 410))

# on response scale: for understanding effects on the individual observations
predict(m_gam_ar, type = "response") %>% #' = `brms::posterior_predict()`
  plot_preds(scale = "response") +
  ylim(c(0, 410))

#' can include process error with `process_error = TRUE`
predict(m_gam_ar, type = "response", process_error = TRUE) %>%
  plot_preds(scale = "response") +
  ylim(c(0, 410))

# the different scales with count data:
# - link: values are all real numbers (+ or -); scale is additive
# - expected: values are > 0; scale is multiplicative
# - response: values are > 0; scale is multiplicative

#' link-scale partial effects are centered around 0
#' allows to add intercept term to the partial effects
#' intercept = est. mean response, averaged across all smooths, on link scale
draw(m_gam_ar$trend_mgcv_model)
plot(forecast(m_gam_ar, type = "link"))
head(as.vector(forecast(m_gam_ar, type = "link")$forecasts$series1))

#' if `link = "log"` (does not apply if `link = "logit"`):
#' expected-scale partial effects are centered around 1
#' show the relative change in response with the predictor
#' allows to multiply intercept term by the partial effects
#' intercept = est. mean response, averaged across all smooths, on expected scale
draw(m_gam_ar$trend_mgcv_model, fun = exp) # relative change in passengers
draw(m_gam_ar$trend_mgcv_model, # partial change in passengers
     fun = \(y) exp(coef(m_gam_ar$trend_mgcv_model)["(Intercept)"]) *
       exp(y))
plot(forecast(m_gam_ar, type = "expected"))
head(as.vector(forecast(m_gam_ar, type = "expected")$forecasts$series1))

#' *NOTE:* in `{mgcv}` and `{gratia}`, predictions on the "response" scale are
#'         actually on the expected scale, since they only include uncertainty
#'         in the mean
#'         in `{mvgam}`, predictions on the response scale include uncertainty
#'         at the observation level, rather than just the mean
plot(forecast(m_gam_ar, type = "response"))
head(as.vector(forecast(m_gam_ar, type = "response")$forecasts$series1))

# useful for checking if the model has a good fit: can it simulate data well?
pp_check(m_gam_ar, type = "ribbon", ndraws = 100)
pp_check(m_gam_ar, type = "intervals", ndraws = 100)
pp_check(m_gam_ar, type = "scatter", ndraws = 9)
pp_check(m_gam_ar, type = "scatter_avg", ndraws = 100)
pp_check(m_gam_ar, type = "hist", ndraws = 8, bins = 10)
pp_check(m_gam_ar, type = "dens_overlay", ndraws = 100)
pp_check(m_gam_ar, type = "ecdf_overlay", ndraws = 100)

pp_check(m_gam_ar, type = "stat", ndraws = 1000, stat = "mean", binwidth = 1)
pp_check(m_gam_ar, type = "stat", ndraws = 1000, stat = "median", binwidth = 2.5)
pp_check(m_gam_ar, type = "stat", ndraws = 1000, stat = "sd", binwidth = 2.5)
pp_check(m_gam_ar, type = "stat", ndraws = 1000, stat = "var", binwidth = 250)

pp_check(m_gam_ar, type = "stat_2d", ndraws = 1000, stat = c("mean", "sd"))

# comparing models ----
layout(1:2)
acf(air_passengers$passengers)
pacf(air_passengers$passengers)
layout(1)

m_bad <- mvgam(formula = passengers ~ 0,
               trend_formula = ~ 1,
               trend_model = AR(p = 1), # AR(1) model
               noncentred = TRUE, # use a noncentered AR(1) model
               family = poisson(link = "log"),
               data = data_train, # calculate forecast while fitting
               newdata = data_test,
               chains = 4,
               burnin = 750,
               samples = 500,
               parallel = TRUE,
               silent = 2)

layout(matrix(c(1, 1:3), ncol = 2, byrow = TRUE))
plot(m_bad, type = "forecast") # the forecast is an AR(1) around the mean
plot(loo(m_bad), diagnostic = "k")
plot(loo(m_bad), diagnostic = "ESS") #' same as `diagnostic = "n_eff"`
# points closer to the estimated mean have better scores
layout(1)

# both models predict decently well for past data
plot(hindcast(m_gam_ar)) / plot(hindcast(m_bad))

# but they do not model the data the same way
plot_predictions(m_gam_ar, by = "time") / # GAM model "understands" the trends
  plot_predictions(m_bad, by = "time") # the AR model always assumes stationarity

# the naive model forecasts quite badly (reverts to the long-term mean)
layout(1:2)
plot_mvgam_fc(m_bad)
plot_mvgam_fc(m_gam_ar)
layout(1)

#' expected log posterior density: `mean(log(lik(y|model)))`
#' higher ELPD is better, but ELPD can be affected substantially by outliers
#' `m_gam` is best
#' `m_gam_ar` is worse
#' `m_bad` is clearly much worse
#' but ELPD can give too much importance to data in the tails because of it uses
#' the logged likelihood
#' as we'll see later, the score for `m_gam_ar` is penalized too much due to a
#' single uncertain prediction with `k_psis > 0.7`
#' `m_gam_ar` is actually the best model of the three
loo_compare(m_gam_ar, m_bad, m_gam)

#' can calculate ELPD quickly and easily with `{mvgam}`
#' requires `type = "link"`
score(forecast(m_gam_ar, type = "link"), score = "elpd")$series1 %>% head()

tibble(
  model_name = c("m_bad", "m_gam", "m_gam_ar"),
  elpd_data = map(model_name, \(m_n) {
    m <- get(m_n)
    bind_cols(data_test,
              score(forecast(m, type = "link"), score = "elpd")$series1)
  })) %>%
  unnest(elpd_data) %>%
  rename(ELPD = score) %>%
  ggplot(aes(eval_horizon, ELPD, color = model_name)) +
  geom_line(lwd = 1) +
  xlab("Forecast horizon (months)") +
  scale_color_highcontrast() +
  theme(legend.position = "top")

# predicting from new data (can't predict from bad model: no terms to plot!)
plot_predictions(m_gam_ar, condition = "month", points = 0.5)
plot_predictions(m_gam_ar, condition = "year", points = 0.5)

# rates of change: useful for link and expected scales, not for response scale
newd_slopes <- tibble(time = 1:100,
                      year = mean(data_train$year),
                      month = seq(0, 12, length.out = length(time)),
                      lag_12_passengers = mean(data_train$lag_12_passengers))

plot_grid(
  draw(m_gam_ar$trend_mgcv_model, select = "s(month)",
       data = tibble(lag_12_passengers = 0,
                     month = seq(0, 12, length.out = 400),
                     year = 0)),
  plot_slopes(m_gam_ar, variables = "month", by = "month", type = "link",
              newdata = newd_slopes) +
    geom_hline(yintercept = 0, linetype = "dashed"),
  ncol = 1)

plot_slopes(m_gam_ar, variables = "month", by = "month", type = "expected",
            newdata = newd_slopes) +
  geom_hline(yintercept = 0, linetype = "dashed")

# rates of change are too dramatic on the response scale because observations
# are too stochastic
plot_slopes(m_gam_ar, variables = "month", by = "month", type = "response",
            newdata = newd_slopes) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-1000, 1000))
