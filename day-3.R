source('packages.R') # attach necessary packages

## Forecasting from dynamic models ---

##' applications of state-space models:
##' Kalman filter and Apollo missions: `https://doi.org/10.1109/MCS.2010.936465`

air_passengers <-
  tibble(time = 1:length(AirPassengers),
         dec_date = time(AirPassengers) %>% as.numeric(),
         year = floor(dec_date),
         month = round((dec_date - year) * 12) + 1, # to ensure January is 1
         passengers = as.numeric(AirPassengers)) # in thousands

data_train <- filter(air_passengers, year <= 1955)
data_test <- filter(air_passengers, year > 1955)

m_gam_ar <- mvgam(formula = passengers ~ 0, # no error in observation process
                  trend_formula = ~
                    s(year, k = 5, bs = 'tp') +
                    s(month, k = 10, bs = 'cc'),
                  trend_model = AR(p = 1), # AR(1) model
                  noncentred = TRUE, # use a noncentered AR(1) model
                  knots = list(month = c(0.5, 12.5)),
                  family = poisson(),
                  data = data_train,
                  newdata = data_test,
                  chains = 4,
                  burnin = 750,
                  samples = 500,
                  control = list(max_treedepth = 20, adapt_delta = 0.95),
                  parallel = TRUE)

##' Dynamic `mvgam` models contain draws for many quantities, all stored as MCMC
##' draws in an object of class `stanfit` in the `model_output` slot:
##' - `β` coefficients for linear predictor terms (called `b`)
##' - Family-specific shape/scale parameters:
##'    - `ϕ` for Negative Binomial,
##'    - `σ_obs` for Normal / LogNormal
##' - Trend-specific parameters:
##'    - `α` and `ρ` for GP trends,
##'    - `σ` and `ar1` for AR trends
##' - In-sample posterior predictions: `ypred`
##' - In-sample posterior trend estimates: `trend`

class(m_gam_ar$model_output)
m_gam_ar$model_output@model_pars # names of model parameters
m_gam_ar$model_output@par_dims # dimensions of parameter vectors/matrices
m_gam_ar$model_output # summary table of parameter samples
m_gam_ar$model_output@sim$samples[[1]][1:10, 1:8] # df of samples for chain 1

# view posterior draws of the trend
plot(m_gam_ar, type = 'forecast')
plot(forecast(m_gam_ar))
plot(forecast(m_gam_ar), realisations = TRUE)

# random draws from the posterior
plot(m_gam_ar, type = 'trend', realisations = TRUE, n_realisations = 10) +
  geom_vline(xintercept = max(data_train), lty = 'dashed')

# draws summarized to credible intervals
plot(m_gam_ar, type = 'trend') +
  geom_vline(xintercept = nrow(data_train), lty = 'dashed')

## generate forectasts for additional new data
data_test
new_data <- tibble(time = max(data_train$time) + 1:(12 * (2027 - 1963)),
                   dec_date = max(data_train$dec_date) + time/12,
                   year = floor(dec_date),
                   month = round((dec_date - year) * 12) + 1)
new_data
max(new_data$dec_date)
plot_mvgam_fc(m_gam_ar, newdata = new_data)

#' `plot(forecast(m_gam_ar, newdata = new_data))` fails with error:
#' `arguments imply differing number of rows: 792, 852`
##' function is adding the test data twice: `nrow(data_test)` is 60

## Point-based forecast evaluation

## Probabilistic forecast evaluation

## Bayesian posterior predictive checks

#' ADD SSMs?
##' `O_t ~ MVN(Y_proc, s_obs)`
##' `Y_proc ~ MVN(mu_proc, s_proc)`
##' `Mu_proc = GAM`
