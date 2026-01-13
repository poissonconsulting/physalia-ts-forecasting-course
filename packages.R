if(FALSE) { # installations and checks that don't need to be run each time
  # check that R is (relatively) up to date
  if(version$major < 4) {
    stop(paste0('Please update R. Your current version (', version$major, '.',
                version$minor), ') is too old.')
  } else if (as.numeric(version$minor) <= 5) {
    warning(paste0('Your R version (', version$major, '.', version$minor,
                   ') is out of date. You may want to update it.'))
  }
  
  #' update any installed `R` packages
  update.packages(ask = FALSE, checkBuilt = TRUE)
  
  #' install the development version of `{brms}` and its dependencies
  install.packages("remotes")
  remotes::install_github("paul-buerkner/brms", dependencies = TRUE)
  
  #' install `{mvgam}` and a few other packages we will use for plotting
  install.packages(c("mvgam", "gratia", "tidybayes"))
  
  #' install and check `Stan`
  #' *NOTE:* `{cmdstanr}` is not on CRAN, so it's more up to date because it's
  #'         easier to update often. it's also more lightweight than `{rstan}`
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
  
  library('cmdstanr')
  check_cmdstan_toolchain(fix = TRUE)
  install_cmdstan(cores = 2)
  
  #' check that `{mvgam}` is working properly
  library('mvgam')
  simdat <- sim_mvgam()
  mod <- mvgam(y ~ s(season, bs = 'cc', k = 5) +
                 s(time, series, bs = 'fs', k = 8),
               data = simdat$data_train)
  plot(mod$mgcv_model)
  mod$adapt_delta
}

# attach necessary packages
library('dplyr')     #' for data wrangling and tidying
library('tidyr')     #' for data wrangling and tidying
library('lubridate') #' for working with dates
library('ggplot2')   #' for fancy plots
library('mvgam')     #' for fitting models; uses `mgcv`, `brms` and `cmdstan`

theme_set(theme_bw() + theme(text = element_text(face = 'bold')))
