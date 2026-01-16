if(FALSE) {
  source('setup.R') #' run once before running `packages.R` or any other scripts
}

# attach necessary packages
library('dplyr')     #' for data wrangling and tidying
library('tidyr')     #' for data wrangling and tidying
library('lubridate') #' for working with dates
library('ggplot2')   #' for fancy plots
library('mvgam')     #' for fitting models; uses `mgcv`, `brms` and `cmdstan`

theme_set(theme_bw() + theme(text = element_text(face = 'bold')))
