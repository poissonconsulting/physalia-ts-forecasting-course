if(FALSE) {
  source('setup.R') #' run once before running `packages.R` or any other scripts
}

# attach necessary packages
library('openxlsx')  #' for importing xlsx files
library('gratia')    #' for plotting GAMs
library('ctmm')      #' for continuous-time stochastic models
library('dplyr')     #' for data wrangling and tidying
library('tidyr')     #' for data wrangling and tidying
library('purrr')     #' for functional programming
library('lubridate') #' for working with dates
library('ggplot2')   #' for fancy plots
library('mvgam')     #' for fitting models; uses `mgcv`, `brms` and `cmdstan`
library('khroma')    #' for color palettes
library('cowplot')   #' for multi-panel plots

theme_set(theme_bw() + theme(text = element_text(face = 'bold')))
