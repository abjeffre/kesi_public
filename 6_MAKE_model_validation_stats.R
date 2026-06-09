#########################################
######### Bio ECON  #####################
library(readr)
library(abind)
library(rethinking)
library(posterior)
library(ggplot2)
library(cowplot)

set_project_wd <- function(folder){
  user=Sys.info()[[6]]
  if(user=="jeffrey_andrews") setwd(paste0("C:/Users/jeffrey_andrews/OneDrive/Documents/", folder))
  else if(user=="Jeff") setwd(paste0("C:/Users/Jeff/OneDrive/Documents/", folder))
  else if(user == 'jeffr')  setwd(paste0("C:/Users/jeffr/OneDrive/Documents/", folder))
  else if(user == 'unknown')  setwd(paste0("~/", folder))
}
set_project_wd("Bio_econ")
source("/functions/utility.R")


################################################################
########## RUN STAN MODELS FOR CLIMATE ON EARNINGS #############

source("/code/sweeps/weather_on_earnings_recover.R")
source("/code/sweeps/earnings_on_kesi_recover.R")
