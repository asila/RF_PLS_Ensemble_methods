# Script for fitting PLS, RF and ensemble models
#1. Begin by sourcing necessary wrapper functions. Ensure you have internet connection.
#2. Read MIR data
#3. Read reference data

# install.packages(c("dplyr","downloader"), dependencies = TRUE)

library(dplyr)

library(downloader)

# Data set up
dir.create("./Calibrations")

setwd("./Calibrations")

download("https://www.dropbox.com/s/76beos74ufvgiyu/iiss_cal.zip?dl=1", "iiss_cal.zip", mode="wb")

unzip("iiss_cal.zip", overwrite = T)

raw <- read.csv("./IISS_mir spectra.csv") # Read HTS-xt spectral or Alpha or MPA tables

ref<-read.csv("./IISS_cn data.csv") # Read reference data

colnames(ref) <- c("SSN",colnames(ref[,-1]))

#Source wrapper functions for fitting models

rf_pls_url <- "https://raw.githubusercontent.com/asila/RF_PLS_Ensemble_methods/master/RF_PLS_optimal.R"

source_url(rf_pls_url)

# Create results outputpath

plswd <- "PLS_calibration_models"

dir.create(plswd)

wd <- plswd

set.seed(8892)

#Randomly split training and validation reference samples

m <- round(0.3*nrow(ref))

test <- sample(1:nrow(ref),m)

k <- which(raw$SSN %in% ref$SSN)

hout<-ref[test,]

calibrate(wd,raw,ref,hout, method="RF", process = "derivative") # Use PLS or RF regression methods 

