# Script for fitting PLS, RF and ensemble models
#1. Begin by sourcing necessary wrapper functions. Ensure you have internet connection.
#2. Read MIR data
#3. Read reference data

# install.packages(c("dplyr","downloader"), dependencies = TRUE)

library(dplyr)

library(downloader)

# Data set up
dir.create("~/Calibrations")

setwd("~/Calibrations")

download("https://www.dropbox.com/s/h4wwc0e5v3747fw/Workshop_data.zip?dl=1", "Workshop_data.zip", mode="wb")

unzip("Workshop_data.zip", overwrite = T)

raw <- read.csv("./Workshop_data/AfSIS_MIR_htsxt.csv") # Read HTS-xt spectral or Alpha or MPA tables

ref<-read.csv("./Workshop_data/AfSIS_reference.csv") # Read reference data

colnames(ref) <- c("SSN",colnames(ref[,-1]))

#Source wrapper functions for fitting models

rf_pls_url <- "https://raw.githubusercontent.com/asila/RF_PLS_Ensemble_methods/master/RF_PLS_optimal.R"

source_url(rf_pls_url)

# Create results outputpath

plswd <- "PLS_calibration_models")

dir.create(plswd)

wd <- plswd

set.seed(8892)

#Randomly split training and validation reference samples

m <- round(0.3*nrow(ref))

test <- sample(1:nrow(ref),m)

k <- which(raw$SSN %in% ref$SSN)

hout<-ref[test,]

calibrate(wd,raw,ref,hout, method="PLS", process = "none") # Use PLS or RF regression methods 

