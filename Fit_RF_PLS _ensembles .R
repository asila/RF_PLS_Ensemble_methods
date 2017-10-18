#PLS or RF
#1. Setwd
#2. Read MIR data
#3. Read reference data

library(dplyr)

# download soil data
dir.create("~/AfSIS")

download("https://www.dropbox.com/s/h4wwc0e5v3747fw/Workshop_data.zip?dl=1", "Workshop_data.zip", mode="wb")

unzip("Workshop_data.zip", overwrite=T)

raw <- read.csv("./Workshop_data/AfSIS_MIR_htsxt.csv") # Read HTS-xt spectral or Alpha or MPA tables

ref<-read.csv("./Workshop_data/AfSIS_reference.csv") # Read reference data

colnames(ref) <- c("SSN",colnames(ref[,-1]))

source('~/Dropbox/temp/share/calibrations/RF_PLS_optimal.R', chdir = TRUE)

# Create results outputpath

dir.create(paste0("Infrared_calibration_models_",Sys.time()))

wd <- paste0("Infrared_calibration_models_",Sys.time())

set.seed(8892)

#Randomly split training and validation reference samples

m <- round(0.3*nrow(ref))

test <- sample(1:nrow(ref),m)

k <- which(raw$SSN %in% ref$SSN)

hout<-ref[test,]

calibrate(wd,raw,ref,hout, method="PLS") # Use PLS or RF regression methods 

#Reset working directory

setwd("~/AfSIS") # Root directory

dir.create(paste0("Infrared_calibration_models_",Sys.time()))

wd <- paste0("Infrared_calibration_models_",Sys.time())

calibrate(wd,raw,ref,hout, method="RF")

## Ensemble
#source('~/local.github/calibrations/Ensembles/0_ensemble.R', chdir = TRUE)
source("/Users/andrewsila/local.github/calibrations/Ensembles/0_ensemble.R")

setwd("~/AfSIS") # Root directory

dir.create(paste0("Ensemble_Infrared_calibration_models_",Sys.time()))

wd <- paste0("Ensemble_Infrared_calibration_models_",Sys.time())

set.seed(809892)

ensemble(wd,raw,ref,hout = 0.3)

