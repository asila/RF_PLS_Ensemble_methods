# Script for making predictions from random forest models
# Author: Andrew Sila
# email: a.sila@cgiar.org
#
# Mantainer: Gard okello
# email: g.okello@cgiar.org
#
#
# Begin by installing the required packages if not already done
#
# install.packages(c("soil.spec","caret","downlaoder")) # remove hashtag before install.packages to do the installation
 
library(soil.spec)

library(caret)

library(downloader)

dir.create("./RF_predictions")

setwd("./RF_predictions")

download("https://www.dropbox.com/s/k717jqr2pef1h2r/Data.zip?dl=1", "Data.zip",mode = "wb")

unzip("Data.zip")

raw <- read.csv("mir spectra.csv")# New spectral table can be used here!

codes <- read.csv("sample codes.csv") # Read file with field sampling details

wavenumbers <- as.numeric(gsub("w","",colnames(raw[,-1])))

colnames(raw) <- c("SSN", wavenumbers)

raw0 <- as.matrix(raw[,-1])

spec <- trans(raw0)$trans

#Read RF models

dir.create("./RF_models")

setwd("./RF_models")

download("https://www.dropbox.com/s/60xoafveccovqrq/Models.zip?dl=1", "Models.zip",mode = "wb")

unzip("Models.zip")

#pr <- substr(colnames(spec),1,1)

#colnames(spec) <- paste0(pr,colnames(spec))

#Read all models and predict

mods <- list.files(pattern = ".Rds")

para <- gsub(".Rds","", mods)

pred.all <- NULL

for (m in 1:length(mods)){
	
	rfm.m <- readRDS(mods[m])
	
	mod.m <- round(predict(rfm.m, spec),3)
	
	pred.all <- cbind(pred.all,mod.m)
}

colnames(pred.all) <- para

pred.ssn <- cbind(as.vector(raw[,1]),pred.all)

colnames(pred.ssn) <- c("SSN", colnames(pred.ssn[,-1]))

pred.ssn[1:9,]

#To know where the predicted is saved use-  getwd(); otherwise set a new working directory path

write.table(pred.ssn, file = "MIR_predicted_results.csv", sep = ",", row.names = FALSE)

#Read the predictions back and merge with sample codes
 pred.ssn <- read.csv("/Volumes/TBO/VS/Ug_vs2 mir predicted.csv")
 
 # Link the MIR predictions to the field details table
  
 predicted <- merge(codes, pred.ssn)
 
 #To know where the predicted matched with field details is saved use-  getwd(); otherwise set a new working directory path

 write.table(predicted, file = "MIR_predicted_results_with_field_codes.csv", sep =",", row.names = FALSE)
