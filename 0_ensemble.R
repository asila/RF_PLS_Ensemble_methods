# Script for fitting calibration models for infrared spectroscopy data. ----------------------------------------------------------------
# Author: Andrew Sila , April 2016. --------------------------------

# Begin by checking the required packages are installed. --------------------------------
is.installed <- function(anypkg){
	
  is.element(anypkg, installed.packages()[,1])
  
}

required.packages <- c("caret","soil.spec","doParallel","reshape")

installp <- which(!is.installed(required.packages) == TRUE)

#Install missing packages

if (length (installp) > 0){
	
install.packages(required.packages[installp])
}

library(caret)

library(soil.spec)

library(reshape)

library(doParallel)

registerDoParallel()

getDoParWorkers()

ensemble <- function(wd,infrared.data,reference.data,hout){
			
  # Start by setting working directory. 
  
  setwd(wd)
  
  si <- infrared.data
  
  ref <- reference.data
  
  hout <- hout

# 'NIR: Si spectra

prefx <- substr(colnames(si),1,1)[900]

wavenumbers <- round(as.numeric(gsub(prefx,"",colnames(si[,-1]))),1)

colnames(si) <- c("SSN",wavenumbers)

#preprocess sensor data

sid <- as.data.frame(trans(si[,-1])$trans)

si <- as.data.frame(cbind(as.vector(si[,1]),sid))

colnames(si) <- c("SSN", colnames(si[,-1]))


# combine transformed with reference
si.r <- na.omit(merge(ref,si))

#Select a random set for calibration and validation
if(hout!=0){
m <- round(hout*nrow(si.r))

set.seed(98752)
intest <- sample(1:nrow(si.r),m)

cal <- si.r [-intest, ]

val <- si.r [intest, ]
}

if(hout==0){
cal <- si.r 

val <- si.r
}

## check potential labels for response variable

ref.hd <- colnames(ref)[-1] # Remove metadata

for (q in 1 : length(ref.hd)){
	
lt <- cal[,ref.hd[q]]

lv <- val[,ref.hd[q]]



# Soil spectral features

mirt <- cal[,-c(1:ncol(ref))] # soil MIR

mirv <- val[,-c(1:ncol(ref))] # soil MIR

mir.a <- si[,-1] # all mir spectra for prediction 


# RF models ---------------------------------------------------------------

library(doParallel)

library(randomForest)




# Start doParallel to parallelize model fitting

mc <- makeCluster(detectCores())

registerDoParallel(mc)




# Control setup

set.seed(1385321)

tc <- trainControl(method = "cv", allowParallel = T)




# Tuning parameters

#tg <- expand.grid(mtry=seq(10, 150, by=10))

tg <- expand.grid(mtry=seq(10, 80, by=10))



# Fit model

mir.rfo <- train(mirt, lt,

                 preProc = c("center", "scale"),

                 method = "rf",

                 ntree = 101,

                 tuneGrid = tg,

                 trControl = tc)

print(mir.rfo)

rfo_mir <- predict(mir.rfo, mirv) ## predict validation set

rfo_mir.a <- predict(mir.rfo, mir.a) ## predict validation set

# Save the model 

saveRDS(mir.rfo, file = paste0(ref.hd[q],"_RFO.RDS"))


rm("mir.rfo")




stopCluster(mc)

#detach("package:randomForest", unload=TRUE)




# GBM models --------------------------------------------------------------

library(plyr)

library(gbm)




# Start doParallel to parallelize model fitting

mc <- makeCluster(detectCores())

registerDoParallel(mc)




# Control setup

set.seed(1385321)

tc <- trainControl(method = "repeatedcv", repeats=5, allowParallel = T)




# Tuning parameters

tg <- expand.grid(.n.trees=seq(10, 200, by=10), 

                  .interaction.depth = 10,

                  .shrinkage = 0.1,

                  .n.minobsinnode = 10)




# Fit model

mir.gbm <- train(mirt, lt, 

                 method = "gbm",

                 trControl = tc,

                 tuneGrid = tg)

print(mir.gbm)

gbm_mir <- predict(mir.gbm, mirv) ## predict validation set

gbm_mir.a <- predict(mir.gbm, mir.a) ## predict full set


# Save the model 

saveRDS(mir.gbm, file = paste0(ref.hd[q],"_gbm.RDS"))


rm("mir.gbm")




stopCluster(mc)

detach("package:gbm", unload=TRUE)




# PLS models --------------------------------------------------------------

library(pls)




# Start doParallel to parallelize model fitting

mc <- makeCluster(detectCores())

registerDoParallel(mc)




# Control setup

set.seed(1385321)

tc <- trainControl(method = "repeatedcv", repeats = 5, allowParallel = TRUE)




# Fit models

mir.pls <- train(mirt, log(ifelse(lt==0,NA,lt)),

                 preProc = c("center", "scale"),

                 method = "pls",

                 tuneGrid = expand.grid(ncomp=seq(2, 20, by=1)),

                 trControl = tc)

print(mir.pls)

pls_mir <- exp(predict(mir.pls, mirv)) ## predict validation set

pls_mir.a <- exp(predict(mir.pls, mir.a)) ## predict full set


# Save the model 

saveRDS(mir.pls ,file = paste0(ref.hd[q],"_pls.RDS"))


rm("mir.pls")




stopCluster(mc)

detach("package:pls", unload=TRUE)

# bartMachine models ------------------------------------------------------

options(java.parameters = "-Xmx5g")

library(bartMachine)




# Start doParallel to parallelize model fitting

mc <- makeCluster(detectCores())

registerDoParallel(mc)




# Control setup

tc <- trainControl(method = "cv", returnResamp = "all", allowParallel = T)




# Fit model

mir.bar <- train(mirt, lt,

                 method = "bartMachine", 

                 preProc = c("center", "scale"),

                 trControl = tc,

                 tuneLength = 2,
                 
                 serialize = TRUE,

                 seed = 123)

print(mir.bar)

bar_mir <- predict(mir.bar, mirv)

bar_mir.a <- predict(mir.bar, mir.a) ## predict full set


# Save the model 

saveRDS(mir.bar, file = paste0(ref.hd[q],"_bart.RDS"))


rm("mir.bar")


gcinfo(TRUE)

stopCluster(mc)

detach("package:bartMachine", unload=TRUE)


# cubist models 
library(Cubist)

# Start doParallel to parallelize model fitting

mc <- makeCluster(detectCores())

registerDoParallel(mc)
# Control setup

set.seed(1385321)

ctc <- cubistControl(unbiased = TRUE, rules = NA, sample = 0.3)



# Fit models
#committees <- round(c(0.25*nrow(mirt),0.35*nrow(mirt),0.45*nrow(mirt),0.5*nrow(mirt),0.75*nrow(mirt)),0)

mir.cubist <- train(mirt, lt,

                 preProc = c("center", "scale"),

                 method = "cubist")

                 #tuneGrid = expand.grid(.committees = committees,.neighbors = 0))
                 
                 #.committees = committees)

                 #trControl = ctc)

print(mir.cubist)

cubist_mir <- predict(mir.cubist, mirv) ## predict validation set

cubist_mir.a <- predict(mir.cubist, mir.a) ## predict full set

# Save the model 

saveRDS(mir.cubist, file = paste0(ref.hd[q],"_cubist.RDS"))

rm("mir.cubist")




stopCluster(mc)

detach("package:Cubist", unload=TRUE)

# Model stacking setup ----------------------------------------------------

pmirv <- as.data.frame(cbind(lv, rfo_mir, gbm_mir, pls_mir,bar_mir,cubist_mir))


pmirv.a <- as.data.frame(cbind(rfo_mir.a, gbm_mir.a, pls_mir.a,bar_mir.a,cubist_mir.a))

names(pmirv) <- c("L", "RFO", "GBM", "PLS","BART","CUBIST")

names(pmirv.a) <- c("RFO", "GBM", "PLS","BART","CUBIST")



# Remove extraneous objects from memory -----------------------------------

# rm(list=setdiff(ls(), pmirv"))




# Model stacking ----------------------------------------------------------

library(glmnet)




# Start doParallel to parallelize model fitting

mc <- makeCluster(detectCores())

registerDoParallel(mc)




# Control setup

set.seed(1385321)

tc <- trainControl(method = "cv", allowParallel = T)




# MIR model stack

set.seed(1385321)

mir.ens <- train(L ~ ., data = pmirv,

                 method = "glmnet",

                 family = "gaussian",

                 trControl = tc)

print(mir.ens)

saveRDS(mir.ens, file = paste0(ref.hd[q],".RDS"))

ens_mir <- as.data.frame(predict(mir.ens, pmirv))

names(ens_mir) <- c("ENS")

ens_mir.a <- as.data.frame(predict(mir.ens, pmirv.a))

names(ens_mir.a) <- c("ENS")

pmirv <- cbind(pmirv, ens_mir)

pmirv.a <- cbind(pmirv.a, ens_mir.a)



stopCluster(mc)




# Write data files --------------------------------------------------------

write.csv(pmirv, paste0(ref.hd[q],"_pmirv.csv"), row.names=F)

write.csv(pmirv.a, paste0(ref.hd[q],"_pmirv.a.csv"), row.names=F)

png(file = paste0(ref.hd[q],".png"), height = 400, width = 600)
# Prediction plots --------------------------------------------------------

par(mfrow=c(2,3), mar=c(5,4.5,1,1))




# MIR predictions # note that x & y axis limits will need to be adjusted

lmin <- 0

lmax <- max(pmirv$L)

plot(L ~ RFO, pmirv, xlim=c(min(pmirv$L,pmirv$RFO), max(pmirv$L,pmirv$RFO)), ylim=c(min(pmirv$L,pmirv$RFO), max(pmirv$L,pmirv$RFO)), xlab = "RFO prediction", ylab = "Observed", cex.lab=1.3)

abline(c(0,1), col="red")

plot(L ~ GBM, pmirv, xlim=c(min(pmirv$L,pmirv$GBM), max(pmirv$L,pmirv$GBM)), ylim=c(min(pmirv$L,pmirv$GBM), max(pmirv$L,pmirv$GBM)), xlab = "GBM prediction", ylab = "Observed", cex.lab=1.3)

abline(c(0,1), col="red")

plot(L ~ PLS, pmirv, xlim=c(min(pmirv$L,pmirv$PLS), max(pmirv$L,pmirv$PLS)), ylim=c(min(pmirv$L,pmirv$PLS), max(pmirv$L,pmirv$PLS)), xlab = "PLS prediction", ylab = "Observed", cex.lab=1.3)

abline(c(0,1), col="red")

plot(L ~ CUBIST, pmirv, xlim=c(min(pmirv$L,pmirv$CUBIST), max(pmirv$L,pmirv$CUBIST)), ylim=c(min(pmirv$L,pmirv$CUBIST), max(pmirv$L,pmirv$CUBIST)), xlab = "CUBIST prediction", ylab = "Observed", cex.lab=1.3)

abline(c(0,1), col="red")

plot(L ~ BART, pmirv, xlim=c(min(pmirv$L,pmirv$BART), max(pmirv$L,pmirv$BART)), ylim=c(min(pmirv$L,pmirv$BART), max(pmirv$L,pmirv$BART)), xlab = "BART prediction", ylab = "Observed", cex.lab=1.3)

abline(c(0,1), col="red")

# Ensemble predictions 

plot(L ~ ENS, pmirv, xlim=c(min(pmirv$L,pmirv$ENS), max(pmirv$L,pmirv$ENS)), ylim=c(min(pmirv$L,pmirv$ENS), max(pmirv$L,pmirv$ENS)), xlab = "Model ensemble prediction", ylab = "Observed", cex.lab=1.3,col = "blue",pch = 16)

abline(c(0,1), col="red")
dev.off()
}
}