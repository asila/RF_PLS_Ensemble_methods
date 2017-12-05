library(dplyr)
pred <- read.csv("~/studies/India/ICRAF_data/data/Sentinels/RF_calibration_models2017-11-24 08:33:20/All predictions.csv")

ref <- read.csv("")

colnames(ref) <- c("SSN",colnames(ref[,-1]))

ref <- ref%>% group_by(SSN) %>% summarise_each(funs(mean), -SSN)

ref <- ref[,-4]

#Second round
ref <- read.csv("~/studies/India/Africa_India/AfSIS_India_excludes_outliers1.csv")


rp <- merge(ref, pred, by ="SSN")

colnames(rp) <- gsub(".x", "_measured", colnames(rp))

colnames(rp) <- gsub(".y", "_predicted", colnames(rp))

c#olnames(rp) <- gsub("EC.S", "EC",colnames(rp))

vs <- colnames(ref[,-1])

out.all <- NULL

rpp <- NULL

for (i in 1:length(vs)){

s <- grep(vs[i],colnames(rp))

rpi <- rp[,s]

res <- na.omit(rpi[,2] - rpi[,1])

res.m <-  mean(res)

res.d <- sd(res)

out <- which(res < (res.m - 3*res.d)|res > (res.m + 3*res.d))

rpi[out,1] <- ""

rpp <- cbind(rpp, rpi[,1])

out.all <- c(out.all, out)

}

#u.out <- unique(out.all)

rpp <- cbind(as.vector(rp[,1]), rpp)

colnames(rpp) <- c("SSN", vs)


#Save a new table with reference data which excludes reference values for samples identified with large outliers

write.table(rpp, file = "~/studies/India/Africa_India/AfSIS_India_excludes_outliers2.csv", sep = ",", row.names = FALSE)

rpp <- read.csv("~/studies/India/Africa_India/AfSIS_India_excludes_outliers2.csv")

str(rpp)