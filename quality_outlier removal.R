setwd("D:\\2017_data\\MLINGANO\\mlingano_ref\\result\\pls")

list.files()

pred <- read.csv("All predictions.csv")

ref<-read.csv("D:/2017_data/MLINGANO/mlingano_ref/All Mli_ref combined_renamed_ssns.csv")[,-c(2:4)]

colnames(ref) <- gsub("EC.S", "EC",colnames(ref))

rp <- merge(ref, pred, by ="SSN")


colnames(rp) <- gsub(".x", "_measured", colnames(rp))

colnames(rp) <- gsub(".y", "_predicted", colnames(rp))

colnames(rp) <- gsub("EC.S", "EC",colnames(rp))

vs <- colnames(ref[,-1])

out.all <- NULL

rpp <- NULL

for (i in 1:length(vs)){

s <- grep(vs[i],colnames(rp))

rpi <- rp[,s]

res <- rpi[,2] - rpi[,1]

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

write.table(rpp, file = "D:/2017_data/MLINGANO/mlingano_ref/All Mli_ref combined_renamed_ssns_excludes_outliers.csv", sep = ",", row.names = FALSE)
u.out