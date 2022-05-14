setwd("/Users/staceyvann/Desktop/2022Spring/AppliedSequencing/Final")

hi_ctrl1_t <- read.table("counts/PRAMEhi-ctrl1_featurecounts", header = FALSE, sep = "", dec = ".") 
colnames(hi_ctrl1_t) <- hi_ctrl1_t[1,]
hi_ctrl1_t <- hi_ctrl1_t[-1,]
hi_ctrl1_t <- hi_ctrl1_t[,c(1,7)]
colnames(hi_ctrl1_t)[2] <- "hi_ctrl1" 
hi_ctrl1_t[,1] <- gsub("\\..*", "", hi_ctrl1_t[,1])
hi_ctrl1_rm0 <- hi_ctrl1_t[!(apply(hi_ctrl1_t, 1, function(y) any(y == 0))),]

hi_ctrl2_t <- read.table("counts/PRAMEhi-ctrl2_featurecounts", header = FALSE, sep = "", dec = ".") 
colnames(hi_ctrl2_t) <- hi_ctrl2_t[1,]
hi_ctrl2_t <- hi_ctrl2_t[-1,]
hi_ctrl2_t <- hi_ctrl2_t[,c(1,7)]
colnames(hi_ctrl2_t)[2] <- "hi_ctrl2" 
hi_ctrl2_t[,1] <- gsub("\\..*", "", hi_ctrl2_t[,1])
hi_ctrl2_rm0 <- hi_ctrl2_t[!(apply(hi_ctrl2_t, 1, function(y) any(y == 0))),]

hi_ctrl3_t <- read.table("counts/PRAMEhi-ctrl3_featurecounts", header = FALSE, sep = "", dec = ".") 
colnames(hi_ctrl3_t) <- hi_ctrl3_t[1,]
hi_ctrl3_t <- hi_ctrl3_t[-1,]
hi_ctrl3_t <- hi_ctrl3_t[,c(1,7)]
colnames(hi_ctrl3_t)[2] <- "hi_ctrl3" 
hi_ctrl3_t[,1] <- gsub("\\..*", "", hi_ctrl3_t[,1])
hi_ctrl3_rm0 <- hi_ctrl3_t[!(apply(hi_ctrl3_t, 1, function(y) any(y == 0))),]

hi_siRNA1_t <- read.table("counts/PRAMEhi-siRNA1_featurecounts", header = FALSE, sep = "", dec = ".") 
colnames(hi_siRNA1_t) <- hi_siRNA1_t[1,]
hi_siRNA1_t <- hi_siRNA1_t[-1,]
hi_siRNA1_t <- hi_siRNA1_t[,c(1,7)]
colnames(hi_siRNA1_t)[2] <- "hi_siRNA1" 
hi_siRNA1_t[,1] <- gsub("\\..*", "", hi_siRNA1_t[,1])
hi_siRNA1_rm0 <- hi_siRNA1_t[!(apply(hi_siRNA1_t, 1, function(y) any(y == 0))),]

hi_siRNA2_t <- read.table("counts/PRAMEhi-siRNA2_featurecounts", header = FALSE, sep = "", dec = ".") 
colnames(hi_siRNA2_t) <- hi_siRNA2_t[1,]
hi_siRNA2_t <- hi_siRNA2_t[-1,]
hi_siRNA2_t <- hi_siRNA2_t[,c(1,7)]
colnames(hi_siRNA2_t)[2] <- "hi_siRNA2" 
hi_siRNA2_t[,1] <- gsub("\\..*", "", hi_siRNA2_t[,1])
hi_siRNA2_rm0 <- hi_siRNA2_t[!(apply(hi_siRNA2_t, 1, function(y) any(y == 0))),]

hi_siRNA3_t <- read.table("counts/PRAMEhi-siRNA3_featurecounts", header = FALSE, sep = "", dec = ".") 
colnames(hi_siRNA3_t) <- hi_siRNA3_t[1,]
hi_siRNA3_t <- hi_siRNA3_t[-1,]
hi_siRNA3_t <- hi_siRNA3_t[,c(1,7)]
colnames(hi_siRNA3_t)[2] <- "hi_siRNA3" 
hi_siRNA3_t[,1] <- gsub("\\..*", "", hi_siRNA3_t[,1])
hi_siRNA3_rm0 <- hi_siRNA3_t[!(apply(hi_siRNA3_t, 1, function(y) any(y == 0))),]

PRAMEhi <- merge(hi_ctrl1_rm0,hi_ctrl2_rm0, all = TRUE)
PRAMEhi <- merge(PRAMEhi,hi_ctrl3_rm0, all = TRUE)
PRAMEhi <- merge(PRAMEhi,hi_siRNA1_rm0, all = TRUE)
PRAMEhi <- merge(PRAMEhi,hi_siRNA2_rm0, all = TRUE)
PRAMEhi <- merge(PRAMEhi,hi_siRNA3_rm0, all = TRUE)

PRAMEhi <- cbind(PRAMEhi$Geneid,PRAMEhi[,5:7],PRAMEhi[,2:4])
colnames(PRAMEhi) <- c("Geneid","hi_ctrl1","hi_ctrl2","hi_ctrl3","hi_siRNA1","hi_siRNA2","hi_siRNA3")

write.table(PRAMEhi, "/Users/staceyvann/Desktop/2022Spring/AppliedSequencing/Final/PRAMEhi.csv",row.names = FALSE,col.names = TRUE,sep = ",")

