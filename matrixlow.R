setwd("/Users/staceyvann/Desktop/2022Spring/AppliedSequencing/Final")

lo_ctrl1_t <- read.table("counts/PRAMElo-ctrl1_featurecounts", header = FALSE, sep = "", dec = ".") 
colnames(lo_ctrl1_t) <- lo_ctrl1_t[1,]
lo_ctrl1_t <- lo_ctrl1_t[-1,]
lo_ctrl1_t <- lo_ctrl1_t[,c(1,7)]
colnames(lo_ctrl1_t)[2] <- "lo_ctrl1" 
lo_ctrl1_t[,1] <- gsub("\\..*", "", lo_ctrl1_t[,1])
lo_ctrl1_rm0 <- lo_ctrl1_t[!(apply(lo_ctrl1_t, 1, function(y) any(y == 0))),]

lo_ctrl2_t <- read.table("counts/PRAMElo-ctrl2_featurecounts", header = FALSE, sep = "", dec = ".") 
colnames(lo_ctrl2_t) <- lo_ctrl2_t[1,]
lo_ctrl2_t <- lo_ctrl2_t[-1,]
lo_ctrl2_t <- lo_ctrl2_t[,c(1,7)]
colnames(lo_ctrl2_t)[2] <- "lo_ctrl2" 
lo_ctrl2_t[,1] <- gsub("\\..*", "", lo_ctrl2_t[,1])
lo_ctrl2_rm0 <- lo_ctrl2_t[!(apply(lo_ctrl2_t, 1, function(y) any(y == 0))),]

lo_ctrl3_t <- read.table("counts/PRAMElo-ctrl3_featurecounts", header = FALSE, sep = "", dec = ".") 
colnames(lo_ctrl3_t) <- lo_ctrl3_t[1,]
lo_ctrl3_t <- lo_ctrl3_t[-1,]
lo_ctrl3_t <- lo_ctrl3_t[,c(1,7)]
colnames(lo_ctrl3_t)[2] <- "lo_ctrl3" 
lo_ctrl3_t[,1] <- gsub("\\..*", "", lo_ctrl3_t[,1])
lo_ctrl3_rm0 <- lo_ctrl3_t[!(apply(lo_ctrl3_t, 1, function(y) any(y == 0))),]

lo_oe1_t <- read.table("counts/PRAMElo-oe1_featurecounts", header = FALSE, sep = "", dec = ".") 
colnames(lo_oe1_t) <- lo_oe1_t[1,]
lo_oe1_t <- lo_oe1_t[-1,]
lo_oe1_t <- lo_oe1_t[,c(1,7)]
colnames(lo_oe1_t)[2] <- "lo_oe1" 
lo_oe1_t[,1] <- gsub("\\..*", "", lo_oe1_t[,1])
lo_oe1_rm0 <- lo_oe1_t[!(apply(lo_oe1_t, 1, function(y) any(y == 0))),]

lo_oe2_t <- read.table("counts/PRAMElo-oe2_featurecounts", header = FALSE, sep = "", dec = ".") 
colnames(lo_oe2_t) <- lo_oe2_t[1,]
lo_oe2_t <- lo_oe2_t[-1,]
lo_oe2_t <- lo_oe2_t[,c(1,7)]
colnames(lo_oe2_t)[2] <- "lo_oe2" 
lo_oe2_t[,1] <- gsub("\\..*", "", lo_oe2_t[,1])
lo_oe2_rm0 <- lo_oe2_t[!(apply(lo_oe2_t, 1, function(y) any(y == 0))),]

lo_oe3_t <- read.table("counts/PRAMElo-oe3_featurecounts", header = FALSE, sep = "", dec = ".") 
colnames(lo_oe3_t) <- lo_oe3_t[1,]
lo_oe3_t <- lo_oe3_t[-1,]
lo_oe3_t <- lo_oe3_t[,c(1,7)]
colnames(lo_oe3_t)[2] <- "lo_oe3" 
lo_oe3_t[,1] <- gsub("\\..*", "", lo_oe3_t[,1])
lo_oe3_rm0 <- lo_oe3_t[!(apply(lo_oe3_t, 1, function(y) any(y == 0))),]


PRAMElo <- merge(lo_ctrl1_rm0,lo_ctrl2_rm0, all = TRUE)
PRAMElo <- merge(PRAMElo,lo_ctrl3_rm0, all = TRUE)
PRAMElo <- merge(PRAMElo,lo_oe1_rm0, all = TRUE)
PRAMElo <- merge(PRAMElo,lo_oe2_rm0, all = TRUE)
PRAMElo <- merge(PRAMElo,lo_oe3_rm0, all = TRUE)

write.table(PRAMElo, "/Users/staceyvann/Desktop/2022Spring/AppliedSequencing/Final/PRAMElo.csv",row.names = FALSE,col.names = TRUE,sep = ",")

