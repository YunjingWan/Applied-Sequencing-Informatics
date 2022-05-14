library(dplyr)
library(enrichplot)
library(openxlsx)
library(ggvenn)
library(qpcR)
library(ggplot2)
library(ggrepel)

setwd("/Users/staceyvann/Desktop/2022Spring/AppliedSequencing/Final")

lo <- read.csv("/Users/staceyvann/Desktop/2022Spring/AppliedSequencing/Final/DEGdeseq2_lo.csv")
hi <- read.csv("/Users/staceyvann/Desktop/2022Spring/AppliedSequencing/Final/DEGdeseq2_hi.csv")
lo <- lo[,c(1,10,14)]
hi <- hi[,c(1,10,14)]
hi_down <- subset(hi, Differential_Gene == "Down-Regulated")
hi_up <- subset(hi, Differential_Gene == "Up-Regulated")
lo_down <- subset(lo, Differential_Gene == "Down-Regulated")
lo_up <- subset(lo, Differential_Gene == "Up-Regulated")
lo_up[9,2] <- "SNORD3A"

x = list(hi_up = hi_up$symbol, lo_down = lo_down$symbol)
y = list(hi_down = hi_down$symbol, lo_up = lo_up$symbol)

x1 <- ggvenn(
  x,
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)
y1 <- ggvenn(
  y,
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)

x1+y1

#matched genes in detail
matched <- match(hi_down$symbol,lo_up$symbol,nomatch = NA)
posmatched <- na.omit(lo_up[matched,])
