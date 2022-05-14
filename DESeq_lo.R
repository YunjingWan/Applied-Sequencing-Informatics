library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("vsn")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("genefilter")
library("biomaRt")
library("IHW")
library("ggplot2")
library("apeglm")
library(ggrepel)
library("clusterProfiler")
library(dplyr)
library(enrichplot)


setwd("/Users/staceyvann/Desktop/2022Spring/AppliedSequencing/Final")
CountTable <- read.table("/Users/staceyvann/Desktop/2022Spring/AppliedSequencing/Final/PRAMElo.csv", header=TRUE,sep = ",") 
rownames(CountTable) <- CountTable$Geneid
CountTable <- CountTable[,-1]
CountTable[is.na(CountTable)] <- 0 
#check PRAME expression levels
counttable <- CountTable
counttable$geneid <- mapIds(org.Hs.eg.db, keys=row.names(counttable), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
PRAME <- counttable[which(counttable$geneid == "PRAME"),]
PRAME <- PRAME[,-7]
PRAME <- rbind(colnames(PRAME),PRAME)
PRAME <- t(PRAME)
colnames(PRAME) <- c("sample","counts")
PRAME <- as.data.frame(PRAME)
PRAME$counts <- as.integer(PRAME$counts)
p <- ggplot(data=PRAME,aes(x=sample,y=counts)) + geom_bar(stat="identity",fill="steelblue") +
  ggtitle("PRAME expression levels among samples")
p

CountTable <- CountTable[,c(-1,-4)]
samples <- data.frame(sample = c("lo_ctrl2","lo_ctrl3","lo_oe2","lo_oe3"),
                      batch = c("2","3","2","3"),
                      condition = c("Ctrl","Ctrl","oe","oe"))
#samples <- data.frame(sample = c("lo_ctrl1","lo_ctrl2","lo_ctrl3","lo_oe1","lo_oe2","lo_oe3"),
#                      batch = c("1","2","3","1","2","3"),
#                      condition = c("Ctrl","Ctrl","Ctrl","oe","oe","oe"))

samples$batch <- factor(samples$batch,levels=c("1","2","3"))
samples$condition <- factor(samples$condition,levels=c("Ctrl","oe"))

dds <- DESeqDataSetFromMatrix(countData = CountTable, colData=samples, design=~batch+condition)
#dds <- DESeqDataSetFromMatrix(countData = CountTable, colData=samples, design=~condition)
dds = DESeq(dds)
counts(dds)->raw_counts
raw_counts<-as.data.frame(counts(dds))

#Create a normalized matrix
norm_counts = counts(dds, normalized = TRUE)

### PRELIMINARY ANALYSES ###
# The first steps in your analysis should focus on better understanding the relationship of the datasets being studied. This can
# be simply achieved by generating a PCA plot showing the relationship of your samples.
# First we transform our raw count data using a variance stabilizing transformation (VST) that roughly mirrors how DeSeq2 models the data.
vsd1 <- varianceStabilizingTransformation(dds, blind=FALSE)

# Then we plot a PCA, grouping and coloring our datasets according to batch
plotPCA(vsd1, c("condition","batch"))
data_vsd1 <- plotPCA(vsd1, "condition",returnData=TRUE)

# we can also attempt to replicate the batch effect correction performed by DeSeq2 using the limma::removeBatchEffect function
vsd2 <- varianceStabilizingTransformation(dds, blind=FALSE)
assay(vsd2) <- limma::removeBatchEffect(assay(vsd2), vsd2$batch)
plotPCA(vsd2, c("condition","batch"))
data_vsd2 <- plotPCA(vsd2, "condition",returnData=TRUE)


percentVar <- round(100 * attr(data_vsd1, "percentVar"))
pca1<- ggplot(data_vsd1, aes(PC1, PC2, color=condition,label=name)) +
  theme_bw() + 
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  geom_point(size=5)+
  geom_label_repel(show.legend = FALSE)+
  scale_color_manual(values=c("#E8681D","#41B7C4"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(text = element_text(size=15))+
  theme(legend.position="bottom") +
  theme(axis.text = element_text( color = "black", size = 15))+
  theme(legend.title=element_blank())+
  scale_y_continuous(expand = expansion(add = 2))+
  ggtitle("PCA before batch correction")

percentVar <- round(100 * attr(data_vsd2, "percentVar"))
pca2<-ggplot(data_vsd2, aes(PC1, PC2, color=condition,label=name)) +
  theme_bw() + 
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  geom_point(size=5)+
  geom_label_repel(show.legend = FALSE)+
  scale_color_manual(values=c("#E8681D","#41B7C4","purple"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(text = element_text(size=15))+
  theme(legend.position="bottom") +
  theme(axis.text = element_text( color = "black", size = 15))+
  theme(legend.title=element_blank())+
  scale_y_continuous(expand = expansion(add = 2))+
  ggtitle("PCA after batch correction")

pca1/pca2


#if you have multiple comparisons you can put them in line 120 
comparisons = rbind(c("Ctrl","oe"))

adj_p_val = 0.1
abs_log2fc = 0


for(i in 1:nrow(comparisons)){
  res = results(dds, contrast = c("condition", comparisons[i,2], comparisons[i,1]))
  res_sig = res[res[,6] < adj_p_val & !is.na(res[,6]) & abs(res[,2]) >abs_log2fc,]
  assign(paste(comparisons[i,2],"_", comparisons[i,1], sep =""), res)
  assign(paste(comparisons[i,2],"_", comparisons[i,1],"_sig", sep =""), res_sig)
}


#dim(`Ctrl_oe_sig`)


Ctrl_oe_df<-as.data.frame(counts(dds))

# Here we modify our output data to include two additional columns that contain the baseMeans (a proxy for counts)
# This is useful for downstream filtering of lowly expressed genes

baseMeanCtrl = rowMeans(counts(dds,normalized=TRUE)[,dds$condition == "Ctrl"])
baseMeanoe = rowMeans(counts(dds,normalized=TRUE)[,dds$condition == "oe"])


res1 = cbind(as.data.frame(res), baseMeanCtrl, baseMeanoe)
res1$padj <- p.adjust(res1$pvalue,method= "BH")

# Here we add two further columns, the gene symbol (common name) and entrez ID - both of which may be useful downstream
res1$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res1), column="SYMBOL", keytype="ENSEMBL", multiVals="first") # MAPS GENE IDs
res1$entrez <- mapIds(org.Hs.eg.db, keys=row.names(res1), column="ENTREZID", keytype="ENSEMBL", multiVals="first")

# Finally we write the complete results object to an outfile
write.csv(res1, "/Users/staceyvann/Desktop/2022Spring/AppliedSequencing/Final/DEG_lo.csv", row.names=TRUE)

sigRNA<-as.data.frame(res1)
sigRNA$log10.pvalue<-(-1*log10(sigRNA$padj))
sigRNA$Significant.Gene<-"No"

row2<-which(sigRNA$padj< adj_p_val & !is.na(sigRNA$padj) & abs(sigRNA$log2FoldChange)>abs_log2fc)
sigRNA$Significant.Gene[row2]<-"Yes"
sigRNA$Differential_Gene<-sigRNA$Significant.Gene
sigRNA$Differential_Gene[which(sigRNA$padj<adj_p_val & sigRNA$log2FoldChange>abs_log2fc)]<-"Up-Regulated"
sigRNA$Differential_Gene[which(sigRNA$padj<adj_p_val & sigRNA$log2FoldChange< -abs_log2fc)]<-"Down-Regulated"


#RNA volcano

#make sure to adjust xlim and ylim based off of the data. Run hist(sigRNA$log2FoldChange)  or range to get an idea

hist(sigRNA$log10.pvalue) 
hist(sigRNA$log2FoldChange) 

# pick a set of genes to label:
sigRNA$label<-"No"
sigRNA$label[which(sigRNA$padj< adj_p_val & abs(sigRNA$log2FoldChange)>abs_log2fc)]<-"Yes"


#generate volcano plot
#pdf("Volcano.pdf", width=8, height=4,useDingbats=FALSE)
pRNA <- ggplot(sigRNA, aes(log2FoldChange, log10.pvalue))+
  geom_point(aes(colour = Differential_Gene),size=1)+
  scale_colour_manual(values=c("blue", "black","red"))+
  xlim(-10,10)+
  ylim(0,20)+
  xlab(expression(paste(log[2], 'FC(Ctrl/oe)')))+
  ylab(expression(paste(-log[10], '(Q-value)'))) +
  ggtitle("Differentially expressed genes in PRAMElo over-expression experiment")+
  theme_bw() + 
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(text = element_text(size=10)) +
  theme(axis.text = element_text( color = "black", size = 10))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line=element_line()) +
  theme(legend.position="none")+
  geom_vline(xintercept=-abs_log2fc,linetype="dotted")+
  geom_vline(xintercept=abs_log2fc,linetype="dotted")+
  geom_text_repel(
    data = subset(sigRNA, label=="Yes"),
    aes(label = sigRNA$symbol[which(sigRNA$label=="Yes")]),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  )
pRNA

#

Significant_table<-sigRNA[which(sigRNA$Significant.Gene=="Yes"),]


#Up/Down regulated genes
up_reg <- filter(sigRNA, Differential_Gene == "Up-Regulated")
down_reg <- filter(sigRNA, Differential_Gene == "Down-Regulated")
dim(up_reg)
dim(down_reg)
DEG_sig <- rbind(up_reg, down_reg)
write.csv(DEG_sig, "/Users/staceyvann/Desktop/2022Spring/AppliedSequencing/Final/DEGdeseq2_lo.csv", row.names=TRUE)

#GO and KEGG analysis
#GO and KEGG analysis 

go <- enrichGO(
  gene = up_reg$entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = T
)
kegg <- enrichKEGG(
  gene = up_reg$entrez,
  keyType = "kegg",
  pvalueCutoff = 0.05
)
dotplot(go, title = "GO analysis of upregulated genes in PRAME over-expression")

dotplot(kegg, title = "KEGG analysis of upregulated genes PRAME over-expression")

go <- enrichGO(
  gene = down_reg$entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = T
)
kegg <- enrichKEGG(
  gene = down_reg$entrez,
  keyType = "kegg",
  pvalueCutoff = 0.05
)
dotplot(go, title = "GO analysis of down-regulated genes PRAME over-expression")

dotplot(kegg, title = "KEGG analysis of down-regulated genes PRAME over-expression")
