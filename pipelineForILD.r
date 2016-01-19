library(CLL)
setwd("/home/owht/桌面/KIZ/shuju from web/ILD/gse21411/")
rawData<-read.table(filenames=list.celfiles())
library(gcrma)
rawDataGcrma<-gcrma(rawData)
library(graph)
eset<-exprs(rawDataGcrma)
pearson_cor<-cor(eset)
dist.lower<-as.dist(1-pearson_cor)
hc<-hclust(dist.lower,"ave")
plot(hc,main="The Cluster Dendrogram of GSE21411's Gene Expression Data")

#then 86,88,94 delete!!!!

colNamesForGSE21369<-c(533882:533885,533887,533889:533893,533895:533910)
colNamesForGSE21369<-paste(colNamesForGSE21369,"CEL",sep=".")
GSE21369ForPlot2000genes<-GSE21369[1:2000,]
pheatmap(GSE21369ForPlot2000genes,main = "Expression Hearmap of 2000 Genes in GSE21369 Selected Samples",legend_labels = "Scaled log2 intensity")
