library(nlme)
library(pbapply)
library(parallel)
library(magrittr)
library(data.table)
load("IQR10_meth135.RData")
load("total135_IQR10_meth_anno.RData")
load("FPKM135_ens.RData")
source("methExp_each.R")
sampleinfo<-read.table(file ="meth_135_sampleinfo",header = T)

total135_meth_anno$loc<-paste(total135_meth_anno$seqnames,
                              total135_meth_anno$start,
                              total135_meth_anno$end,sep = ".")
meth_in_anno<-meth_135IQR10[rownames(meth_135IQR10) %in%total135_meth_anno$loc,]

##precorrect the sampleinfo
rownames(sampleinfo)<-sampleinfo$sub
sampleinfo<-sampleinfo[colnames(meth_135IQR10),]
sampleinfo$sub<-as.factor(sampleinfo$sub)

###pbapply
cl<-makeForkCluster(4)
a<-pblapply(1:nrow(total135_meth_anno),FUN = methExp_each,
            relation = total135_meth_anno,
            meth = meth_in_anno,
            exp = fpkm_135_total_ens,
            sampleinfo = sampleinfo,
            cl = cl)
stopCluster(cl)
a<-do.call(rbind,a)
save(a,file = "p_slope.RData")











