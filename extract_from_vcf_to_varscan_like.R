########################################################
#a simple script to extract the RD AD infomation
#RD is the reads that support the ref, AD is for the alt
########################################################
library(stringr)
setwd("The_PATH_to_vcf_files")
library(VariantAnnotation)
filelist<-list.files()
for (i in 1:length(filelist)){
  tempName<-filelist[i]
  aliquote<-str_split_fixed(tempName,"\\.",3)[2]
  vcf<-readVcf(tempName)
  rr<-rowRanges(vcf)
  location<-rr@ranges@start
  ref<-rr@elementMetadata@listData$REF@metadata
  AD<-geno(vcf)$AD
  RD<-geno(vcf)$RD
  QUAL<-rr@elementMetadata$FILTER
  totoal_info<-data.frame(AD,RD,location,QUAL)
  colnames(totoal_info)<-c("AD","RD","loc","FILTER")
  aaa<-str_split_fixed(rownames(totoal_info),pattern = "[/_]",3)
  totoal_info$ref<-aaa[,2]
  totoal_info$alt<-aaa[,3]
  newPath<-str_c("./csv/",aliquote)
  write.csv(x=totoal_info,file = newPath)
}
