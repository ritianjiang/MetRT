#This script can process the GTF get from BLAT; Multi-scaffolds
#alignment will be trimmed and transcript will be produced;
#Wang haotian
#20180122

library(stringr)
setwd("/home/owht/KIZ/shuju_from_web/Bowhead_whale/BowheadWhale_gtf")
gtf<-read.table("PCG_bowheadwhalr_11_intron.gtf",
                col.names = c("chr","tool","site","start","end",
                           "score","strand","zero","gene"),
                stringsAsFactors = F)
gtf$site<-"CDS"
gtf$gene<-str_split(gtf$gene,pattern = ";",n=2,simplify = T)[,1]

#we get thegtf from blat; its start is 0;
gtf$start<-gtf$start+1 #Caution please!!
#if the gtf is not from blat, delete the line above.

gtf_tras<-data.frame(chr=0,tool="BLAT",site="transcript",start=0,
                     end=0,score=100,strand = "f",zero=".",
                     attribute=unique(gtf$gene),stringsAsFactors = F)

TrimGTF<-file("PCG_bowheadwhalr_11_intron.gtf.trimmed.plustrans","w")
for(i in 1:nrow(gtf_tras)){
  tempdf<-gtf[gtf$gene == gtf_tras[i,]$attribute,]
  if(length(unique(tempdf$chr))>1){
      stdsca<-table(tempdf$chr) %>% data.frame()
      stdsca<-stdsca[stdsca$Freq==max(stdsca$Freq),]$Var1[1] %>% as.character()
      tempdf<-tempdf[tempdf$chr == stdsca,]
  }
  tempdf$gene<-paste("gene_id",
                     paste0("\"",tempdf[1,]$gene,"\";"),
                     "transcript_id",
                     paste0("\"",tempdf[1,]$gene,".1\";"))
  write.table(tempdf,file=TrimGTF,quote = F,col.names = F,
              row.names = F,sep="\t",append = T)
  gtf_tras[i,]$chr<-tempdf[1,]$chr
  gtf_tras[i,]$start<-tempdf[1,]$start
  gtf_tras[i,]$end<-tempdf[nrow(tempdf),]$end
  gtf_tras[i,]$strand<-tempdf[1,]$strand
  gtf_tras[i,]$attribute<-tempdf[1,]$gene
  write.table(gtf_tras[i,],file=TrimGTF,quote = F,col.names = F,
              row.names = F,sep="\t",append = T)
}
close(TrimGTF)
#write.table(gtf_tras,file="merged.trans",quote=F,col.names = F,row.names = F,sep="\t")

# setwd("/home/owht/KIZ/shuju_from_web/Bowhead_whale/BowheadWhale_gtf")
# pcg_mrna<-read.table("PCG_BowheadW_withmrna.gtf",stringsAsFactors = F)
# 
# exon=0
# for(i in 1:nrow(pcg_mrna)){
#   if(pcg_mrna[i,]$V3 == "exon"){
#     exon=exon+1;
#     id<-str_split_fixed(pcg_mrna[i,]$V9,pattern = ";",n=2)[1]
#     pcg_mrna[i,]$V9<-paste("gene_id",
#                            paste0("\"",id,"\"",";"),
#                            "transcript_id",
#                            paste0("\"",id,"_trans","\";"),
#                            "exon_number",
#                            paste0("\"",exon,"\"",";"))
#   }
#   else{
#     exon=0;
#     pcg_mrna[i,]$V9<-paste("gene_id",
#                            paste0("\"",id,"\"",";"),
#                            "transcript_id",
#                            paste0("\"",id,"_trans",";"))
#   }
# }
# write.table(pcg_mrna,file="PCG_Bowhead_withmrna2.gtf",col.names = F,row.names = F,quote=F,sep='\t')