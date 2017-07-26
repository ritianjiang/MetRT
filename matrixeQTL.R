#test for matrix eQTL
library(reshape)
library(TCGAbiolinks)
GDCquery_Maf(tumor = "KIRC",save.csv = T, directory = "GDCdata")
setwd("/home/owht/下载/data_from_web/human/KIRC/gdac.broadinstitute.org_KIRC.Mutation_Packager_Raw_Calls.Level_3.2016012800.0.0/GDCdata/")
KIRCSomatic<-read.csv("TCGA.KIRC.somatic.csv",header = T,stringsAsFactors = F)
KIRCSomatic$Count<-rep(0,nrow(KIRCSomatic))
KIRCSomatic[KIRCSomatic$allele1 != KIRCSomatic$allele2,]$Count<-1
KIRCMM<-cast(KIRCSomatic,location~TCGA.Barcode);KIRCMM[is.na(KIRCMM)]<-0
colN<-colnames(KIRCMM)
colN<-data.frame(colN)
colN$colN<-as.character(colN$colN)
colN$V2<-rep(0,nrow(colN))
for (i in 2:nrow(colN)){
  colN[i,]$V2<-(strsplit(x = colN[i,]$colN,split = "A-0"))[[1]][1]}
colN[1,]$V2<-"snpid"

colnames(KIRCMM)<-c(colN$V2)
write.table(KIRCMM,file = "/home/owht/下载/data_from_web/temp/TEST/KIRCMM.csv",sep = "\t",row.names = F)

KIRCMM<-read.table( "/home/owht/下载/data_from_web/temp/TEST/KIRCMM.csv",header = T)
KIRCEx<-KIRCgenomicMatrix[,colnames(KIRCgenomicMatrix) %in% colnames(KIRCMM)]
KIRCEx[,1]<-rownames(KIRCEx);colnames(KIRCEx)[1]<-"geneid";
write.table(KIRCEx,file = "/home/owht/下载/data_from_web/temp/TEST/KIRCEx.csv",row.names = F,sep = "\t")

KIRCsnpid<-KIRCMM$snpid
KIRCMM1<-KIRCMM[,colnames(KIRCMM) %in% colnames(KIRCEx)]
write.table(KIRCMM1,file = "/home/owht/下载/data_from_web/temp/TEST/KIRCMM1.csv",sep = "\t",row.names = F)

KIRCMM1<-KIRCMM[,colnames(KIRCEx)[-1]]
KIRCMM1<-data.frame(snpid = KIRCsnpid,KIRCMM1)
write.table(KIRCMM1,file = "/home/owht/下载/data_from_web/temp/TEST/KIRCMM1.csv",sep = "\t",row.names = F)

library(MatrixEQTL)

## Location of the package with the data files.
base.dir = "/home/owht/下载/data_from_web/temp/TEST/";
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste(base.dir, "KIRCMM1.csv", sep="");

# Gene expression file name
expression_file_name = paste(base.dir, "KIRCEx.csv", sep="");

output_file_name = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold = 1e-2;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)

## Plot the histogram of all p-values

plot(me)
