library(Biobase)
library(MPR)
pdf("NB.RILs.pdf")
#################################################################################
#step.0 load data
## load SNP alleles at putative SNP sites
#data(snpData)
#snpData <- as.matrix(read.table("/rhome/cjinfeng/software/tools/MPR/MPR_data/snpData_chr05.txt"))
snpData <- as.matrix(read.table("NB.RILs.dbSNP.SNPs.RILs"))

## load overlapping SNP alleles from low-depth sequences of one parent
#data(markerData)
#markerData <- as.matrix(read.table("/rhome/cjinfeng/software/tools/MPR/MPR_data/markerData_chr05.txt",row.names=1,header=T))[,1]
markerData <- as.matrix(read.table("NB.RILs.dbSNP.SNPs.Markers",row.names=1,header=T))[,1]

#################################################################################
#step.1 MPR
set.seed(123)

if(0){ ## set to 0 if want to make a whole-genome prediction and set to 1 if want to test program
## select 30 markers randomly
markers <- sample(names(markerData)[10:50],20)

## select SNP sites which contain the 30 markers
ids <- match(markers,rownames(snpData))
str(myBaseData <- snpData[min(ids):max(ids),])
}else{
myBaseData <- snpData #use all SNP data not a subsample
}

if(0){
## global MPR aiding with marker data 
allele.MPR <- globalMPRByMarkers(myBaseData,markers=markerData,numTry=3,numBaseStep=50,
             numBaseCandidateStep=100,numMarkerStep=10,useMedianToFindKnown=TRUE,
             maxIterate=150,maxNStep=3,scoreMin=0.8,verbose=TRUE)

#allele.MPR <- allele.MPR[-which(apply(allele.MPR,1,function(x)all(is.na(x)))),]
}else{
allele.MPR <- as.matrix(read.table("NB.RILs.dbSNP.SNPs.sub.parents"))
}
