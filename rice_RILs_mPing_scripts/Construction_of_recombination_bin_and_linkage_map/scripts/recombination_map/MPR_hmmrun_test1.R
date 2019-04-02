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


if (0){
##################################################################################
#step.2 refine MPR
## then you need to refine the MPR results
set.seed(123);system.time(all.res <- globalMPRRefine(myBaseData,alleleA=na.omit(
             allele.MPR[,1]),numGroup=47,groupSort=TRUE,numPerm=20,numTry=3,
             numBaseStep=50,numBaseCandidateStep=100,numKnownStep=30,
             numKnownCandidateStep=50,useMedianToFindKnown=TRUE,maxIterate=150,
             maxNStep=3,scoreMin=0.8,saveMidData=TRUE,verbose=TRUE))

## summarize results using Bayesian inference
perm <- 10
res <- all.res$midData[[perm]]
table(geno.res <- genotypeCallsBayes(res$call,errorRate=1e-11,eps=1e-10,maxIterate=100,
             verbose=FALSE)$type)

## reconstruct parental genotypes
allele.MPR <- res$allele
allele.MPR[geno.res==2|rowMin(res$call)>=perm,] <- NA
allele.MPR[geno.res==3,] <- allele.MPR[geno.res==3,c(2,1)]

## predicted alleles of parent 1
table(is.na(alleleA <- allele.MPR[,1]))
}

#####################################################################################
#step.3 genotype RILs by HMM
snpSet <- sort(rownames(na.omit(allele.MPR)))
ids <- match(snpSet,rownames(myBaseData))
table(is.na(geno.data <- base2Geno(myBaseData[ids,],allele.MPR[ids,])))

geno.data.cr <- geno.data
SNPbyChr <- split(1:nrow(geno.data),substr(rownames(geno.data),1,2))
for(chr in names(SNPbyChr)){
  cat("\r",chr)
  ids <- SNPbyChr[[chr]]
  geno.data.chr <- correctGeno(geno.data[ids,],correct.FUN=correctFUNHMM,hmmFUN=hmm.vitFUN.rils,
                               geno.probability=c(0.4975, 0.4975,0.005),transitionFUN=phy2get.haldane.rils,
                               emissionFUN=makeEmissionFUN(errorRate=0.0106))
  geno.data.cr[ids,] <- geno.data.chr
}


##write basedata, genotype and HMM corrected genotype data into file
write.table(allele.MPR,"MPR.allele.MPR",sep="\t",col.names=NA)
write.table(myBaseData,"MPR.base.data",sep="\t",col.names=NA)
write.table(geno.data,"MPR.geno.data",sep="\t",col.names=NA)
write.table(geno.data.cr,"MPR.geno.data.HMMcr",sep="\t",col.names=NA)
