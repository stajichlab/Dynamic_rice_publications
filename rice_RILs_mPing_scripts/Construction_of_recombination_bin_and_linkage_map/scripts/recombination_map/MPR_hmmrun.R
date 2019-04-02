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
                               emissionFUN=makeEmissionFUN(errorRate=0.014))
  geno.data.cr[ids,] <- geno.data.chr
}


##write basedata, genotype and HMM corrected genotype data into file
write.table(allele.MPR,"MPR.allele.MPR",sep="\t",col.names=NA)
write.table(myBaseData,"MPR.base.data",sep="\t",col.names=NA)
write.table(geno.data,"MPR.geno.data",sep="\t",col.names=NA)
write.table(geno.data.cr,"MPR.geno.data.HMMcr",sep="\t",col.names=NA)

#####################################################################################
#step.4 genotype to bin map
str(geno.data.subset <- geno.data.cr[rowSums(!is.na(geno.data.cr))>0,])
snpSite <- as.numeric(substr(rownames(geno.data.subset),1,10))
geno.data.bin <- genoToBin(geno.data.subset,base.position=snpSite,corrected=TRUE,size=250e3,num=5,fillSmallNA=TRUE,minBinsize=5e3,heterozygote=FALSE)
## genptype bins with missing data
str(geno.bin <- geno.data.bin[[2]])
write.table(geno.bin,"MPR.geno.bin",sep="\t",col.names=NA)

#################################################################################
#step.5 fill missing data in bin map
library(qtl)
#/rhome/cjinfeng/software/tools/MPR/MPR_data/phenoData.txt
#phenoData<-as.matrix(read.table("/rhome/cjinfeng/software/tools/MPR/MPR_data/phenoData.txt",header=T,row.name=1))
phenoData<-as.matrix(read.table("/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.first",header=T,sep="\t",row.name=1))

## Attention!!!!! the name of SNP MUST BE started with ten numeric characters, the first two characters for chromosome and the other eight charactes for physical position of the SNP, otherwise you have to modify the code of geno2Cross (the row assigning value to myGeno.site)
summary(myCrossData <- geno2Cross(geno.bin,phenoData))
cross.data <- myCrossData

myMap <- est.map(cross.data, error.prob=0.00, map.function=c("haldane","kosambi","c-f","morgan")[1], m=0, p=0, maxit=1000, tol=1e-6, sex.sp=FALSE, verbose=FALSE, omit.noninformative=TRUE)
summary(myMap)
plot.map(myMap)

cross.data <- replace.map(cross.data, myMap)
cross.data <- calc.genoprob(cross.data,error.prob=0.00)
cross.fill <- fill.geno(cross.data,method=c("imp","argmax")[2], error.prob=0.00,map.function="haldane");
cross.fill <- calc.genoprob(cross.fill,error.prob=0.00)

id.pheno <- 1
## res <- cim(cross.data,pheno.col=id.pheno,n.marcovar=3,window=10,method="em",imp.method=c("imp", "argmax")[2], error.prob=0.01)
res <- scanone(cross.fill,pheno.col=id.pheno,method="em")
summary(res,threshold=4, lodcolumn=1)
i <- which.max(res$lod)
res[(i-3):(i+3),]
#plot(res,chr=5) ## vtype="S",chr=5,xlim=c(25,30)
plot(res)

str(geno.fill <- t(pull.geno(cross.fill))-1)
colnames(geno.fill) <- colnames(geno.bin)

## Final Bin Map
str(geno.uniq <- mergeBinMap(geno.fill))
## uniq bin map for R/qtl
cross.uniq <- geno2Cross(geno.uniq,phenoData)
uniq.Map <- est.map(cross.uniq, error.prob=0.00, map.function=c("haldane","kosambi","c-f","morgan")[1], m=0, p=0, maxit=1000, tol=1e-6, sex.sp=FALSE, verbose=FALSE, omit.noninformative=TRUE)
cross.uniq <- replace.map(cross.uniq, uniq.Map)
cross.uniq <- calc.genoprob(cross.uniq,error.prob=0.00)

summary(uniq.Map)
plot.map(uniq.Map)

##write genotype bin, fill bin and combined bin data into files
write.cross(cross.data,"csv","MPR.cross")
write.cross(cross.fill,"qtlcart","MPR.cross.fill")
write.cross(cross.uniq,"qtlcart","MPR.cross.uniq")
#write.table(geno.bin,"MPR.geno.bin",sep="\t",col.names=NA)
write.table(geno.fill,"MPR.geno.bin.fill",sep="\t",col.names=NA)
write.table(geno.uniq,"MPR.geno.bin.uniq",sep="\t",col.names=NA)

#################################################################################
#step.6 comparision of bin map between orginal and filled map
if(0){
par(mfrow=c(2,1),mar=c(4.1,4.1,2.1,1.1))
sel.cols <- 1:10

## comparisons
## raw
geno.colors <- geno.data[,sel.cols]+1;geno.colors[is.na(geno.colors)] <- 0
plot(matrix(1:nrow(geno.data),nrow=nrow(geno.data),ncol=length(sel.cols)),
     matrix(rep(sel.cols,each=nrow(geno.data)),nrow=nrow(geno.data),ncol=length(sel.cols)),col=geno.colors,pch="|",
     xlab='',ylab="RIL index",mar=c(1.1,2.1,2.1,1.1))
## corrected
geno.colors <- geno.data.cr[,sel.cols]+1;geno.colors[is.na(geno.colors)] <- 0
plot(matrix(1:nrow(geno.data),nrow=nrow(geno.data),ncol=length(sel.cols)),
     matrix(rep(sel.cols,each=nrow(geno.data)),nrow=nrow(geno.data),ncol=length(sel.cols)),col=geno.colors,pch="|",
     xlab="SNP index",ylab="RIL index",mar=c(4.1,2.1,2.1,1.1))
}
dev.off()

sessionInfo()
## R version 2.5.0 (2007-04-23) 
## i686-pc-linux-gnu 

## locale:
## LC_CTYPE=en_US.UTF-8;LC_NUMERIC=C;LC_TIME=en_US.UTF-8;LC_COLLATE=en_US.UTF-8;LC_MONETARY=en_US.UTF-8;LC_MESSAGES
## =en_US.UTF-8;LC_PAPER=en_US.UTF-8;LC_NAME=C;LC_ADDRESS=C;LC_TELEPHONE=C;LC_MEASUREMENT=en_US.UTF-8;LC_IDENTIFICATION=C

## attached base packages:
## [1] "splines"   "tools"     "stats"     "graphics"  "grDevices" "utils"    
## [7] "datasets"  "methods"   "base"     

## other attached packages:
##        MPR genefilter   survival    Biobase 
##      "0.1"   "1.14.1"     "2.34"   "1.14.1"

