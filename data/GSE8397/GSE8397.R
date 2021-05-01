GSE8397

library(limma)
library(R.utils)
library(affy)
library(tidyr)
library(hgu133acdf)
library(RColorBrewer)
library(hgu133a.db)
library(calibrate)
library(GEOquery)



#############
#download the matrix file from GEO and only keep the case, age and gender columns
#############
library(GEOquery)
gse8397 <- getGEO("GSE8397",GSEMatrix=TRUE)  #gets the matrix file with information about samples from the GEO webiste
pData <- pData(phenoData(gse8397[[1]]))[,c(1,10)]
pData <- pData %>%
	separate(characteristics_ch1, c("age", "sex"), "; gender: ")
pData <- pData[ grep("Medial", pData$title),]
####delete these two are are both potential for Alz
pData <- pData[ -grep("MS155", pData$title),]
pData <- pData[ -grep("PDC1", pData$title),]
pData$title <- gsub("[[:digit:]]+", "", pData$title)
pData$title <- gsub("Medial substantia nigra Parkinson's disease case  - A chip", "PD", pData$title) 
pData$title <- gsub("Medial substantia nigra Parkinson's disease case  -A chip", "PD", pData$title) 
pData$title <- gsub("Medial substantia nigra control case  - A chip", "control", pData$title) 
pData$title <- gsub("Medial substantia nigra control case   - A chip", "control", pData$title) 
pData$age <- gsub("age: ", "", pData$age) 
colnames(pData)[1] <- "case"  #change coloumn names


##########
#switch case and control so case info is first (required for metaUnion package)
##########
GSE1 = pData[7:21,] 
GSE2 = pData[1:6,]
pData <- rbind(GSE1, GSE2)

###########
#read in the cel files into raw data 
###########
#untar("GSE8397_RAW.tar")
#cels = list.files(pattern = "cel")
#sapply(paste(cels, sep="/"), gunzip)
cels = list.files(pattern = "cel")
i <- cels[7:21]
p <- cels[1:6]
cels <- c(i, p)
raw.data=ReadAffy(filenames = cels, phenoData=pData, cdfname="hgu133acdf")


############
#plot the prenomalised plots
############

###
#plot boxplots of each sample
###
library(RColorBrewer)
brewer.cols <- brewer.pal (6, "Set1")
png("Unprocessed log scale probe intensities.png")
boxplot(raw.data, col = brewer.cols, ylab = "Unprocessed log (base 2)scale Probe Intensities", xlab = "Array Names", main = "Unprocessed log(base2) scale Probe Intensities of each Array")
dev.off()

#Construct density plots
png("Unprocessed density plot.png")
hist (raw.data, main = "Density plot of log(2) probe intensities", col = brewer.cols, lty = 1, xlab = "Log (base 2) Intensities", lwd = 3)
samp.leg.names <- colnames(raw.data)
legend (12,0.6, cex = 0.5, legend = samp.leg.names, lty = 1,col = brewer.cols, lwd = 3)
dev.off()

#plot the MA plots for each sample and saved in one pdf
pdf("MAplots.pdf")
MAplot(raw.data)
dev.off()
#-created 29 plots off all of the arays - -Quality problems are most apparent from an MA-plot in cases where the loess smoother oscillates wildly or if the variability of the M values appears to be greater in one or more arrays relative to the others-  Bolstad BM,  Collin F,  Brettschneider J, et al. Gentleman R,  Carey VJ,  Huber W,  Irizarry RA,  Dudoit S. Quality assessment of Affymetrix GeneChip Data, Bioinformatics and Computational Biology Solutions Using R and Bioconductor , 2005New YorkSpringerpg. 36 .



############
#normalise data using RMA
############
data.rma.norm = rma(raw.data)
data.rma.exprs = exprs(data.rma.norm)
write.table(data.rma.exprs, file = "GSE8397_rma_data.txt", quote = F)

############
#plot normalised data
############
png("RMA normalised log scale probe intensities.png")
boxplot(data.rma.exprs, col = brewer.cols, ylab = "RMA normalised log (base 2)scale Probe Intensities", xlab = "Array Names", main = "RMA normalised log(base2) scale Probe Intensities of each Array")
dev.off()
nrow(data.rma.exprs) # 22283


################
#make MAS5 present absent calls and exclude probes with over 15% absent
###############
calls.mas5 <- mas5calls(raw.data)
calls.exprs <- exprs(calls.mas5)
calls.sum <- rowSums(calls.exprs == "A")    #returns a 1 in matrix if probe is absent
calls.sum = calls.sum[rownames(data.rma.exprs)]    #filter out the probes that have been filtered out before
present.call <- calls.sum[calls.sum < (0.15*nrow(pData))] #sums up the ones in a row and has to be below 10% of samples
data.rma.norm1 = data.rma.norm[names(present.call),]
data.rma.exprs1 = data.rma.exprs[names(present.call),]
nrow(data.rma.exprs1) #7207

############
#filtering by mean - with the bottom 5% being removed
############
myVar <- exprs(data.rma.norm1)
myVar <- rowMeans(myVar)
myVar1 <- myVar[myVar > quantile(myVar, 0.05) ]  
data.rma.norm2 = data.rma.norm1[names(myVar1),]
data.rma.exprs2 = data.rma.exprs1[names(myVar1),]
nrow(data.rma.exprs2) #6846



#############
#filter out any without or duplicate genes and sort those with muliple probes
#############
probes=row.names(data.rma.exprs2)     #get a list of probes that are in the raw.data
Symbols = unlist(mget(probes, hgu133aSYMBOL, ifnotfound=NA))   #get a list of symbols that match the probes of raw.data
Entrez_IDs = unlist(mget(probes, hgu133aENTREZID, ifnotfound=NA))   #get a list of entrez IDs that match the probes of raw.data
data.rma.exprs.annotate =cbind(probes,Symbols,Entrez_IDs,data.rma.exprs2) #binds probes, entrez Ids, symbols and expression set together


remove = function(x) {
	comma = grep(",", rownames(x))    #find any row with comma in them
        if(length(comma)>0){
        	x = x[-comma,]
        }
	space = grep(" ", rownames(x))    #find any row with space in them
        if(length(space)>0){
        	x = x[-space,]
        }
	semicolon = grep(";", rownames(x))    #find any row with semicolon in them
        if(length(semicolon)>0){
        	x = x[-semicolon,]
	}
	colon = grep(":", rownames(x))    #find any row with colon in them
        if(length(colon)>0){
        	x = x[-colon,]
        }
	i <- is.na(x[,3])
	l <- x[!i,]      # remove any probes that don't have a mapped gene
	id=grep("AFFX",rownames(l))
	m = l[-id,]      #remove any probes that are controlr AFFX probes
}



data.rma.exprs.sorted = remove(data.rma.exprs.annotate)

data.rma.exprs.varmed = data.rma.exprs2[rownames(data.rma.exprs.sorted),]
data.rma.norm.varmed = data.rma.norm2[rownames(data.rma.exprs.sorted),]
nrow(data.rma.exprs.varmed) #6329



################
#create a coloumn of probe expression average for control and case 
################
##--averaging the control and case samples to make coloumns (will be used in MA plot later)
control.row <- which(pData$case == "control")
control.row = rownames(pData[control.row,])
control.row <- data.rma.exprs.varmed[,control.row]
control.row <- rowSums(control.row)/ncol(control.row)


case.row <- which(pData$case == "PD")
case.row = rownames(pData[case.row,])
case.row <- data.rma.exprs.varmed[,case.row]
case.row <- rowSums(case.row)/ncol(case.row)

###############
#plot RMA of control and against the PD
###############
A = (case.row + control.row) / 2
M = case.row - control.row
plot(A, M)
png("MA plot of PD and Control.png")
plot(A, M, main="MA plot of PD vs Control", pch=19, cex=0.2)
dev.off()




###############
#create design factor
###############
samples <- data.rma.norm.varmed$case   
samples <- as.factor(samples)
gender <- data.rma.norm.varmed$sex
gender <- as.factor(gender)  #factor of if male or not
age <- data.rma.norm.varmed$age
age <- as.numeric(age) 
design <- model.matrix(~0 + samples+ gender +age)  #design is a model matrix of intercept, case or PD, gender and age
colnames(design)[1] <- "Control"    #renames the first column to Control
colnames(design)[2] <- "PD"    #renames the second column to PD



################
#limma to create DEGs
#################
library(limma)

limma = function(x, design) {
	fit <- lmFit(x ,design)
	contrast.matrix <- makeContrasts(PD - Control, levels=design)
	fit2 <- contrasts.fit(fit, contrast.matrix)
	fit2 <- eBayes(fit2)
	limma.df = fit2$df.total
	x <- topTable(fit2, coef=1, adjust="BH",number=nrow(x),sort.by="none") #coef 1 is used as that is the comparision i want(also only one coef)
	cbind(x,limma.df)
}


var.medlimma = limma(data.rma.exprs.varmed, design)
dd <- var.medlimma[(var.medlimma$adj.P.Val <0.05),]
###(1869 DEG under 0.05 adj p val in here)
probes.Dif.Exprs = row.names(var.medlimma)
Symbol = unlist(mget(probes.Dif.Exprs, hgu133aSYMBOL, ifnotfound=NA))
Entrez.Gene = unlist(mget(probes.Dif.Exprs, hgu133aENTREZID, ifnotfound=NA))
var.medlimma = data.frame(var.medlimma,Entrez.Gene, Symbol)
########





################
#annotate the files
################
probes.limma = row.names(data.rma.exprs.varmed)
Symbol = unlist(mget(probes.limma, hgu133aSYMBOL, ifnotfound=NA))
Entrez.Gene = unlist(mget(probes.limma, hgu133aENTREZID, ifnotfound=NA))
limma.rma = data.frame(data.rma.exprs.varmed,Entrez.Gene, Symbol)



###############
#for genes mapped to multiple probes, keep the probe with the highest effect size and discard rest
###############

t  = var.medlimma$t #t statistic
n1 = 15 # number of cases
n2 = 6 #number of controls
m  = var.medlimma$limma.df  #degrees of freedom as calulated by limma
nprime = (n1 * n2)/(n1 + n2)
d  = t/sqrt(nprime)  #calculate the moderated effect size
cm = (gamma(m/2)) / (gamma((m-1)/2) * sqrt(m/2))
dprime = cm * d
y.varmed = cbind(dprime, limma.rma)
uniqGeneId <- function(x) {
   x = x[order(x$Entrez.Gene, abs(x$dprime), decreasing = TRUE), ]
   entrezID = unique(x$Entrez.Gene)
   id = match(entrezID, x$Entrez.Gene)
   x = x[id[!is.na(id)], ]
   x
}

limma.rma <- uniqGeneId(y.varmed)
limma.rma1 <- limma.rma[,-1]
nrow(limma.rma1) #4898

#-------------------------------


###
#create data frame for use in metaUnion list
###
x = c(15,6)
names(x) <- c("case", "control")

df = data.frame(x)
df <- data.frame(t(df))
rownames(df) <- 1


###
#final file
### 

GSE8397 <- list(expression.mat = limma.rma1, samples.info = df, design = design)
save(GSE8397, file = "GSE8397.Rdata")


##############
#plot volcano plots
##############

library(calibrate)

png('rma.varmed_volcano.png')
with(var.medlimma, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", xlim=c(-2.5,2)))
with(subset(var.medlimma, adj.P.Val<.05 ), points(logFC, -log10(P.Value), pch=20, col="red"))
with(subset(var.medlimma, abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col="orange"))
with(subset(var.medlimma, adj.P.Val<.05 & abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col="green"))
with(subset(var.medlimma, adj.P.Val<.05 & abs(logFC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.8))
abline(h=1,col=3,lty=3)
abline(h=1.30103,col=3,lty=3)
abline(v=0.5849625 ,col=3,lty=3)
abline(v=0.3219281 ,col=3,lty=3)
abline(v=-0.5849625 ,col=3,lty=3)
abline(v=-0.3219281 ,col=3,lty=3)
dev.off()
