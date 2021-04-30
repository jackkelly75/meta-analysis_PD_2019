load(file = "GSE7621.Rdata", verbose = FALSE)
load(file = "GSE8397.Rdata", verbose = FALSE)
load(file = "GSE20141.Rdata", verbose = FALSE)
load(file = "GSE20163.Rdata", verbose = FALSE)
load(file = "GSE20164.Rdata", verbose = FALSE)
load(file = "GSE20292.Rdata", verbose = FALSE)
load(file = "GSE20333.Rdata", verbose = FALSE)

trem <- list(GSE7621 = GSE7621, GSE8397 = GSE8397, GSE20141 = GSE20141,GSE20164 = GSE20164, GSE20292 = GSE20292, GSE20163 = GSE20163, GSE20333= GSE20333)

library(metaUnion)
#then replace the metaAnalysis function with the edited version
trace(metaAnalysis, edit = TRUE)

#run meta-analysis
metanaly <- metaAnalysis(trem, uniqGeneSelMethod="dprime", calWithLimma=TRUE, combinedPval=FALSE, filterData=FALSE)

#add average FC and logFC to the results
metanaly <- data.matrix(metanaly)
temp <- NA
for (i in 1:nrow(metanaly)){
	temp[i] <- mean(metanaly[i,3:9], na.rm = T)
}
metanaly <- as.data.frame(metanaly)
metanaly$PD_Avg_Log_FC <- temp

for (i in 1:nrow(metanaly)){
	metanaly$PD_Avg_FC[i] <- 2^ metanaly$PD_Avg_Log_FC[i]
}

#save file of all genes
all_genes <- metanaly
write.table(all_genes, file = "all_genes.txt", sep = "\t", quote = F, row.names = F)
save(all_genes, file = "all_genes.Rdata")

#find genes that have an FDR corrected pvalue < 0.05
DEGs <- all_genes[all_genes$metaPvalFDR < 0.05,]
write.table(DEGs, file = "DEGs.txt", sep = "\t", quote = F, row.names = F)
save(DEGs, file = "DEGs.Rdata")

#
upreg_DEGs <- DEGs[DEGs$PD_Avg_Log_FC > 0,]
downreg_DEGs <- DEGs[DEGs$PD_Avg_Log_FC < 0,]

write.table(upreg_DEGs, file = "upreg_DEGs.txt", sep = "\t", quote = F, row.names = F)
write.table(downreg_DEGs, file = "downreg_DEGs.txt", sep = "\t", quote = F, row.names = F)