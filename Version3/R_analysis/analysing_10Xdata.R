library(scater)
library(ggplot2)
library(biomaRt)
library(tibble)
library(data.table)

load("../multigenome_10X8x10K_V3/10X8x10K.hsap.full.SCE.Robj")
X10x8x10K.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

X10x8x10K.hsap.count <- as.data.frame(counts(X10x8x10K.hsap))
#X10x8x10K.hsap.count.mapped <- mapIDs(X10x8x10K.hsap.count, "hsap")
X10x8x10K.hsap.metadata <- as.data.frame(colData(X10x8x10K.hsap))
X10x8x10K.hsap.metadata[is.na(X10x8x10K.hsap.metadata)] <- 0
X10x8x10K.hsap.tmappedReads <- rowSums(X10x8x10K.hsap.metadata[, c(2:4,6)])
X10x8x10K.hsap.mappedReadsPerc <- (X10x8x10K.hsap.tmappedReads/X10x8x10K.hsap.metadata$nTReads)*100
X10x8x10K.hsap.metadata <- add_column(X10x8x10K.hsap.metadata, mappedReads= X10x8x10K.hsap.tmappedReads, mappedReadsPerc= X10x8x10K.hsap.mappedReadsPerc, .after = 1 )
X10x8x10K.hsap.cells <- which(X10x8x10K.hsap.metadata$Species == "Human")
X10x8x10K.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(X10x8x10K.hsap.count[,X10x8x10K.hsap.cells])), colData = X10x8x10K.hsap.metadata[X10x8x10K.hsap.cells,c(1:15)])

#X10x8x10K.hsap.is.mito <- grep("^MT-", rownames(counts(X10x8x10K.hsap.new)))
#X10x8x10K.hsap.new <- calculateQCMetrics(X10x8x10K.hsap.new, feature_controls=list(Mt=X10x8x10K.hsap.is.mito))
X10x8x10K.hsap.nTReads.drop <- which(isOutlier(X10x8x10K.hsap.new$mappedReads, nmads=2, type="lower", log=TRUE))
X10x8x10K.hsap.mappedPct.drop <- which(X10x8x10K.hsap.new$mappedReadsPerc < 65)
#X10x8x10K.hsap.mito.drop <- which(isOutlier(X10x8x10K.hsap.new$log10_total_counts_Mt, nmads=2, type="higher"))
#plot(density(colData(X10x8x10K.hsap.new)[-X10x8x10K.hsap.mito.drop, "log10_total_counts_Mt"]))
X10x8x10K.hsap.final <- X10x8x10K.hsap.new[,-c(X10x8x10K.hsap.nTReads.drop, X10x8x10K.hsap.mappedPct.drop)]
X10x8x10K.hsap.fin.metadata <- colData(X10x8x10K.hsap.final)
save(X10x8x10K.hsap.fin.metadata, file = "10X8x10K_hsap_fin_metadata.RData")

hsap.nread.numi.df <- data.frame(X10x8x10K.hsap.fin.metadata[, c("nTReads", "nUMIs", "Library")], tech = rep("X10x8x10K", nrow(X10x8x10K.hsap.fin.metadata)))
pdf("10X8x10K_hsap_nReadsVSnUMIs.pdf")
ggplot(hsap.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_point(aes(color=Library)) + geom_smooth(aes(color=tech)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))# + geom_line(aes(color=tech))#
dev.off()
hsap.nread.genes.df <- data.frame(X10x8x10K.hsap.fin.metadata[, c("nTReads", "nGenes", "Library")], tech = rep("X10x8x10K", nrow(X10x8x10K.hsap.fin.metadata)))
pdf("10X8x10K_hsap_nReadsVSnGenes.pdf")
ggplot(hsap.nread.genes.df, aes(x=nTReads, y=nGenes, group=tech)) + geom_point(aes(color=Library)) + geom_smooth(aes(color=tech)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))# + geom_line(aes(color=tech))#
dev.off()