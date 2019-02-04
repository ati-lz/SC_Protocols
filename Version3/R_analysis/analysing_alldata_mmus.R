library(scater)
library(ggplot2)
library(biomaRt)
library(tibble)
library(data.table)
#require(minpack.lm)

# Functions ####
load("Biomart_mmus_mapping_table.RData")
load("Biomart_mmus_mapping_table.RData")
mapIDs <- function (count.mat, species){
  ens_ids <- rownames(count.mat)
  mat.rownames <- as.character(lapply(ens_ids, function (x) unlist(strsplit(x, split = ".", fixed = T))[1]))
  mat.dup.rows <- which(duplicated(mat.rownames))
  if (length(mat.dup.rows) != 0){
    mat.rownames <- mat.rownames[-mat.dup.rows]
    count.mat <- count.mat[-mat.dup.rows,]
  }
  if (species == "hsap"){
    hsap.genes.with.ids <- intersect(mat.rownames, rownames(hsap.GID.mapping.final))
    hsap.rownames.mapped <- mat.rownames
    names(hsap.rownames.mapped) <- hsap.rownames.mapped
    hsap.rownames.mapped[hsap.genes.with.ids] <- hsap.GID.mapping.final[hsap.genes.with.ids, "hgnc_symbol"]
    rownames(count.mat) <- hsap.rownames.mapped
  }
  else if (species == "mmus"){
    mmus.genes.with.ids <- intersect(mat.rownames, rownames(mmus.GID.mapping.final))
    mmus.rownames.mapped <- mat.rownames
    names(mmus.rownames.mapped) <- mmus.rownames.mapped
    mmus.rownames.mapped[mmus.genes.with.ids] <- mmus.GID.mapping.final[mmus.genes.with.ids, "mgi_symbol"]
    rownames(count.mat) <- mmus.rownames.mapped
  }
  return(count.mat)
  
}
#### end ####


# MARSseq mmus Preparation ====
print("MARSseq is Running...")
load("../SCE_Robjects/MARSseq.mmus.full.SCE.Robj")
MARSseq.mmus <- full.SCE.mmus
rm(full.SCE.mmus)

MARSseq.mmus.count <- as.data.frame(counts(MARSseq.mmus))
MARSseq.mmus.count.mapped <- mapIDs(MARSseq.mmus.count, "mmus")
MARSseq.mmus.metadata <- as.data.frame(colData(MARSseq.mmus))
MARSseq.mmus.metadata[is.na(MARSseq.mmus.metadata)] <- 0
MARSseq.mmus.tmappedReads <- rowSums(MARSseq.mmus.metadata[, c(2:4,6,7)])
MARSseq.mmus.mappedReadsPerc <- (MARSseq.mmus.tmappedReads/MARSseq.mmus.metadata$nTReads)*100
MARSseq.mmus.metadata <- add_column(MARSseq.mmus.metadata, mappedReads= MARSseq.mmus.tmappedReads, mappedReadsPerc= MARSseq.mmus.mappedReadsPerc, .after = 1 )
MARSseq.mmus.cells <- which(MARSseq.mmus.metadata$Species == "Mouse")
MARSseq.mmus.new <- SingleCellExperiment(assays = list(counts = as.matrix(MARSseq.mmus.count.mapped[,MARSseq.mmus.cells])), colData = MARSseq.mmus.metadata[MARSseq.mmus.cells,c(1:15)])
rm(MARSseq.mmus.count)
rm(MARSseq.mmus.count.mapped)
rm(MARSseq.mmus.metadata)
MARSseq.mmus.is.mito <- grep("^mt-", rownames(counts(MARSseq.mmus.new)))
MARSseq.mmus.new <- calculateQCMetrics(MARSseq.mmus.new, feature_controls=list(Mt=MARSseq.mmus.is.mito))
MARSseq.mmus.nTReads.drop <- which(isOutlier(MARSseq.mmus.new$mappedReads, nmads=3, type="lower", log=TRUE))
MARSseq.mmus.mappedPct.drop <- which(MARSseq.mmus.new$mappedReadsPerc < 65)
MARSseq.mmus.mito.drop <- which(isOutlier(MARSseq.mmus.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(MARSseq.mmus.new)[-MARSseq.mmus.mito.drop, "log10_total_counts_Mt"]))
MARSseq.mmus.final <- MARSseq.mmus.new[,-c(MARSseq.mmus.nTReads.drop, MARSseq.mmus.mappedPct.drop,MARSseq.mmus.mito.drop)]
MARSseq.mmus.fin.metadata <- colData(MARSseq.mmus.final)
#### end ####

# CELseq2 mmus Preparation ====
print("CELseq2 is Running...")
load("../SCE_Robjects/CELseq2.mmus.full.SCE.Robj")
CELseq2.mmus <- full.SCE.mmus
rm(full.SCE.mmus)

CELseq2.mmus.count <- as.data.frame(counts(CELseq2.mmus))
CELseq2.mmus.count.mapped <- mapIDs(CELseq2.mmus.count, "mmus")
CELseq2.mmus.metadata <- as.data.frame(colData(CELseq2.mmus))
CELseq2.mmus.metadata[is.na(CELseq2.mmus.metadata)] <- 0
CELseq2.mmus.tmappedReads <- rowSums(CELseq2.mmus.metadata[, c(2:4,6,7)])
CELseq2.mmus.mappedReadsPerc <- (CELseq2.mmus.tmappedReads/CELseq2.mmus.metadata$nTReads)*100
CELseq2.mmus.metadata <- add_column(CELseq2.mmus.metadata, mappedReads= CELseq2.mmus.tmappedReads, mappedReadsPerc= CELseq2.mmus.mappedReadsPerc, .after = 1 )
CELseq2.mmus.cells <- which(CELseq2.mmus.metadata$Species == "Mouse")
CELseq2.mmus.new <- SingleCellExperiment(assays = list(counts = as.matrix(CELseq2.mmus.count.mapped[,CELseq2.mmus.cells])), colData = CELseq2.mmus.metadata[CELseq2.mmus.cells,c(1:15)])
rm(CELseq2.mmus.count)
rm(CELseq2.mmus.count.mapped)
#rm(CELseq2.mmus.metadata)
CELseq2.mmus.is.mito <- grep("^mt-", rownames(counts(CELseq2.mmus.new)))
CELseq2.mmus.new <- calculateQCMetrics(CELseq2.mmus.new, feature_controls=list(Mt=CELseq2.mmus.is.mito))
CELseq2.mmus.nTReads.drop <- which(isOutlier(CELseq2.mmus.new$mappedReads, nmads=3, type="lower", log=TRUE))
CELseq2.mmus.mappedPct.drop <- which(CELseq2.mmus.new$mappedReadsPerc < 65)
CELseq2.mmus.mito.drop <- which(isOutlier(CELseq2.mmus.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(CELseq2.mmus.new)[-CELseq2.mmus.mito.drop, "log10_total_counts_Mt"]))
CELseq2.mmus.final <- CELseq2.mmus.new[,-c(CELseq2.mmus.nTReads.drop, CELseq2.mmus.mappedPct.drop,CELseq2.mmus.mito.drop)]
CELseq2.mmus.fin.metadata <- colData(CELseq2.mmus.final)
#### end ####

# QUARTZseq mmus Preparation ====
print("QUARTZseq is Running...")
load("../SCE_Robjects/QUARTZseq.mmus.full.SCE.Robj")
QUARTZseq.mmus <- full.SCE.mmus
rm(full.SCE.mmus)

QUARTZseq.mmus.count <- as.data.frame(counts(QUARTZseq.mmus))
QUARTZseq.mmus.count.mapped <- mapIDs(QUARTZseq.mmus.count, "mmus")
QUARTZseq.mmus.metadata <- as.data.frame(colData(QUARTZseq.mmus))
QUARTZseq.mmus.metadata[is.na(QUARTZseq.mmus.metadata)] <- 0
QUARTZseq.mmus.tmappedReads <- rowSums(QUARTZseq.mmus.metadata[, c(2:4,6,7)])
QUARTZseq.mmus.mappedReadsPerc <- (QUARTZseq.mmus.tmappedReads/QUARTZseq.mmus.metadata$nTReads)*100
QUARTZseq.mmus.metadata <- add_column(QUARTZseq.mmus.metadata, mappedReads= QUARTZseq.mmus.tmappedReads, mappedReadsPerc= QUARTZseq.mmus.mappedReadsPerc, .after = 1 )
QUARTZseq.mmus.cells <- which(QUARTZseq.mmus.metadata$Species == "Mouse")
QUARTZseq.mmus.new <- SingleCellExperiment(assays = list(counts = as.matrix(QUARTZseq.mmus.count.mapped[,QUARTZseq.mmus.cells])), colData = QUARTZseq.mmus.metadata[QUARTZseq.mmus.cells,c(1:15)])
rm(QUARTZseq.mmus.count)
rm(QUARTZseq.mmus.count.mapped)
#rm(QUARTZseq.mmus.metadata)
QUARTZseq.mmus.is.mito <- grep("^mt-", rownames(counts(QUARTZseq.mmus.new)))
QUARTZseq.mmus.new <- calculateQCMetrics(QUARTZseq.mmus.new, feature_controls=list(Mt=QUARTZseq.mmus.is.mito))
QUARTZseq.mmus.nTReads.drop <- which(isOutlier(QUARTZseq.mmus.new$mappedReads, nmads=3, type="lower", log=TRUE))
QUARTZseq.mmus.mappedPct.drop <- which(QUARTZseq.mmus.new$mappedReadsPerc < 65)
QUARTZseq.mmus.mito.drop <- which(isOutlier(QUARTZseq.mmus.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(QUARTZseq.mmus.new)[-QUARTZseq.mmus.mito.drop, "log10_total_counts_Mt"]))
QUARTZseq.mmus.final <- QUARTZseq.mmus.new[,-c(QUARTZseq.mmus.nTReads.drop, QUARTZseq.mmus.mappedPct.drop,QUARTZseq.mmus.mito.drop)]
QUARTZseq.mmus.fin.metadata <- colData(QUARTZseq.mmus.final)
#### end ####

# Dropseq mmus Preparation ====
print("Dropseq is Running...")
load("../SCE_Robjects/Dropseq.mmus.full.SCE.1000cellsPerPool.Robj")
Dropseq.mmus <- full.SCE.mmus
rm(full.SCE.mmus)

Dropseq.mmus.count <- as.data.frame(counts(Dropseq.mmus))
Dropseq.mmus.count.mapped <- mapIDs(Dropseq.mmus.count, "mmus")
Dropseq.mmus.metadata <- as.data.frame(colData(Dropseq.mmus))
Dropseq.mmus.metadata[is.na(Dropseq.mmus.metadata)] <- 0
Dropseq.mmus.tmappedReads <- rowSums(Dropseq.mmus.metadata[, c(2:4,6,7)])
Dropseq.mmus.mappedReadsPerc <- (Dropseq.mmus.tmappedReads/Dropseq.mmus.metadata$nTReads)*100
Dropseq.mmus.metadata <- add_column(Dropseq.mmus.metadata, mappedReads= Dropseq.mmus.tmappedReads, mappedReadsPerc= Dropseq.mmus.mappedReadsPerc, .after = 1 )
Dropseq.mmus.cells <- which(Dropseq.mmus.metadata$Species == "Mouse")
Dropseq.mmus.new <- SingleCellExperiment(assays = list(counts = as.matrix(Dropseq.mmus.count.mapped[,Dropseq.mmus.cells])), colData = Dropseq.mmus.metadata[Dropseq.mmus.cells,c(1:15)])
rm(Dropseq.mmus.count)
rm(Dropseq.mmus.count.mapped)

Dropseq.mmus.is.mito <- grep("^mt-", rownames(counts(Dropseq.mmus.new)))
Dropseq.mmus.new <- calculateQCMetrics(Dropseq.mmus.new, feature_controls=list(Mt=Dropseq.mmus.is.mito))
Dropseq.mmus.nTReads.drop <- which(isOutlier(Dropseq.mmus.new$mappedReads, nmads=3, type="lower", log=TRUE))
Dropseq.mmus.mappedPct.drop <- which(Dropseq.mmus.new$mappedReadsPerc < 65)
Dropseq.mmus.mito.drop <- which(isOutlier(Dropseq.mmus.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(Dropseq.mmus.new)[-Dropseq.mmus.mito.drop, "log10_total_counts_Mt"]))
Dropseq.mmus.final <- Dropseq.mmus.new[,-c(Dropseq.mmus.nTReads.drop, Dropseq.mmus.mappedPct.drop,Dropseq.mmus.mito.drop)]
Dropseq.mmus.fin.metadata <- colData(Dropseq.mmus.final)
#### end ####

# SCRBseq mmus Preparation ====
print("SCRBseq is Running...")
load("../SCE_Robjects/SCRBseq.mmus.full.SCE.Robj")
SCRBseq.mmus <- full.SCE.mmus
rm(full.SCE.mmus)

SCRBseq.mmus.count <- as.data.frame(counts(SCRBseq.mmus))
SCRBseq.mmus.count.mapped <- mapIDs(SCRBseq.mmus.count, "mmus")
SCRBseq.mmus.metadata <- as.data.frame(colData(SCRBseq.mmus))
SCRBseq.mmus.metadata[is.na(SCRBseq.mmus.metadata)] <- 0
SCRBseq.mmus.tmappedReads <- rowSums(SCRBseq.mmus.metadata[, c(2:4,6,7)])
SCRBseq.mmus.mappedReadsPerc <- (SCRBseq.mmus.tmappedReads/SCRBseq.mmus.metadata$nTReads)*100
SCRBseq.mmus.metadata <- add_column(SCRBseq.mmus.metadata, mappedReads= SCRBseq.mmus.tmappedReads, mappedReadsPerc= SCRBseq.mmus.mappedReadsPerc, .after = 1 )
SCRBseq.mmus.cells <- which(SCRBseq.mmus.metadata$Species == "Mouse")
SCRBseq.mmus.new <- SingleCellExperiment(assays = list(counts = as.matrix(SCRBseq.mmus.count.mapped[,SCRBseq.mmus.cells])), colData = SCRBseq.mmus.metadata[SCRBseq.mmus.cells,c(1:15)])

SCRBseq.mmus.is.mito <- grep("^mt-", rownames(counts(SCRBseq.mmus.new)))
SCRBseq.mmus.new <- calculateQCMetrics(SCRBseq.mmus.new, feature_controls=list(Mt=SCRBseq.mmus.is.mito))
SCRBseq.mmus.nTReads.drop <- which(isOutlier(SCRBseq.mmus.new$mappedReads, nmads=3, type="lower", log=TRUE))
SCRBseq.mmus.mappedPct.drop <- which(SCRBseq.mmus.new$mappedReadsPerc < 65)
SCRBseq.mmus.mito.drop <- which(isOutlier(SCRBseq.mmus.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(SCRBseq.mmus.new)[-SCRBseq.mmus.mito.drop, "log10_total_counts_Mt"]))
SCRBseq.mmus.final <- SCRBseq.mmus.new[,-c(SCRBseq.mmus.nTReads.drop, SCRBseq.mmus.mappedPct.drop,SCRBseq.mmus.mito.drop)]
SCRBseq.mmus.fin.metadata <- colData(SCRBseq.mmus.final)
#### end ####


# SeqwellV2 mmus Preparation ====
print("SeqwellV2 is Running...")
load("../SCE_Robjects/SeqWellV2.mmus.full.SCE.Robj")
SeqwellV2.mmus <- full.SCE.mmus
rm(full.SCE.mmus)

SeqwellV2.mmus.count <- as.data.frame(counts(SeqwellV2.mmus))
SeqwellV2.mmus.count.mapped <- mapIDs(SeqwellV2.mmus.count, "mmus")
SeqwellV2.mmus.metadata <- as.data.frame(colData(SeqwellV2.mmus))
SeqwellV2.mmus.metadata[is.na(SeqwellV2.mmus.metadata)] <- 0
SeqwellV2.mmus.tmappedReads <- rowSums(SeqwellV2.mmus.metadata[, c(2:4,6,7)])
SeqwellV2.mmus.mappedReadsPerc <- (SeqwellV2.mmus.tmappedReads/SeqwellV2.mmus.metadata$nTReads)*100
SeqwellV2.mmus.metadata <- add_column(SeqwellV2.mmus.metadata, mappedReads= SeqwellV2.mmus.tmappedReads, mappedReadsPerc= SeqwellV2.mmus.mappedReadsPerc, .after = 1 )
SeqwellV2.mmus.cells <- which(SeqwellV2.mmus.metadata$Species == "Mouse")
SeqwellV2.mmus.new <- SingleCellExperiment(assays = list(counts = as.matrix(SeqwellV2.mmus.count.mapped[,SeqwellV2.mmus.cells])), colData = SeqwellV2.mmus.metadata[SeqwellV2.mmus.cells,c(1:15)])

SeqwellV2.mmus.is.mito <- grep("^mt-", rownames(counts(SeqwellV2.mmus.new)))
SeqwellV2.mmus.new <- calculateQCMetrics(SeqwellV2.mmus.new, feature_controls=list(Mt=SeqwellV2.mmus.is.mito))
SeqwellV2.mmus.nTReads.drop <- which(isOutlier(SeqwellV2.mmus.new$mappedReads, nmads=3, type="lower", log=TRUE))
SeqwellV2.mmus.mappedPct.drop <- which(SeqwellV2.mmus.new$mappedReadsPerc < 65)
SeqwellV2.mmus.mito.drop <- which(isOutlier(SeqwellV2.mmus.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(SeqwellV2.mmus.new)[-SeqwellV2.mmus.mito.drop, "log10_total_counts_Mt"]))
SeqwellV2.mmus.final <- SeqwellV2.mmus.new[,-c(SeqwellV2.mmus.nTReads.drop, SeqwellV2.mmus.mappedPct.drop,SeqwellV2.mmus.mito.drop)]
SeqwellV2.mmus.fin.metadata <- colData(SeqwellV2.mmus.final)
#### end ####


# SeqwellV1 mmus Preparation ====
print("SeqwellV1 is Running...")
load("../SCE_Robjects/SeqWellV1.mmus.full.SCE.Robj")
SeqwellV1.mmus <- full.SCE.mmus
rm(full.SCE.mmus)

SeqwellV1.mmus.count <- as.data.frame(counts(SeqwellV1.mmus))
SeqwellV1.mmus.count.mapped <- mapIDs(SeqwellV1.mmus.count, "mmus")
SeqwellV1.mmus.metadata <- as.data.frame(colData(SeqwellV1.mmus))
SeqwellV1.mmus.metadata[is.na(SeqwellV1.mmus.metadata)] <- 0
SeqwellV1.mmus.tmappedReads <- rowSums(SeqwellV1.mmus.metadata[, c(2:4,6)])
SeqwellV1.mmus.mappedReadsPerc <- (SeqwellV1.mmus.tmappedReads/SeqwellV1.mmus.metadata$nTReads)*100
SeqwellV1.mmus.metadata <- add_column(SeqwellV1.mmus.metadata, mappedReads= SeqwellV1.mmus.tmappedReads, mappedReadsPerc= SeqwellV1.mmus.mappedReadsPerc, .after = 1 )
SeqwellV1.mmus.cells <- which(SeqwellV1.mmus.metadata$Species == "Mouse")
SeqwellV1.mmus.new <- SingleCellExperiment(assays = list(counts = as.matrix(SeqwellV1.mmus.count.mapped[,SeqwellV1.mmus.cells])), colData = SeqwellV1.mmus.metadata[SeqwellV1.mmus.cells,c(1:15)])

SeqwellV1.mmus.is.mito <- grep("^mt-", rownames(counts(SeqwellV1.mmus.new)))
SeqwellV1.mmus.new <- calculateQCMetrics(SeqwellV1.mmus.new, feature_controls=list(Mt=SeqwellV1.mmus.is.mito))
SeqwellV1.mmus.nTReads.drop <- which(isOutlier(SeqwellV1.mmus.new$mappedReads, nmads=2, type="lower", log=TRUE))
SeqwellV1.mmus.mappedPct.drop <- which(SeqwellV1.mmus.new$mappedReadsPerc < 65)
SeqwellV1.mmus.mito.drop <- which(isOutlier(SeqwellV1.mmus.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(SeqwellV1.mmus.new)[-SeqwellV1.mmus.mito.drop, "log10_total_counts_Mt"]))
SeqwellV1.mmus.final <- SeqwellV1.mmus.new[,-c(SeqwellV1.mmus.nTReads.drop, SeqwellV1.mmus.mappedPct.drop,SeqwellV1.mmus.mito.drop)]
SeqwellV1.mmus.fin.metadata <- colData(SeqwellV1.mmus.final)
#### end ####


# Nuclei10X mmus Preparation ====
print("Nuclei10X is Running...")
load("../SCE_Robjects/Nuclei10X.mmus.full.SCE.Robj")
Nuclei10X.mmus <- full.SCE.mmus
rm(full.SCE.mmus)

Nuclei10X.mmus.count <- as.data.frame(counts(Nuclei10X.mmus))
Nuclei10X.mmus.count.mapped <- mapIDs(Nuclei10X.mmus.count, "mmus")
Nuclei10X.mmus.metadata <- as.data.frame(colData(Nuclei10X.mmus))
Nuclei10X.mmus.metadata[is.na(Nuclei10X.mmus.metadata)] <- 0
Nuclei10X.mmus.tmappedReads <- rowSums(Nuclei10X.mmus.metadata[, c(2:4,6)])
Nuclei10X.mmus.mappedReadsPerc <- (Nuclei10X.mmus.tmappedReads/Nuclei10X.mmus.metadata$nTReads)*100
Nuclei10X.mmus.metadata <- add_column(Nuclei10X.mmus.metadata, mappedReads= Nuclei10X.mmus.tmappedReads, mappedReadsPerc= Nuclei10X.mmus.mappedReadsPerc, .after = 1 )
Nuclei10X.mmus.cells <- which(Nuclei10X.mmus.metadata$Species == "Mouse")
Nuclei10X.mmus.new <- SingleCellExperiment(assays = list(counts = as.matrix(Nuclei10X.mmus.count.mapped[,Nuclei10X.mmus.cells])), colData = Nuclei10X.mmus.metadata[Nuclei10X.mmus.cells,c(1:15)])

Nuclei10X.mmus.is.mito <- grep("^mt-", rownames(counts(Nuclei10X.mmus.new)))
Nuclei10X.mmus.new <- calculateQCMetrics(Nuclei10X.mmus.new, feature_controls=list(Mt=Nuclei10X.mmus.is.mito))
Nuclei10X.mmus.nTReads.drop <- which(isOutlier(Nuclei10X.mmus.new$mappedReads, nmads=2, type="lower", log=TRUE))
Nuclei10X.mmus.mappedPct.drop <- which(Nuclei10X.mmus.new$mappedReadsPerc < 65)
Nuclei10X.mmus.mito.drop <- which(isOutlier(Nuclei10X.mmus.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(Nuclei10X.mmus.new)[-Nuclei10X.mmus.mito.drop, "log10_total_counts_Mt"]))
Nuclei10X.mmus.final <- Nuclei10X.mmus.new[,-c(Nuclei10X.mmus.nTReads.drop, Nuclei10X.mmus.mappedPct.drop,Nuclei10X.mmus.mito.drop)]
Nuclei10X.mmus.fin.metadata <- colData(Nuclei10X.mmus.final)
#### end ####


# ICELL8 mmus Preparation ====
print("ICELL8 is Running...")
load("../SCE_Robjects/ICELL8.mmus.full.SCE.Robj")
ICELL8.mmus <- full.SCE.mmus
rm(full.SCE.mmus)

ICELL8.mmus.count <- as.data.frame(counts(ICELL8.mmus))
ICELL8.mmus.count.mapped <- mapIDs(ICELL8.mmus.count, "mmus")
ICELL8.mmus.metadata <- as.data.frame(colData(ICELL8.mmus))
ICELL8.mmus.metadata[is.na(ICELL8.mmus.metadata)] <- 0
ICELL8.mmus.tmappedReads <- rowSums(ICELL8.mmus.metadata[, c(2:4,6)])
ICELL8.mmus.mappedReadsPerc <- (ICELL8.mmus.tmappedReads/ICELL8.mmus.metadata$nTReads)*100
ICELL8.mmus.metadata <- add_column(ICELL8.mmus.metadata, mappedReads= ICELL8.mmus.tmappedReads, mappedReadsPerc= ICELL8.mmus.mappedReadsPerc, .after = 1 )
ICELL8.mmus.cells <- which(ICELL8.mmus.metadata$Species == "Mouse")
ICELL8.mmus.new <- SingleCellExperiment(assays = list(counts = as.matrix(ICELL8.mmus.count.mapped[,ICELL8.mmus.cells])), colData = ICELL8.mmus.metadata[ICELL8.mmus.cells,c(1:15)])

ICELL8.mmus.is.mito <- grep("^mt-", rownames(counts(ICELL8.mmus.new)))
ICELL8.mmus.new <- calculateQCMetrics(ICELL8.mmus.new, feature_controls=list(Mt=ICELL8.mmus.is.mito))
ICELL8.mmus.nTReads.drop <- which(isOutlier(ICELL8.mmus.new$mappedReads, nmads=2, type="lower", log=TRUE))
ICELL8.mmus.mappedPct.drop <- which(ICELL8.mmus.new$mappedReadsPerc < 65)
ICELL8.mmus.mito.drop <- which(isOutlier(ICELL8.mmus.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(ICELL8.mmus.new)[-ICELL8.mmus.mito.drop, "log10_total_counts_Mt"]))
ICELL8.mmus.final <- ICELL8.mmus.new[,-c(ICELL8.mmus.nTReads.drop, ICELL8.mmus.mappedPct.drop,ICELL8.mmus.mito.drop)]
ICELL8.mmus.fin.metadata <- colData(ICELL8.mmus.final)
                                                                                                       
#### end ####


# ddSEQ mmus Preparation ====
print("ddSEQ is Running...")
load("../SCE_Robjects/ddSEQ.mmus.full.SCE.Robj")
ddSEQ.mmus <- full.SCE.mmus
rm(full.SCE.mmus)

ddSEQ.mmus.count <- as.data.frame(counts(ddSEQ.mmus))
ddSEQ.mmus.count.mapped <- mapIDs(ddSEQ.mmus.count, "mmus")
ddSEQ.mmus.metadata <- as.data.frame(colData(ddSEQ.mmus))
ddSEQ.mmus.metadata[is.na(ddSEQ.mmus.metadata)] <- 0
ddSEQ.mmus.tmappedReads <- rowSums(ddSEQ.mmus.metadata[, c(2:4,6)])
ddSEQ.mmus.mappedReadsPerc <- (ddSEQ.mmus.tmappedReads/ddSEQ.mmus.metadata$nTReads)*100
ddSEQ.mmus.metadata <- add_column(ddSEQ.mmus.metadata, mappedReads= ddSEQ.mmus.tmappedReads, mappedReadsPerc= ddSEQ.mmus.mappedReadsPerc, .after = 1 )
ddSEQ.mmus.cells <- which(ddSEQ.mmus.metadata$Species == "Mouse")
ddSEQ.mmus.new <- SingleCellExperiment(assays = list(counts = as.matrix(ddSEQ.mmus.count.mapped[,ddSEQ.mmus.cells])), colData = ddSEQ.mmus.metadata[ddSEQ.mmus.cells,c(1:15)])

ddSEQ.mmus.is.mito <- grep("^mt-", rownames(counts(ddSEQ.mmus.new)))
ddSEQ.mmus.new <- calculateQCMetrics(ddSEQ.mmus.new, feature_controls=list(Mt=ddSEQ.mmus.is.mito))
ddSEQ.mmus.nTReads.drop <- which(isOutlier(ddSEQ.mmus.new$mappedReads, nmads=2, type="lower", log=TRUE))
ddSEQ.mmus.mappedPct.drop <- which(ddSEQ.mmus.new$mappedReadsPerc < 65)
ddSEQ.mmus.mito.drop <- which(isOutlier(ddSEQ.mmus.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(ddSEQ.mmus.new)[-ddSEQ.mmus.mito.drop, "log10_total_counts_Mt"]))
ddSEQ.mmus.final <- ddSEQ.mmus.new[,-c(ddSEQ.mmus.nTReads.drop, ddSEQ.mmus.mappedPct.drop,ddSEQ.mmus.mito.drop)]
ddSEQ.mmus.fin.metadata <- colData(ddSEQ.mmus.final)
#### end ####

# ddSEQexp1 mmus Preparation ====
print("ddSEQexp1 is Running...")
load("../SCE_Robjects/ddSEQexp1.mmus.full.SCE.Robj")
ddSEQexp1.mmus <- full.SCE.mmus
rm(full.SCE.mmus)

ddSEQexp1.mmus.count <- as.data.frame(counts(ddSEQexp1.mmus))
ddSEQexp1.mmus.count.mapped <- mapIDs(ddSEQexp1.mmus.count, "mmus")
ddSEQexp1.mmus.metadata <- as.data.frame(colData(ddSEQexp1.mmus))
ddSEQexp1.mmus.metadata[is.na(ddSEQexp1.mmus.metadata)] <- 0
ddSEQexp1.mmus.tmappedReads <- rowSums(ddSEQexp1.mmus.metadata[, c(2:4,6)])
ddSEQexp1.mmus.mappedReadsPerc <- (ddSEQexp1.mmus.tmappedReads/ddSEQexp1.mmus.metadata$nTReads)*100
ddSEQexp1.mmus.metadata <- add_column(ddSEQexp1.mmus.metadata, mappedReads= ddSEQexp1.mmus.tmappedReads, mappedReadsPerc= ddSEQexp1.mmus.mappedReadsPerc, .after = 1 )
ddSEQexp1.mmus.cells <- which(ddSEQexp1.mmus.metadata$Species == "Mouse")
ddSEQexp1.mmus.new <- SingleCellExperiment(assays = list(counts = as.matrix(ddSEQexp1.mmus.count.mapped[,ddSEQexp1.mmus.cells])), colData = ddSEQexp1.mmus.metadata[ddSEQexp1.mmus.cells,c(1:15)])

ddSEQexp1.mmus.is.mito <- grep("^mt-", rownames(counts(ddSEQexp1.mmus.new)))
ddSEQexp1.mmus.new <- calculateQCMetrics(ddSEQexp1.mmus.new, feature_controls=list(Mt=ddSEQexp1.mmus.is.mito))
ddSEQexp1.mmus.nTReads.drop <- which(isOutlier(ddSEQexp1.mmus.new$mappedReads, nmads=2, type="lower", log=TRUE))
ddSEQexp1.mmus.mappedPct.drop <- which(ddSEQexp1.mmus.new$mappedReadsPerc < 65)
ddSEQexp1.mmus.mito.drop <- which(isOutlier(ddSEQexp1.mmus.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(ddSEQexp1.mmus.new)[-ddSEQexp1.mmus.mito.drop, "log10_total_counts_Mt"]))
ddSEQexp1.mmus.final <- ddSEQexp1.mmus.new[,-c(ddSEQexp1.mmus.nTReads.drop, ddSEQexp1.mmus.mappedPct.drop,ddSEQexp1.mmus.mito.drop)]
ddSEQexp1.mmus.fin.metadata <- colData(ddSEQexp1.mmus.final)
#### end ####


# C1HT mmus Preparation ====
print("C1HT is Running...")
load("../SCE_Robjects/C1HT.mmus.full.SCE.Robj")
C1HT.mmus <- full.SCE.mmus
rm(full.SCE.mmus)

C1HT.mmus.count <- as.data.frame(counts(C1HT.mmus))
C1HT.mmus.count.mapped <- mapIDs(C1HT.mmus.count, "mmus")
C1HT.mmus.metadata <- as.data.frame(colData(C1HT.mmus))
C1HT.mmus.metadata[is.na(C1HT.mmus.metadata)] <- 0
C1HT.mmus.tmappedReads <- rowSums(C1HT.mmus.metadata[, c(2:4,6)])
C1HT.mmus.mappedReadsPerc <- (C1HT.mmus.tmappedReads/C1HT.mmus.metadata$nTReads)*100
C1HT.mmus.metadata <- add_column(C1HT.mmus.metadata, mappedReads= C1HT.mmus.tmappedReads, mappedReadsPerc= C1HT.mmus.mappedReadsPerc, .after = 1 )
C1HT.mmus.cells <- which(C1HT.mmus.metadata$Species == "Mouse")
C1HT.mmus.new <- SingleCellExperiment(assays = list(counts = as.matrix(C1HT.mmus.count.mapped[,C1HT.mmus.cells])), colData = C1HT.mmus.metadata[C1HT.mmus.cells,c(1:15)])

C1HT.mmus.is.mito <- grep("^mt-", rownames(counts(C1HT.mmus.new)))
C1HT.mmus.new <- calculateQCMetrics(C1HT.mmus.new, feature_controls=list(Mt=C1HT.mmus.is.mito))
C1HT.mmus.nTReads.drop <- which(isOutlier(C1HT.mmus.new$mappedReads, nmads=2, type="lower", log=TRUE))
C1HT.mmus.mappedPct.drop <- which(C1HT.mmus.new$mappedReadsPerc < 65)
C1HT.mmus.mito.drop <- which(isOutlier(C1HT.mmus.new$log10_total_counts_Mt, nmads=2, type="higher"))
#plot(density(colData(C1HT.mmus.new)[-C1HT.mmus.mito.drop, "log10_total_counts_Mt"]))
C1HT.mmus.final <- C1HT.mmus.new[,-c(C1HT.mmus.nTReads.drop, C1HT.mmus.mappedPct.drop,C1HT.mmus.mito.drop)]
C1HT.mmus.fin.metadata <- colData(C1HT.mmus.final)

#### end ####



#10X8x10K mmus Preparation ====
print("10X8x10K is Running...")
load("../SCE_Robjects/10X8x10K.mmus.full.SCE.Robj")
X10x8x10K.mmus <- full.SCE.mmus
rm(full.SCE.mmus)

X10x8x10K.mmus.count <- as.data.frame(counts(X10x8x10K.mmus))
X10x8x10K.mmus.count.mapped <- mapIDs(X10x8x10K.mmus.count, "mmus")
X10x8x10K.mmus.metadata <- as.data.frame(colData(X10x8x10K.mmus))
X10x8x10K.mmus.metadata[is.na(X10x8x10K.mmus.metadata)] <- 0
X10x8x10K.mmus.tmappedReads <- rowSums(X10x8x10K.mmus.metadata[, c(2:4,6)])
X10x8x10K.mmus.mappedReadsPerc <- (X10x8x10K.mmus.tmappedReads/X10x8x10K.mmus.metadata$nTReads)*100
X10x8x10K.mmus.metadata <- add_column(X10x8x10K.mmus.metadata, mappedReads= X10x8x10K.mmus.tmappedReads, mappedReadsPerc= X10x8x10K.mmus.mappedReadsPerc, .after = 1 )
X10x8x10K.mmus.cells <- which(X10x8x10K.mmus.metadata$Species == "Mouse")
X10x8x10K.mmus.new <- SingleCellExperiment(assays = list(counts = as.matrix(X10x8x10K.mmus.count.mapped[,X10x8x10K.mmus.cells])), colData = X10x8x10K.mmus.metadata[X10x8x10K.mmus.cells,c(1:15)])

X10x8x10K.mmus.is.mito <- grep("^mt-", rownames(counts(X10x8x10K.mmus.new)))
X10x8x10K.mmus.new <- calculateQCMetrics(X10x8x10K.mmus.new, feature_controls=list(Mt=X10x8x10K.mmus.is.mito))
X10x8x10K.mmus.nTReads.drop <- which(isOutlier(X10x8x10K.mmus.new$mappedReads, nmads=2, type="lower", log=TRUE))
X10x8x10K.mmus.mappedPct.drop <- which(X10x8x10K.mmus.new$mappedReadsPerc < 65)
X10x8x10K.mmus.mito.drop <- which(isOutlier(X10x8x10K.mmus.new$log10_total_counts_Mt, nmads=2, type="higher"))
#plot(density(colData(X10x8x10K.mmus.new)[-X10x8x10K.mmus.mito.drop, "log10_total_counts_Mt"]))
X10x8x10K.mmus.final <- X10x8x10K.mmus.new[,-c(X10x8x10K.mmus.nTReads.drop, X10x8x10K.mmus.mappedPct.drop,X10x8x10K.mmus.mito.drop)]
X10x8x10K.mmus.fin.metadata <- colData(X10x8x10K.mmus.final)

#### end ####


# X10Scilife mmus Preparation ====
print("X10Scilife is Running...")
load("../SCE_Robjects/10XScilife.mmus.full.SCE.Robj")
X10Scilife.mmus <- full.SCE.mmus
rm(full.SCE.mmus)

X10Scilife.mmus.count <- as.data.frame(counts(X10Scilife.mmus))
X10Scilife.mmus.count.mapped <- mapIDs(X10Scilife.mmus.count, "mmus")
X10Scilife.mmus.metadata <- as.data.frame(colData(X10Scilife.mmus))
X10Scilife.mmus.metadata[is.na(X10Scilife.mmus.metadata)] <- 0
X10Scilife.mmus.tmappedReads <- rowSums(X10Scilife.mmus.metadata[, c(2:4,6)])
X10Scilife.mmus.mappedReadsPerc <- (X10Scilife.mmus.tmappedReads/X10Scilife.mmus.metadata$nTReads)*100
X10Scilife.mmus.metadata <- add_column(X10Scilife.mmus.metadata, mappedReads= X10Scilife.mmus.tmappedReads, mappedReadsPerc= X10Scilife.mmus.mappedReadsPerc, .after = 1 )
X10Scilife.mmus.cells <- which(X10Scilife.mmus.metadata$Species == "Mouse")
X10Scilife.mmus.new <- SingleCellExperiment(assays = list(counts = as.matrix(X10Scilife.mmus.count.mapped[,X10Scilife.mmus.cells])), colData = X10Scilife.mmus.metadata[X10Scilife.mmus.cells,c(1:15)])

X10Scilife.mmus.is.mito <- grep("^mt-", rownames(counts(X10Scilife.mmus.new)))
X10Scilife.mmus.new <- calculateQCMetrics(X10Scilife.mmus.new, feature_controls=list(Mt=X10Scilife.mmus.is.mito))
X10Scilife.mmus.nTReads.drop <- which(isOutlier(X10Scilife.mmus.new$mappedReads, nmads=2, type="lower", log=TRUE))
X10Scilife.mmus.mappedPct.drop <- which(X10Scilife.mmus.new$mappedReadsPerc < 65)
X10Scilife.mmus.mito.drop <- which(isOutlier(X10Scilife.mmus.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(X10Scilife.mmus.new)[-X10Scilife.mmus.mito.drop, "log10_total_counts_Mt"]))
X10Scilife.mmus.final <- X10Scilife.mmus.new[,-c(X10Scilife.mmus.nTReads.drop, X10Scilife.mmus.mappedPct.drop,X10Scilife.mmus.mito.drop)]
X10Scilife.mmus.fin.metadata <- colData(X10Scilife.mmus.final)

#### end ####

# CB1 mmus Preparation ====
print("CB1 is Running...")
load("../SCE_Robjects/1CB.mmus.full.SCE.Robj")
CB1.mmus <- full.SCE.mmus
rm(full.SCE.mmus)

CB1.mmus.count <- as.data.frame(counts(CB1.mmus))
CB1.mmus.count.mapped <- mapIDs(CB1.mmus.count, "mmus")
CB1.mmus.metadata <- as.data.frame(colData(CB1.mmus))
CB1.mmus.metadata[is.na(CB1.mmus.metadata)] <- 0
CB1.mmus.tmappedReads <- rowSums(CB1.mmus.metadata[, c(2:4,6)])
CB1.mmus.mappedReadsPerc <- (CB1.mmus.tmappedReads/CB1.mmus.metadata$nTReads)*100
CB1.mmus.metadata <- add_column(CB1.mmus.metadata, mappedReads= CB1.mmus.tmappedReads, mappedReadsPerc= CB1.mmus.mappedReadsPerc, .after = 1 )
CB1.mmus.cells <- which(CB1.mmus.metadata$Species == "Mouse")
CB1.mmus.new <- SingleCellExperiment(assays = list(counts = as.matrix(CB1.mmus.count.mapped[,CB1.mmus.cells])), colData = CB1.mmus.metadata[CB1.mmus.cells,c(1:15)])

CB1.mmus.is.mito <- grep("^mt-", rownames(counts(CB1.mmus.new)))
CB1.mmus.new <- calculateQCMetrics(CB1.mmus.new, feature_controls=list(Mt=CB1.mmus.is.mito))
CB1.mmus.nTReads.drop <- which(isOutlier(CB1.mmus.new$mappedReads, nmads=2, type="lower", log=TRUE))
CB1.mmus.mappedPct.drop <- which(CB1.mmus.new$mappedReadsPerc < 65)
CB1.mmus.mito.drop <- which(isOutlier(CB1.mmus.new$log10_total_counts_Mt, nmads=2, type="higher"))
#plot(density(colData(CB1.mmus.new)[-CB1.mmus.mito.drop, "log10_total_counts_Mt"]))
CB1.mmus.final <- CB1.mmus.new[,-c(CB1.mmus.nTReads.drop, CB1.mmus.mappedPct.drop,CB1.mmus.mito.drop)]
CB1.mmus.fin.metadata <- colData(CB1.mmus.final)

#### end ####




#library(data.table)
# Read distribution in each sample in each tech ====
#print("Read distribution in each sample in each tech")
#mmus.read.percentages <- data.frame(MARSseq.mmus.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("MARSseq", nrow(MARSseq.mmus.fin.metadata)))
#mmus.read.percentages <- rbind(mmus.read.percentages, data.frame(CELseq2.mmus.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("CELseq2", nrow(CELseq2.mmus.fin.metadata))))
#mmus.read.percentages <- rbind(mmus.read.percentages, data.frame(QUARTZseq.mmus.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("QUARTZseq", nrow(QUARTZseq.mmus.fin.metadata))))
#mmus.read.percentages <- rbind(mmus.read.percentages, data.frame(Dropseq.mmus.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("Dropseq", nrow(Dropseq.mmus.fin.metadata))))
#mmus.read.percentages <- rbind(mmus.read.percentages, data.frame(SCRBseq.mmus.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("SCRBseq", nrow(SCRBseq.mmus.fin.metadata))))
#mmus.read.percentages <- rbind(mmus.read.percentages, data.frame(SeqwellV1.mmus.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("SeqwellV2", nrow(SeqwellV2.mmus.fin.metadata))))
#mmus.read.percentages <- rbind(mmus.read.percentages, data.frame(SeqwellV2.mmus.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("SeqwellV2", nrow(SeqwellV2.mmus.fin.metadata))))
#mmus.read.percentages <- rbind(mmus.read.percentages, data.frame(Nuclei10X.mmus.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("Nuclei10X", nrow(Nuclei10X.mmus.fin.metadata))))
#mmus.read.percentages <- rbind(mmus.read.percentages, data.frame(ICELL8.mmus.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("ICELL8", nrow(ICELL8.mmus.fin.metadata))))
#mmus.read.percentages <- rbind(mmus.read.percentages, data.frame(ddSEQ.mmus.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("ddSEQ", nrow(ddSEQ.mmus.fin.metadata))))
#mmus.read.percentages <- rbind(mmus.read.percentages, data.frame(ddSEQexp1.mmus.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.mmus.fin.metadata))))
#mmus.read.percentages <- rbind(mmus.read.percentages, data.frame(C1HT.mmus.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("C1HT", nrow(C1HT.mmus.fin.metadata))))
#mmus.read.percentages <- rbind(mmus.read.percentages, data.frame(X10x8x10K.mmus.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("X10x8x10K", nrow(X10x8x10K.mmus.fin.metadata))))
#mmus.read.percentages <- rbind(mmus.read.percentages, data.frame(X10Scilife.mmus.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("X10Scilife", nrow(X10Scilife.mmus.fin.metadata))))
#mmus.melted.read.percentages <- melt(mmus.read.percentages, id = c("Library","tech"), measure = c("nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads"))

#pdf("all_techs_V2/mmus_read_distributions_all.pdf")
#ggplot(data=mmus.melted.read.percentages, aes(Library, value, fill= variable)) + geom_bar(stat = "identity", width = 0.5, position = "stack") + facet_grid(~ tech, scales = "free")+theme (axis.text.x = element_text(angle = 90, hjust = 1))
#ggplot(data=mmus.melted.read.percentages, aes(tech, value, fill= variable)) + geom_bar(stat = "identity", width = 0.5, position = "stack") + facet_grid(~ tech, scales = "free")+theme (axis.text.x = element_text(angle = 90, hjust = 1))
#dev.off()





# Total number of reads per tech per lib ====
print("Total number of reads per tech per lib")
mmus.total.reads <- data.frame(MARSseq.mmus.fin.metadata[, c("Library", "nTReads")], tech = rep("MARSseq", nrow(MARSseq.mmus.fin.metadata)))
mmus.total.reads <- rbind(mmus.total.reads, data.frame(CELseq2.mmus.fin.metadata[, c("Library", "nTReads")], tech = rep("CELseq2", nrow(CELseq2.mmus.fin.metadata))))
mmus.total.reads <- rbind(mmus.total.reads, data.frame(QUARTZseq.mmus.fin.metadata[, c("Library", "nTReads")], tech = rep("QUARTZseq", nrow(QUARTZseq.mmus.fin.metadata))))
mmus.total.reads <- rbind(mmus.total.reads, data.frame(Dropseq.mmus.fin.metadata[, c("Library", "nTReads")], tech = rep("Dropseq", nrow(Dropseq.mmus.fin.metadata))))
mmus.total.reads <- rbind(mmus.total.reads, data.frame(SCRBseq.mmus.fin.metadata[, c("Library", "nTReads")], tech = rep("SCRBseq", nrow(SCRBseq.mmus.fin.metadata))))
mmus.total.reads <- rbind(mmus.total.reads, data.frame(SeqwellV1.mmus.fin.metadata[, c("Library", "nTReads")], tech = rep("SeqwellV1", nrow(SeqwellV1.mmus.fin.metadata))))
mmus.total.reads <- rbind(mmus.total.reads, data.frame(SeqwellV2.mmus.fin.metadata[, c("Library", "nTReads")], tech = rep("SeqwellV2", nrow(SeqwellV2.mmus.fin.metadata))))
mmus.total.reads <- rbind(mmus.total.reads, data.frame(Nuclei10X.mmus.fin.metadata[, c("Library", "nTReads")], tech = rep("Nuclei10X", nrow(Nuclei10X.mmus.fin.metadata))))
mmus.total.reads <- rbind(mmus.total.reads, data.frame(ICELL8.mmus.fin.metadata[, c("Library", "nTReads")], tech = rep("ICELL8", nrow(ICELL8.mmus.fin.metadata))))
mmus.total.reads <- rbind(mmus.total.reads, data.frame(ddSEQ.mmus.fin.metadata[, c("Library", "nTReads")], tech = rep("ddSEQ", nrow(ddSEQ.mmus.fin.metadata))))
mmus.total.reads <- rbind(mmus.total.reads, data.frame(ddSEQexp1.mmus.fin.metadata[, c("Library", "nTReads")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.mmus.fin.metadata))))
mmus.total.reads <- rbind(mmus.total.reads, data.frame(C1HT.mmus.fin.metadata[, c("Library", "nTReads")], tech = rep("C1HT", nrow(C1HT.mmus.fin.metadata))))
mmus.total.reads <- rbind(mmus.total.reads, data.frame(X10x8x10K.mmus.fin.metadata[, c("Library", "nTReads")], tech = rep("X10x8x10K", nrow(X10x8x10K.mmus.fin.metadata))))
mmus.total.reads <- rbind(mmus.total.reads, data.frame(X10Scilife.mmus.fin.metadata[, c("Library", "nTReads")], tech = rep("X10Scilife", nrow(X10Scilife.mmus.fin.metadata))))
mmus.total.reads <- rbind(mmus.total.reads, data.frame(CB1.mmus.fin.metadata[, c("Library", "nTReads")], tech = rep("CB1", nrow(CB1.mmus.fin.metadata))))
pdf("all_techs_V2/mmus_total_number_Reads_all_V5.pdf")
#ggplot(data=mmus.total.reads, aes(x=Library, y=log(nTReads), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
#ggplot(data=mmus.total.reads, aes(x=tech, y=log(nTReads), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
ggplot(data=mmus.total.reads, aes(x=tech, y=nTReads, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(1000, 20000, 100000, 500000, 1000000)) + facet_grid(. ~ tech, scales = "free") #+ geom_hline(yintercept=c(1000, 20000, 100000, 500000, 1000000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
#ggplot(data=mmus.total.reads, aes(x=tech, y=nTReads, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(labels= scales::comma, breaks = c(1000, 20000, 100000, 500000, 1000000)) + facet_grid(. ~ tech, scales = "free") + geom_hline(yintercept=c(1000, 20000, 100000, 500000, 1000000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
dev.off()




#Total number of genes per tech per lib ====
print("Total number of genes per tech per lib ")
mmus.total.genes <- data.frame(MARSseq.mmus.fin.metadata[, c("Library", "nGenes")], tech = rep("MARSseq", nrow(MARSseq.mmus.fin.metadata)))
mmus.total.genes <- rbind(mmus.total.genes, data.frame(CELseq2.mmus.fin.metadata[, c("Library", "nGenes")], tech = rep("CELseq2", nrow(CELseq2.mmus.fin.metadata))))
mmus.total.genes <- rbind(mmus.total.genes, data.frame(QUARTZseq.mmus.fin.metadata[, c("Library", "nGenes")], tech = rep("QUARTZseq", nrow(QUARTZseq.mmus.fin.metadata))))
mmus.total.genes <- rbind(mmus.total.genes, data.frame(Dropseq.mmus.fin.metadata[, c("Library", "nGenes")], tech = rep("Dropseq", nrow(Dropseq.mmus.fin.metadata))))
mmus.total.genes <- rbind(mmus.total.genes, data.frame(SCRBseq.mmus.fin.metadata[, c("Library", "nGenes")], tech = rep("SCRBseq", nrow(SCRBseq.mmus.fin.metadata))))
mmus.total.genes <- rbind(mmus.total.genes, data.frame(SeqwellV1.mmus.fin.metadata[, c("Library", "nGenes")], tech = rep("SeqwellV1", nrow(SeqwellV1.mmus.fin.metadata))))
mmus.total.genes <- rbind(mmus.total.genes, data.frame(SeqwellV2.mmus.fin.metadata[, c("Library", "nGenes")], tech = rep("SeqwellV2", nrow(SeqwellV2.mmus.fin.metadata))))
mmus.total.genes <- rbind(mmus.total.genes, data.frame(Nuclei10X.mmus.fin.metadata[, c("Library", "nGenes")], tech = rep("Nuclei10X", nrow(Nuclei10X.mmus.fin.metadata))))
mmus.total.genes <- rbind(mmus.total.genes, data.frame(ICELL8.mmus.fin.metadata[, c("Library", "nGenes")], tech = rep("ICELL8", nrow(ICELL8.mmus.fin.metadata))))
mmus.total.genes <- rbind(mmus.total.genes, data.frame(ddSEQ.mmus.fin.metadata[, c("Library", "nGenes")], tech = rep("ddSEQ", nrow(ddSEQ.mmus.fin.metadata))))
mmus.total.genes <- rbind(mmus.total.genes, data.frame(ddSEQexp1.mmus.fin.metadata[, c("Library", "nGenes")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.mmus.fin.metadata))))
mmus.total.genes <- rbind(mmus.total.genes, data.frame(C1HT.mmus.fin.metadata[, c("Library", "nGenes")], tech = rep("C1HT", nrow(C1HT.mmus.fin.metadata))))
mmus.total.genes <- rbind(mmus.total.genes, data.frame(X10x8x10K.mmus.fin.metadata[, c("Library", "nGenes")], tech = rep("X10x8x10K", nrow(X10x8x10K.mmus.fin.metadata))))
mmus.total.genes <- rbind(mmus.total.genes, data.frame(X10Scilife.mmus.fin.metadata[, c("Library", "nGenes")], tech = rep("X10Scilife", nrow(X10Scilife.mmus.fin.metadata))))
mmus.total.genes <- rbind(mmus.total.genes, data.frame(CB1.mmus.fin.metadata[, c("Library", "nGenes")], tech = rep("CB1", nrow(CB1.mmus.fin.metadata))))
pdf("all_techs_V2/mmus_total_number_Genes_all_V5.pdf")
#ggplot(data=mmus.total.genes, aes(x=Library, y=log(nGenes), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
#ggplot(data=mmus.total.genes, aes(x=tech, y=log(nGenes), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
ggplot(data=mmus.total.genes, aes(x=tech, y=nGenes, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(10,150,500, 1500, 5000)) + facet_grid(. ~ tech, scales = "free") #+ geom_hline(yintercept=c(10,150,500, 1500, 5000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
#ggplot(data=mmus.total.genes, aes(x=tech, y=nGenes, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(labels= scales::comma, breaks = c(10,150,500, 1500, 5000)) + facet_grid(. ~ tech, scales = "free") + geom_hline(yintercept=c(10,150,500, 1500, 5000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
dev.off()


# Total number of UMIs per tech per lib ====
print("Total number of UMIs per tech per lib")
mmus.total.UMIs <- data.frame(MARSseq.mmus.fin.metadata[, c("Library", "nUMIs")], tech = rep("MARSseq", nrow(MARSseq.mmus.fin.metadata)))
mmus.total.UMIs <- rbind(mmus.total.UMIs, data.frame(CELseq2.mmus.fin.metadata[, c("Library", "nUMIs")], tech = rep("CELseq2", nrow(CELseq2.mmus.fin.metadata))))
mmus.total.UMIs <- rbind(mmus.total.UMIs, data.frame(QUARTZseq.mmus.fin.metadata[, c("Library", "nUMIs")], tech = rep("QUARTZseq", nrow(QUARTZseq.mmus.fin.metadata))))
mmus.total.UMIs <- rbind(mmus.total.UMIs, data.frame(Dropseq.mmus.fin.metadata[, c("Library", "nUMIs")], tech = rep("Dropseq", nrow(Dropseq.mmus.fin.metadata))))
mmus.total.UMIs <- rbind(mmus.total.UMIs, data.frame(SCRBseq.mmus.fin.metadata[, c("Library", "nUMIs")], tech = rep("SCRBseq", nrow(SCRBseq.mmus.fin.metadata))))
mmus.total.UMIs <- rbind(mmus.total.UMIs, data.frame(SeqwellV1.mmus.fin.metadata[, c("Library", "nUMIs")], tech = rep("SeqwellV1", nrow(SeqwellV1.mmus.fin.metadata))))
mmus.total.UMIs <- rbind(mmus.total.UMIs, data.frame(SeqwellV2.mmus.fin.metadata[, c("Library", "nUMIs")], tech = rep("SeqwellV2", nrow(SeqwellV2.mmus.fin.metadata))))
mmus.total.UMIs <- rbind(mmus.total.UMIs, data.frame(Nuclei10X.mmus.fin.metadata[, c("Library", "nUMIs")], tech = rep("Nuclei10X", nrow(Nuclei10X.mmus.fin.metadata))))
mmus.total.UMIs <- rbind(mmus.total.UMIs, data.frame(ICELL8.mmus.fin.metadata[, c("Library", "nUMIs")], tech = rep("ICELL8", nrow(ICELL8.mmus.fin.metadata))))
mmus.total.UMIs <- rbind(mmus.total.UMIs, data.frame(ddSEQ.mmus.fin.metadata[, c("Library", "nUMIs")], tech = rep("ddSEQ", nrow(ddSEQ.mmus.fin.metadata))))
mmus.total.UMIs <- rbind(mmus.total.UMIs, data.frame(ddSEQexp1.mmus.fin.metadata[, c("Library", "nUMIs")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.mmus.fin.metadata))))
mmus.total.UMIs <- rbind(mmus.total.UMIs, data.frame(C1HT.mmus.fin.metadata[, c("Library", "nUMIs")], tech = rep("C1HT", nrow(C1HT.mmus.fin.metadata))))
mmus.total.UMIs <- rbind(mmus.total.UMIs, data.frame(X10x8x10K.mmus.fin.metadata[, c("Library", "nUMIs")], tech = rep("X10x8x10K", nrow(X10x8x10K.mmus.fin.metadata))))
mmus.total.UMIs <- rbind(mmus.total.UMIs, data.frame(X10Scilife.mmus.fin.metadata[, c("Library", "nUMIs")], tech = rep("X10Scilife", nrow(X10Scilife.mmus.fin.metadata))))
mmus.total.UMIs <- rbind(mmus.total.UMIs, data.frame(CB1.mmus.fin.metadata[, c("Library", "nUMIs")], tech = rep("CB1", nrow(CB1.mmus.fin.metadata))))
pdf("all_techs_V2/mmus_total_number_UMIs_all_V5.pdf")
#ggplot(data=mmus.total.UMIs, aes(x=Library, y=log(nUMIs), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free")
#ggplot(data=mmus.total.UMIs, aes(x=tech, y=log(nUMIs), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free")
ggplot(data=mmus.total.UMIs, aes(x=tech, y=nUMIs, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(200,2000,20000,200000, 1000000)) + facet_grid(. ~ tech, scales = "free") #+ geom_hline(yintercept=c(200,2000,20000,200000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
#ggplot(data=mmus.total.UMIs, aes(x=tech, y=nUMIs, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(labels= scales::comma, breaks = c(200,2000,20000,200000, 1000000)) + facet_grid(. ~ tech, scales = "free")+ geom_hline(yintercept=c(200,2000,20000,200000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
dev.off()


# nTreads vs nUMI curves ====
print("nTreads vs nUMI curves")
mmus.nread.numi.df <- data.frame(MARSseq.mmus.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("MARSseq", nrow(MARSseq.mmus.fin.metadata)))
mmus.nread.numi.df <- rbind(mmus.nread.numi.df, data.frame(CELseq2.mmus.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("CELseq2", nrow(CELseq2.mmus.fin.metadata))))
mmus.nread.numi.df <- rbind(mmus.nread.numi.df, data.frame(QUARTZseq.mmus.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("QUARTZseq", nrow(QUARTZseq.mmus.fin.metadata))))
mmus.nread.numi.df <- rbind(mmus.nread.numi.df, data.frame(Dropseq.mmus.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("Dropseq", nrow(Dropseq.mmus.fin.metadata))))
mmus.nread.numi.df <- rbind(mmus.nread.numi.df, data.frame(SCRBseq.mmus.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("SCRBseq", nrow(SCRBseq.mmus.fin.metadata))))
mmus.nread.numi.df <- rbind(mmus.nread.numi.df, data.frame(SeqwellV1.mmus.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("SeqwellV1", nrow(SeqwellV1.mmus.fin.metadata))))
mmus.nread.numi.df <- rbind(mmus.nread.numi.df, data.frame(SeqwellV2.mmus.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("SeqwellV2", nrow(SeqwellV2.mmus.fin.metadata))))
mmus.nread.numi.df <- rbind(mmus.nread.numi.df, data.frame(Nuclei10X.mmus.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("Nuclei10X", nrow(Nuclei10X.mmus.fin.metadata))))
mmus.nread.numi.df <- rbind(mmus.nread.numi.df, data.frame(ICELL8.mmus.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("ICELL8", nrow(ICELL8.mmus.fin.metadata))))
mmus.nread.numi.df <- rbind(mmus.nread.numi.df, data.frame(ddSEQ.mmus.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("ddSEQ", nrow(ddSEQ.mmus.fin.metadata))))
mmus.nread.numi.df <- rbind(mmus.nread.numi.df, data.frame(ddSEQexp1.mmus.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.mmus.fin.metadata))))
mmus.nread.numi.df <- rbind(mmus.nread.numi.df, data.frame(C1HT.mmus.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("C1HT", nrow(C1HT.mmus.fin.metadata))))
mmus.nread.numi.df <- rbind(mmus.nread.numi.df, data.frame(X10x8x10K.mmus.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("X10x8x10K", nrow(X10x8x10K.mmus.fin.metadata))))
mmus.nread.numi.df <- rbind(mmus.nread.numi.df, data.frame(X10Scilife.mmus.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("X10Scilife", nrow(X10Scilife.mmus.fin.metadata))))
mmus.nread.numi.df <- rbind(mmus.nread.numi.df, data.frame(CB1.mmus.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("CB1", nrow(CB1.mmus.fin.metadata))))

pdf("all_techs_V2/mmus_nReadsVSnUMIs_all_V5.pdf")
#ggplot(mmus.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(aes(color=tech))+ xlim(0, 2000000)# + geom_line(aes(color=tech))#
ggplot(mmus.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(method="lm", formula=y~x, fill="grey",aes(color=tech))+ xlim(0, 100000)# + geom_line(aes(color=tech))#
#ggplot(mmus.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(method="lm", formula=y~log(x), fill="grey", aes(color=tech))+ xlim(0, 2000000)# + geom_line(aes(color=tech))#
#ggplot(mmus.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(method='nlsLM',formula=y ~ a*x^b, se=FALSE, method.args = list(start = list(a=1,b=1)), aes(color=tech))+ xlim(0, 1500000)# + geom_line(aes(color=tech))#
#ggplot(mmus.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(method='nls',formula=y ~ SSasympOff(x, A, lrc, c0), se=FALSE, aes(color=tech))+ xlim(0, 1500000)# + geom_line(aes(color=tech))#

dev.off()



# nTreads vs nGenes curves ====
print("nTreads vs nGenes curves")
mmus.nread.ngene.df <- data.frame(MARSseq.mmus.fin.metadata[, c("nTReads", "nGenes")], tech = rep("MARSseq", nrow(MARSseq.mmus.fin.metadata)))
mmus.nread.ngene.df <- rbind(mmus.nread.ngene.df, data.frame(CELseq2.mmus.fin.metadata[, c("nTReads", "nGenes")], tech = rep("CELseq2", nrow(CELseq2.mmus.fin.metadata))))
mmus.nread.ngene.df <- rbind(mmus.nread.ngene.df, data.frame(QUARTZseq.mmus.fin.metadata[, c("nTReads", "nGenes")], tech = rep("QUARTZseq", nrow(QUARTZseq.mmus.fin.metadata))))
mmus.nread.ngene.df <- rbind(mmus.nread.ngene.df, data.frame(Dropseq.mmus.fin.metadata[, c("nTReads", "nGenes")], tech = rep("Dropseq", nrow(Dropseq.mmus.fin.metadata))))
mmus.nread.ngene.df <- rbind(mmus.nread.ngene.df, data.frame(SCRBseq.mmus.fin.metadata[, c("nTReads", "nGenes")], tech = rep("SCRBseq", nrow(SCRBseq.mmus.fin.metadata))))
mmus.nread.ngene.df <- rbind(mmus.nread.ngene.df, data.frame(SeqwellV1.mmus.fin.metadata[, c("nTReads", "nGenes")], tech = rep("SeqwellV1", nrow(SeqwellV1.mmus.fin.metadata))))
mmus.nread.ngene.df <- rbind(mmus.nread.ngene.df, data.frame(SeqwellV2.mmus.fin.metadata[, c("nTReads", "nGenes")], tech = rep("SeqwellV2", nrow(SeqwellV2.mmus.fin.metadata))))
mmus.nread.ngene.df <- rbind(mmus.nread.ngene.df, data.frame(Nuclei10X.mmus.fin.metadata[, c("nTReads", "nGenes")], tech = rep("Nuclei10X", nrow(Nuclei10X.mmus.fin.metadata))))
mmus.nread.ngene.df <- rbind(mmus.nread.ngene.df, data.frame(ICELL8.mmus.fin.metadata[, c("nTReads", "nGenes")], tech = rep("ICELL8", nrow(ICELL8.mmus.fin.metadata))))
mmus.nread.ngene.df <- rbind(mmus.nread.ngene.df, data.frame(ddSEQ.mmus.fin.metadata[, c("nTReads", "nGenes")], tech = rep("ddSEQ", nrow(ddSEQ.mmus.fin.metadata))))
mmus.nread.ngene.df <- rbind(mmus.nread.ngene.df, data.frame(ddSEQexp1.mmus.fin.metadata[, c("nTReads", "nGenes")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.mmus.fin.metadata))))
mmus.nread.ngene.df <- rbind(mmus.nread.ngene.df, data.frame(C1HT.mmus.fin.metadata[, c("nTReads", "nGenes")], tech = rep("C1HT", nrow(C1HT.mmus.fin.metadata))))
mmus.nread.ngene.df <- rbind(mmus.nread.ngene.df, data.frame(X10x8x10K.mmus.fin.metadata[, c("nTReads", "nGenes")], tech = rep("X10x8x10K", nrow(X10x8x10K.mmus.fin.metadata))))
mmus.nread.ngene.df <- rbind(mmus.nread.ngene.df, data.frame(X10Scilife.mmus.fin.metadata[, c("nTReads", "nGenes")], tech = rep("X10Scilife", nrow(X10Scilife.mmus.fin.metadata))))
mmus.nread.ngene.df <- rbind(mmus.nread.ngene.df, data.frame(CB1.mmus.fin.metadata[, c("nTReads", "nGenes")], tech = rep("CB1", nrow(CB1.mmus.fin.metadata))))
pdf("all_techs_V2/mmus_nReadsVSnGenes_all_V5.pdf")
#ggplot(mmus.nread.ngene.df, aes(x=nGenes, y=nTReads, group=tech)) + geom_point(aes(color=tech)) + geom_smooth(aes(color=tech))+ xlim(0, 18000)# + geom_line(aes(color=tech))#
ggplot(mmus.nread.ngene.df, aes(x=nTReads, y=nGenes, group=tech))  + geom_smooth(method="lm", formula=y~x, fill="grey",aes(color=tech))+ xlim(0, 100000)# + geom_line(aes(color=tech))#
#ggplot(mmus.nread.ngene.df, aes(x=nTReads, y=nGenes, group=tech)) + geom_smooth(method="lm", formula=y~log(x), fill="grey", aes(color=tech))+ xlim(0, 2000000)# + geom_line(aes(color=tech))#
#ggplot(mmus.nread.ngene.df, aes(x=nTReads, y=nGenes, group=tech)) + geom_smooth(method='nlsLM',formula=y ~ a*x^b, se=FALSE, method.args = list(start = list(a=1,b=1)), aes(color=tech))+ xlim(0, 1500000)# + geom_line(aes(color=tech))#
#ggplot(mmus.nread.ngene.df, aes(x=nTReads, y=nGenes, group=tech)) + geom_smooth(method='nls',formula=y ~ SSasympOff(x, A, lrc, c0), se=FALSE, aes(color=tech))+ xlim(0, 1500000)# + geom_line(aes(color=tech))#
dev.off()