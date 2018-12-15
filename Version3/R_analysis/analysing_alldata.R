library(scater)
library(ggplot2)
library(biomaRt)
library(tibble)
library(data.table)


# Functions ####
load("Biomart_hsap_mapping_table.RData")
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


# MARSseq hsap Preparation ====
print("MARSseq is Running...")
load("../SCE_Robjects/MARSseq.hsap.full.SCE.Robj")
MARSseq.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

MARSseq.hsap.count <- as.data.frame(counts(MARSseq.hsap))
MARSseq.hsap.count.mapped <- mapIDs(MARSseq.hsap.count, "hsap")
MARSseq.hsap.metadata <- as.data.frame(colData(MARSseq.hsap))
MARSseq.hsap.metadata[is.na(MARSseq.hsap.metadata)] <- 0
MARSseq.hsap.tmappedReads <- rowSums(MARSseq.hsap.metadata[, c(2:4,6,7)])
MARSseq.hsap.mappedReadsPerc <- (MARSseq.hsap.tmappedReads/MARSseq.hsap.metadata$nTReads)*100
MARSseq.hsap.metadata <- add_column(MARSseq.hsap.metadata, mappedReads= MARSseq.hsap.tmappedReads, mappedReadsPerc= MARSseq.hsap.mappedReadsPerc, .after = 1 )
MARSseq.hsap.cells <- which(MARSseq.hsap.metadata$Species == "Human")
MARSseq.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(MARSseq.hsap.count.mapped[,MARSseq.hsap.cells])), colData = MARSseq.hsap.metadata[MARSseq.hsap.cells,c(1:15)])
rm(MARSseq.hsap.count)
rm(MARSseq.hsap.count.mapped)
rm(MARSseq.hsap.metadata)
MARSseq.hsap.is.mito <- grep("^MT-", rownames(counts(MARSseq.hsap.new)))
MARSseq.hsap.new <- calculateQCMetrics(MARSseq.hsap.new, feature_controls=list(Mt=MARSseq.hsap.is.mito))
MARSseq.hsap.nTReads.drop <- which(isOutlier(MARSseq.hsap.new$mappedReads, nmads=3, type="lower", log=TRUE))
MARSseq.hsap.mappedPct.drop <- which(MARSseq.hsap.new$mappedReadsPerc < 65)
MARSseq.hsap.mito.drop <- which(isOutlier(MARSseq.hsap.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(MARSseq.hsap.new)[-MARSseq.hsap.mito.drop, "log10_total_counts_Mt"]))
MARSseq.hsap.final <- MARSseq.hsap.new[,-c(MARSseq.hsap.nTReads.drop, MARSseq.hsap.mappedPct.drop,MARSseq.hsap.mito.drop)]
MARSseq.hsap.fin.metadata <- colData(MARSseq.hsap.final)
#### end ####

# CELseq2 hsap Preparation ====
print("CELseq2 is Running...")
load("../SCE_Robjects/CELseq2.hsap.full.SCE.Robj")
CELseq2.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

CELseq2.hsap.count <- as.data.frame(counts(CELseq2.hsap))
CELseq2.hsap.count.mapped <- mapIDs(CELseq2.hsap.count, "hsap")
CELseq2.hsap.metadata <- as.data.frame(colData(CELseq2.hsap))
CELseq2.hsap.metadata[is.na(CELseq2.hsap.metadata)] <- 0
CELseq2.hsap.tmappedReads <- rowSums(CELseq2.hsap.metadata[, c(2:4,6,7)])
CELseq2.hsap.mappedReadsPerc <- (CELseq2.hsap.tmappedReads/CELseq2.hsap.metadata$nTReads)*100
CELseq2.hsap.metadata <- add_column(CELseq2.hsap.metadata, mappedReads= CELseq2.hsap.tmappedReads, mappedReadsPerc= CELseq2.hsap.mappedReadsPerc, .after = 1 )
CELseq2.hsap.cells <- which(CELseq2.hsap.metadata$Species == "Human")
CELseq2.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(CELseq2.hsap.count.mapped[,CELseq2.hsap.cells])), colData = CELseq2.hsap.metadata[CELseq2.hsap.cells,c(1:15)])
rm(CELseq2.hsap.count)
rm(CELseq2.hsap.count.mapped)
#rm(CELseq2.hsap.metadata)
CELseq2.hsap.is.mito <- grep("^MT-", rownames(counts(CELseq2.hsap.new)))
CELseq2.hsap.new <- calculateQCMetrics(CELseq2.hsap.new, feature_controls=list(Mt=CELseq2.hsap.is.mito))
CELseq2.hsap.nTReads.drop <- which(isOutlier(CELseq2.hsap.new$mappedReads, nmads=3, type="lower", log=TRUE))
CELseq2.hsap.mappedPct.drop <- which(CELseq2.hsap.new$mappedReadsPerc < 65)
CELseq2.hsap.mito.drop <- which(isOutlier(CELseq2.hsap.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(CELseq2.hsap.new)[-CELseq2.hsap.mito.drop, "log10_total_counts_Mt"]))
CELseq2.hsap.final <- CELseq2.hsap.new[,-c(CELseq2.hsap.nTReads.drop, CELseq2.hsap.mappedPct.drop,CELseq2.hsap.mito.drop)]
CELseq2.hsap.fin.metadata <- colData(CELseq2.hsap.final)
#### end ####

# QUARTZseq hsap Preparation ====
print("QUARTZseq is Running...")
load("../SCE_Robjects/QUARTZseq.hsap.full.SCE.Robj")
QUARTZseq.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

QUARTZseq.hsap.count <- as.data.frame(counts(QUARTZseq.hsap))
QUARTZseq.hsap.count.mapped <- mapIDs(QUARTZseq.hsap.count, "hsap")
QUARTZseq.hsap.metadata <- as.data.frame(colData(QUARTZseq.hsap))
QUARTZseq.hsap.metadata[is.na(QUARTZseq.hsap.metadata)] <- 0
QUARTZseq.hsap.tmappedReads <- rowSums(QUARTZseq.hsap.metadata[, c(2:4,6,7)])
QUARTZseq.hsap.mappedReadsPerc <- (QUARTZseq.hsap.tmappedReads/QUARTZseq.hsap.metadata$nTReads)*100
QUARTZseq.hsap.metadata <- add_column(QUARTZseq.hsap.metadata, mappedReads= QUARTZseq.hsap.tmappedReads, mappedReadsPerc= QUARTZseq.hsap.mappedReadsPerc, .after = 1 )
QUARTZseq.hsap.cells <- which(QUARTZseq.hsap.metadata$Species == "Human")
QUARTZseq.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(QUARTZseq.hsap.count.mapped[,QUARTZseq.hsap.cells])), colData = QUARTZseq.hsap.metadata[QUARTZseq.hsap.cells,c(1:15)])
rm(QUARTZseq.hsap.count)
rm(QUARTZseq.hsap.count.mapped)
#rm(QUARTZseq.hsap.metadata)
QUARTZseq.hsap.is.mito <- grep("^MT-", rownames(counts(QUARTZseq.hsap.new)))
QUARTZseq.hsap.new <- calculateQCMetrics(QUARTZseq.hsap.new, feature_controls=list(Mt=QUARTZseq.hsap.is.mito))
QUARTZseq.hsap.nTReads.drop <- which(isOutlier(QUARTZseq.hsap.new$mappedReads, nmads=3, type="lower", log=TRUE))
QUARTZseq.hsap.mappedPct.drop <- which(QUARTZseq.hsap.new$mappedReadsPerc < 65)
QUARTZseq.hsap.mito.drop <- which(isOutlier(QUARTZseq.hsap.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(QUARTZseq.hsap.new)[-QUARTZseq.hsap.mito.drop, "log10_total_counts_Mt"]))
QUARTZseq.hsap.final <- QUARTZseq.hsap.new[,-c(QUARTZseq.hsap.nTReads.drop, QUARTZseq.hsap.mappedPct.drop,QUARTZseq.hsap.mito.drop)]
QUARTZseq.hsap.fin.metadata <- colData(QUARTZseq.hsap.final)
#### end ####

# Dropseq hsap Preparation ====
print("Dropseq is Running...")
load("../SCE_Robjects/Dropseq.hsap.full.SCE.1000cellsPerPool.Robj")
Dropseq.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

Dropseq.hsap.count <- as.data.frame(counts(Dropseq.hsap))
Dropseq.hsap.count.mapped <- mapIDs(Dropseq.hsap.count, "hsap")
Dropseq.hsap.metadata <- as.data.frame(colData(Dropseq.hsap))
Dropseq.hsap.metadata[is.na(Dropseq.hsap.metadata)] <- 0
Dropseq.hsap.tmappedReads <- rowSums(Dropseq.hsap.metadata[, c(2:4,6,7)])
Dropseq.hsap.mappedReadsPerc <- (Dropseq.hsap.tmappedReads/Dropseq.hsap.metadata$nTReads)*100
Dropseq.hsap.metadata <- add_column(Dropseq.hsap.metadata, mappedReads= Dropseq.hsap.tmappedReads, mappedReadsPerc= Dropseq.hsap.mappedReadsPerc, .after = 1 )
Dropseq.hsap.cells <- which(Dropseq.hsap.metadata$Species == "Human")
Dropseq.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(Dropseq.hsap.count.mapped[,Dropseq.hsap.cells])), colData = Dropseq.hsap.metadata[Dropseq.hsap.cells,c(1:15)])
rm(Dropseq.hsap.count)
rm(Dropseq.hsap.count.mapped)

Dropseq.hsap.is.mito <- grep("^MT-", rownames(counts(Dropseq.hsap.new)))
Dropseq.hsap.new <- calculateQCMetrics(Dropseq.hsap.new, feature_controls=list(Mt=Dropseq.hsap.is.mito))
Dropseq.hsap.nTReads.drop <- which(isOutlier(Dropseq.hsap.new$mappedReads, nmads=3, type="lower", log=TRUE))
Dropseq.hsap.mappedPct.drop <- which(Dropseq.hsap.new$mappedReadsPerc < 65)
Dropseq.hsap.mito.drop <- which(isOutlier(Dropseq.hsap.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(Dropseq.hsap.new)[-Dropseq.hsap.mito.drop, "log10_total_counts_Mt"]))
Dropseq.hsap.final <- Dropseq.hsap.new[,-c(Dropseq.hsap.nTReads.drop, Dropseq.hsap.mappedPct.drop,Dropseq.hsap.mito.drop)]
Dropseq.hsap.fin.metadata <- colData(Dropseq.hsap.final)
#### end ####

# SCRBseq hsap Preparation ====
print("SCRBseq is Running...")
load("../SCE_Robjects/SCRBseq.hsap.full.SCE.Robj")
SCRBseq.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

SCRBseq.hsap.count <- as.data.frame(counts(SCRBseq.hsap))
SCRBseq.hsap.count.mapped <- mapIDs(SCRBseq.hsap.count, "hsap")
SCRBseq.hsap.metadata <- as.data.frame(colData(SCRBseq.hsap))
SCRBseq.hsap.metadata[is.na(SCRBseq.hsap.metadata)] <- 0
SCRBseq.hsap.tmappedReads <- rowSums(SCRBseq.hsap.metadata[, c(2:4,6,7)])
SCRBseq.hsap.mappedReadsPerc <- (SCRBseq.hsap.tmappedReads/SCRBseq.hsap.metadata$nTReads)*100
SCRBseq.hsap.metadata <- add_column(SCRBseq.hsap.metadata, mappedReads= SCRBseq.hsap.tmappedReads, mappedReadsPerc= SCRBseq.hsap.mappedReadsPerc, .after = 1 )
SCRBseq.hsap.cells <- which(SCRBseq.hsap.metadata$Species == "Human")
SCRBseq.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(SCRBseq.hsap.count.mapped[,SCRBseq.hsap.cells])), colData = SCRBseq.hsap.metadata[SCRBseq.hsap.cells,c(1:15)])

SCRBseq.hsap.is.mito <- grep("^MT-", rownames(counts(SCRBseq.hsap.new)))
SCRBseq.hsap.new <- calculateQCMetrics(SCRBseq.hsap.new, feature_controls=list(Mt=SCRBseq.hsap.is.mito))
SCRBseq.hsap.nTReads.drop <- which(isOutlier(SCRBseq.hsap.new$mappedReads, nmads=3, type="lower", log=TRUE))
SCRBseq.hsap.mappedPct.drop <- which(SCRBseq.hsap.new$mappedReadsPerc < 65)
SCRBseq.hsap.mito.drop <- which(isOutlier(SCRBseq.hsap.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(SCRBseq.hsap.new)[-SCRBseq.hsap.mito.drop, "log10_total_counts_Mt"]))
SCRBseq.hsap.final <- SCRBseq.hsap.new[,-c(SCRBseq.hsap.nTReads.drop, SCRBseq.hsap.mappedPct.drop,SCRBseq.hsap.mito.drop)]
SCRBseq.hsap.fin.metadata <- colData(SCRBseq.hsap.final)
#### end ####


# SeqwellV2 hsap Preparation ====
print("SeqwellV2 is Running...")
load("../SCE_Robjects/SeqWellV2.hsap.full.SCE.Robj")
SeqwellV2.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

SeqwellV2.hsap.count <- as.data.frame(counts(SeqwellV2.hsap))
SeqwellV2.hsap.count.mapped <- mapIDs(SeqwellV2.hsap.count, "hsap")
SeqwellV2.hsap.metadata <- as.data.frame(colData(SeqwellV2.hsap))
SeqwellV2.hsap.metadata[is.na(SeqwellV2.hsap.metadata)] <- 0
SeqwellV2.hsap.tmappedReads <- rowSums(SeqwellV2.hsap.metadata[, c(2:4,6,7)])
SeqwellV2.hsap.mappedReadsPerc <- (SeqwellV2.hsap.tmappedReads/SeqwellV2.hsap.metadata$nTReads)*100
SeqwellV2.hsap.metadata <- add_column(SeqwellV2.hsap.metadata, mappedReads= SeqwellV2.hsap.tmappedReads, mappedReadsPerc= SeqwellV2.hsap.mappedReadsPerc, .after = 1 )
SeqwellV2.hsap.cells <- which(SeqwellV2.hsap.metadata$Species == "Human")
SeqwellV2.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(SeqwellV2.hsap.count.mapped[,SeqwellV2.hsap.cells])), colData = SeqwellV2.hsap.metadata[SeqwellV2.hsap.cells,c(1:15)])

SeqwellV2.hsap.is.mito <- grep("^MT-", rownames(counts(SeqwellV2.hsap.new)))
SeqwellV2.hsap.new <- calculateQCMetrics(SeqwellV2.hsap.new, feature_controls=list(Mt=SeqwellV2.hsap.is.mito))
SeqwellV2.hsap.nTReads.drop <- which(isOutlier(SeqwellV2.hsap.new$mappedReads, nmads=3, type="lower", log=TRUE))
SeqwellV2.hsap.mappedPct.drop <- which(SeqwellV2.hsap.new$mappedReadsPerc < 65)
SeqwellV2.hsap.mito.drop <- which(isOutlier(SeqwellV2.hsap.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(SeqwellV2.hsap.new)[-SeqwellV2.hsap.mito.drop, "log10_total_counts_Mt"]))
SeqwellV2.hsap.final <- SeqwellV2.hsap.new[,-c(SeqwellV2.hsap.nTReads.drop, SeqwellV2.hsap.mappedPct.drop,SeqwellV2.hsap.mito.drop)]
SeqwellV2.hsap.fin.metadata <- colData(SeqwellV2.hsap.final)
#### end ####


# SeqwellV1 hsap Preparation ====
print("SeqwellV1 is Running...")
load("../SCE_Robjects/SeqWellV1.hsap.full.SCE.Robj")
SeqwellV1.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

SeqwellV1.hsap.count <- as.data.frame(counts(SeqwellV1.hsap))
SeqwellV1.hsap.count.mapped <- mapIDs(SeqwellV1.hsap.count, "hsap")
SeqwellV1.hsap.metadata <- as.data.frame(colData(SeqwellV1.hsap))
SeqwellV1.hsap.metadata[is.na(SeqwellV1.hsap.metadata)] <- 0
SeqwellV1.hsap.tmappedReads <- rowSums(SeqwellV1.hsap.metadata[, c(2:4,6)])
SeqwellV1.hsap.mappedReadsPerc <- (SeqwellV1.hsap.tmappedReads/SeqwellV1.hsap.metadata$nTReads)*100
SeqwellV1.hsap.metadata <- add_column(SeqwellV1.hsap.metadata, mappedReads= SeqwellV1.hsap.tmappedReads, mappedReadsPerc= SeqwellV1.hsap.mappedReadsPerc, .after = 1 )
SeqwellV1.hsap.cells <- which(SeqwellV1.hsap.metadata$Species == "Human")
SeqwellV1.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(SeqwellV1.hsap.count.mapped[,SeqwellV1.hsap.cells])), colData = SeqwellV1.hsap.metadata[SeqwellV1.hsap.cells,c(1:15)])

SeqwellV1.hsap.is.mito <- grep("^MT-", rownames(counts(SeqwellV1.hsap.new)))
SeqwellV1.hsap.new <- calculateQCMetrics(SeqwellV1.hsap.new, feature_controls=list(Mt=SeqwellV1.hsap.is.mito))
SeqwellV1.hsap.nTReads.drop <- which(isOutlier(SeqwellV1.hsap.new$mappedReads, nmads=2, type="lower", log=TRUE))
SeqwellV1.hsap.mappedPct.drop <- which(SeqwellV1.hsap.new$mappedReadsPerc < 65)
SeqwellV1.hsap.mito.drop <- which(isOutlier(SeqwellV1.hsap.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(SeqwellV1.hsap.new)[-SeqwellV1.hsap.mito.drop, "log10_total_counts_Mt"]))
SeqwellV1.hsap.final <- SeqwellV1.hsap.new[,-c(SeqwellV1.hsap.nTReads.drop, SeqwellV1.hsap.mappedPct.drop,SeqwellV1.hsap.mito.drop)]
SeqwellV1.hsap.fin.metadata <- colData(SeqwellV1.hsap.final)
#### end ####


# Nuclei10X hsap Preparation ====
print("Nuclei10X is Running...")
load("../SCE_Robjects/Nuclei10X.hsap.full.SCE.Robj")
Nuclei10X.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

Nuclei10X.hsap.count <- as.data.frame(counts(Nuclei10X.hsap))
Nuclei10X.hsap.count.mapped <- mapIDs(Nuclei10X.hsap.count, "hsap")
Nuclei10X.hsap.metadata <- as.data.frame(colData(Nuclei10X.hsap))
Nuclei10X.hsap.metadata[is.na(Nuclei10X.hsap.metadata)] <- 0
Nuclei10X.hsap.tmappedReads <- rowSums(Nuclei10X.hsap.metadata[, c(2:4,6)])
Nuclei10X.hsap.mappedReadsPerc <- (Nuclei10X.hsap.tmappedReads/Nuclei10X.hsap.metadata$nTReads)*100
Nuclei10X.hsap.metadata <- add_column(Nuclei10X.hsap.metadata, mappedReads= Nuclei10X.hsap.tmappedReads, mappedReadsPerc= Nuclei10X.hsap.mappedReadsPerc, .after = 1 )
Nuclei10X.hsap.cells <- which(Nuclei10X.hsap.metadata$Species == "Human")
Nuclei10X.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(Nuclei10X.hsap.count.mapped[,Nuclei10X.hsap.cells])), colData = Nuclei10X.hsap.metadata[Nuclei10X.hsap.cells,c(1:15)])

Nuclei10X.hsap.is.mito <- grep("^MT-", rownames(counts(Nuclei10X.hsap.new)))
Nuclei10X.hsap.new <- calculateQCMetrics(Nuclei10X.hsap.new, feature_controls=list(Mt=Nuclei10X.hsap.is.mito))
Nuclei10X.hsap.nTReads.drop <- which(isOutlier(Nuclei10X.hsap.new$mappedReads, nmads=2, type="lower", log=TRUE))
Nuclei10X.hsap.mappedPct.drop <- which(Nuclei10X.hsap.new$mappedReadsPerc < 65)
Nuclei10X.hsap.mito.drop <- which(isOutlier(Nuclei10X.hsap.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(Nuclei10X.hsap.new)[-Nuclei10X.hsap.mito.drop, "log10_total_counts_Mt"]))
Nuclei10X.hsap.final <- Nuclei10X.hsap.new[,-c(Nuclei10X.hsap.nTReads.drop, Nuclei10X.hsap.mappedPct.drop,Nuclei10X.hsap.mito.drop)]
Nuclei10X.hsap.fin.metadata <- colData(Nuclei10X.hsap.final)
#### end ####


# ICELL8 hsap Preparation ====
print("ICELL8 is Running...")
load("../SCE_Robjects/ICELL8.hsap.full.SCE.Robj")
ICELL8.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

ICELL8.hsap.count <- as.data.frame(counts(ICELL8.hsap))
ICELL8.hsap.count.mapped <- mapIDs(ICELL8.hsap.count, "hsap")
ICELL8.hsap.metadata <- as.data.frame(colData(ICELL8.hsap))
ICELL8.hsap.metadata[is.na(ICELL8.hsap.metadata)] <- 0
ICELL8.hsap.tmappedReads <- rowSums(ICELL8.hsap.metadata[, c(2:4,6)])
ICELL8.hsap.mappedReadsPerc <- (ICELL8.hsap.tmappedReads/ICELL8.hsap.metadata$nTReads)*100
ICELL8.hsap.metadata <- add_column(ICELL8.hsap.metadata, mappedReads= ICELL8.hsap.tmappedReads, mappedReadsPerc= ICELL8.hsap.mappedReadsPerc, .after = 1 )
ICELL8.hsap.cells <- which(ICELL8.hsap.metadata$Species == "Human")
ICELL8.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(ICELL8.hsap.count.mapped[,ICELL8.hsap.cells])), colData = ICELL8.hsap.metadata[ICELL8.hsap.cells,c(1:15)])

ICELL8.hsap.is.mito <- grep("^MT-", rownames(counts(ICELL8.hsap.new)))
ICELL8.hsap.new <- calculateQCMetrics(ICELL8.hsap.new, feature_controls=list(Mt=ICELL8.hsap.is.mito))
ICELL8.hsap.nTReads.drop <- which(isOutlier(ICELL8.hsap.new$mappedReads, nmads=2, type="lower", log=TRUE))
ICELL8.hsap.mappedPct.drop <- which(ICELL8.hsap.new$mappedReadsPerc < 65)
ICELL8.hsap.mito.drop <- which(isOutlier(ICELL8.hsap.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(ICELL8.hsap.new)[-ICELL8.hsap.mito.drop, "log10_total_counts_Mt"]))
ICELL8.hsap.final <- ICELL8.hsap.new[,-c(ICELL8.hsap.nTReads.drop, ICELL8.hsap.mappedPct.drop,ICELL8.hsap.mito.drop)]
ICELL8.hsap.fin.metadata <- colData(ICELL8.hsap.final)
                                                                                                       
#### end ####


# ddSEQ hsap Preparation ====
print("ddSEQ is Running...")
load("../SCE_Robjects/ddSEQ.hsap.full.SCE.Robj")
ddSEQ.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

ddSEQ.hsap.count <- as.data.frame(counts(ddSEQ.hsap))
ddSEQ.hsap.count.mapped <- mapIDs(ddSEQ.hsap.count, "hsap")
ddSEQ.hsap.metadata <- as.data.frame(colData(ddSEQ.hsap))
ddSEQ.hsap.metadata[is.na(ddSEQ.hsap.metadata)] <- 0
ddSEQ.hsap.tmappedReads <- rowSums(ddSEQ.hsap.metadata[, c(2:4,6)])
ddSEQ.hsap.mappedReadsPerc <- (ddSEQ.hsap.tmappedReads/ddSEQ.hsap.metadata$nTReads)*100
ddSEQ.hsap.metadata <- add_column(ddSEQ.hsap.metadata, mappedReads= ddSEQ.hsap.tmappedReads, mappedReadsPerc= ddSEQ.hsap.mappedReadsPerc, .after = 1 )
ddSEQ.hsap.cells <- which(ddSEQ.hsap.metadata$Species == "Human")
ddSEQ.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(ddSEQ.hsap.count.mapped[,ddSEQ.hsap.cells])), colData = ddSEQ.hsap.metadata[ddSEQ.hsap.cells,c(1:15)])

ddSEQ.hsap.is.mito <- grep("^MT-", rownames(counts(ddSEQ.hsap.new)))
ddSEQ.hsap.new <- calculateQCMetrics(ddSEQ.hsap.new, feature_controls=list(Mt=ddSEQ.hsap.is.mito))
ddSEQ.hsap.nTReads.drop <- which(isOutlier(ddSEQ.hsap.new$mappedReads, nmads=2, type="lower", log=TRUE))
ddSEQ.hsap.mappedPct.drop <- which(ddSEQ.hsap.new$mappedReadsPerc < 65)
ddSEQ.hsap.mito.drop <- which(isOutlier(ddSEQ.hsap.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(ddSEQ.hsap.new)[-ddSEQ.hsap.mito.drop, "log10_total_counts_Mt"]))
ddSEQ.hsap.final <- ddSEQ.hsap.new[,-c(ddSEQ.hsap.nTReads.drop, ddSEQ.hsap.mappedPct.drop,ddSEQ.hsap.mito.drop)]
ddSEQ.hsap.fin.metadata <- colData(ddSEQ.hsap.final)
#### end ####

# ddSEQexp1 hsap Preparation ====
print("ddSEQexp1 is Running...")
load("../SCE_Robjects/ddSEQexp1.hsap.full.SCE.Robj")
ddSEQexp1.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

ddSEQexp1.hsap.count <- as.data.frame(counts(ddSEQexp1.hsap))
ddSEQexp1.hsap.count.mapped <- mapIDs(ddSEQexp1.hsap.count, "hsap")
ddSEQexp1.hsap.metadata <- as.data.frame(colData(ddSEQexp1.hsap))
ddSEQexp1.hsap.metadata[is.na(ddSEQexp1.hsap.metadata)] <- 0
ddSEQexp1.hsap.tmappedReads <- rowSums(ddSEQexp1.hsap.metadata[, c(2:4,6)])
ddSEQexp1.hsap.mappedReadsPerc <- (ddSEQexp1.hsap.tmappedReads/ddSEQexp1.hsap.metadata$nTReads)*100
ddSEQexp1.hsap.metadata <- add_column(ddSEQexp1.hsap.metadata, mappedReads= ddSEQexp1.hsap.tmappedReads, mappedReadsPerc= ddSEQexp1.hsap.mappedReadsPerc, .after = 1 )
ddSEQexp1.hsap.cells <- which(ddSEQexp1.hsap.metadata$Species == "Human")
ddSEQexp1.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(ddSEQexp1.hsap.count.mapped[,ddSEQexp1.hsap.cells])), colData = ddSEQexp1.hsap.metadata[ddSEQexp1.hsap.cells,c(1:15)])

ddSEQexp1.hsap.is.mito <- grep("^MT-", rownames(counts(ddSEQexp1.hsap.new)))
ddSEQexp1.hsap.new <- calculateQCMetrics(ddSEQexp1.hsap.new, feature_controls=list(Mt=ddSEQexp1.hsap.is.mito))
ddSEQexp1.hsap.nTReads.drop <- which(isOutlier(ddSEQexp1.hsap.new$mappedReads, nmads=2, type="lower", log=TRUE))
ddSEQexp1.hsap.mappedPct.drop <- which(ddSEQexp1.hsap.new$mappedReadsPerc < 65)
ddSEQexp1.hsap.mito.drop <- which(isOutlier(ddSEQexp1.hsap.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(ddSEQexp1.hsap.new)[-ddSEQexp1.hsap.mito.drop, "log10_total_counts_Mt"]))
ddSEQexp1.hsap.final <- ddSEQexp1.hsap.new[,-c(ddSEQexp1.hsap.nTReads.drop, ddSEQexp1.hsap.mappedPct.drop,ddSEQexp1.hsap.mito.drop)]
ddSEQexp1.hsap.fin.metadata <- colData(ddSEQexp1.hsap.final)
#### end ####


# C1HT hsap Preparation ====
print("C1HT is Running...")
load("../SCE_Robjects/C1HT.hsap.full.SCE.Robj")
C1HT.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

C1HT.hsap.count <- as.data.frame(counts(C1HT.hsap))
C1HT.hsap.count.mapped <- mapIDs(C1HT.hsap.count, "hsap")
C1HT.hsap.metadata <- as.data.frame(colData(C1HT.hsap))
C1HT.hsap.metadata[is.na(C1HT.hsap.metadata)] <- 0
C1HT.hsap.tmappedReads <- rowSums(C1HT.hsap.metadata[, c(2:4,6)])
C1HT.hsap.mappedReadsPerc <- (C1HT.hsap.tmappedReads/C1HT.hsap.metadata$nTReads)*100
C1HT.hsap.metadata <- add_column(C1HT.hsap.metadata, mappedReads= C1HT.hsap.tmappedReads, mappedReadsPerc= C1HT.hsap.mappedReadsPerc, .after = 1 )
C1HT.hsap.cells <- which(C1HT.hsap.metadata$Species == "Human")
C1HT.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(C1HT.hsap.count.mapped[,C1HT.hsap.cells])), colData = C1HT.hsap.metadata[C1HT.hsap.cells,c(1:15)])

C1HT.hsap.is.mito <- grep("^MT-", rownames(counts(C1HT.hsap.new)))
C1HT.hsap.new <- calculateQCMetrics(C1HT.hsap.new, feature_controls=list(Mt=C1HT.hsap.is.mito))
C1HT.hsap.nTReads.drop <- which(isOutlier(C1HT.hsap.new$mappedReads, nmads=2, type="lower", log=TRUE))
C1HT.hsap.mappedPct.drop <- which(C1HT.hsap.new$mappedReadsPerc < 65)
C1HT.hsap.mito.drop <- which(isOutlier(C1HT.hsap.new$log10_total_counts_Mt, nmads=2, type="higher"))
#plot(density(colData(C1HT.hsap.new)[-C1HT.hsap.mito.drop, "log10_total_counts_Mt"]))
C1HT.hsap.final <- C1HT.hsap.new[,-c(C1HT.hsap.nTReads.drop, C1HT.hsap.mappedPct.drop,C1HT.hsap.mito.drop)]
C1HT.hsap.fin.metadata <- colData(C1HT.hsap.final)

#### end ####



#10X8x10K hsap Preparation ====
print("10X8x10K is Running...")
load("../SCE_Robjects/10X8x10K.hsap.full.SCE.Robj")
X10x8x10K.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

X10x8x10K.hsap.count <- as.data.frame(counts(X10x8x10K.hsap))
X10x8x10K.hsap.count.mapped <- mapIDs(X10x8x10K.hsap.count, "hsap")
X10x8x10K.hsap.metadata <- as.data.frame(colData(X10x8x10K.hsap))
X10x8x10K.hsap.metadata[is.na(X10x8x10K.hsap.metadata)] <- 0
X10x8x10K.hsap.tmappedReads <- rowSums(X10x8x10K.hsap.metadata[, c(2:4,6)])
X10x8x10K.hsap.mappedReadsPerc <- (X10x8x10K.hsap.tmappedReads/X10x8x10K.hsap.metadata$nTReads)*100
X10x8x10K.hsap.metadata <- add_column(X10x8x10K.hsap.metadata, mappedReads= X10x8x10K.hsap.tmappedReads, mappedReadsPerc= X10x8x10K.hsap.mappedReadsPerc, .after = 1 )
X10x8x10K.hsap.cells <- which(X10x8x10K.hsap.metadata$Species == "Human")
X10x8x10K.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(X10x8x10K.hsap.count.mapped[,X10x8x10K.hsap.cells])), colData = X10x8x10K.hsap.metadata[X10x8x10K.hsap.cells,c(1:15)])

X10x8x10K.hsap.is.mito <- grep("^MT-", rownames(counts(X10x8x10K.hsap.new)))
X10x8x10K.hsap.new <- calculateQCMetrics(X10x8x10K.hsap.new, feature_controls=list(Mt=X10x8x10K.hsap.is.mito))
X10x8x10K.hsap.nTReads.drop <- which(isOutlier(X10x8x10K.hsap.new$mappedReads, nmads=2, type="lower", log=TRUE))
X10x8x10K.hsap.mappedPct.drop <- which(X10x8x10K.hsap.new$mappedReadsPerc < 65)
X10x8x10K.hsap.mito.drop <- which(isOutlier(X10x8x10K.hsap.new$log10_total_counts_Mt, nmads=2, type="higher"))
#plot(density(colData(X10x8x10K.hsap.new)[-X10x8x10K.hsap.mito.drop, "log10_total_counts_Mt"]))
X10x8x10K.hsap.final <- X10x8x10K.hsap.new[,-c(X10x8x10K.hsap.nTReads.drop, X10x8x10K.hsap.mappedPct.drop,X10x8x10K.hsap.mito.drop)]
X10x8x10K.hsap.fin.metadata <- colData(X10x8x10K.hsap.final)

#### end ####


# X10Scilife hsap Preparation ====
print("X10Scilife is Running...")
load("../SCE_Robjects/10XScilife.hsap.full.SCE.Robj")
X10Scilife.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

X10Scilife.hsap.count <- as.data.frame(counts(X10Scilife.hsap))
X10Scilife.hsap.count.mapped <- mapIDs(X10Scilife.hsap.count, "hsap")
X10Scilife.hsap.metadata <- as.data.frame(colData(X10Scilife.hsap))
X10Scilife.hsap.metadata[is.na(X10Scilife.hsap.metadata)] <- 0
X10Scilife.hsap.tmappedReads <- rowSums(X10Scilife.hsap.metadata[, c(2:4,6)])
X10Scilife.hsap.mappedReadsPerc <- (X10Scilife.hsap.tmappedReads/X10Scilife.hsap.metadata$nTReads)*100
X10Scilife.hsap.metadata <- add_column(X10Scilife.hsap.metadata, mappedReads= X10Scilife.hsap.tmappedReads, mappedReadsPerc= X10Scilife.hsap.mappedReadsPerc, .after = 1 )
X10Scilife.hsap.cells <- which(X10Scilife.hsap.metadata$Species == "Human")
X10Scilife.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(X10Scilife.hsap.count.mapped[,X10Scilife.hsap.cells])), colData = X10Scilife.hsap.metadata[X10Scilife.hsap.cells,c(1:15)])

X10Scilife.hsap.is.mito <- grep("^MT-", rownames(counts(X10Scilife.hsap.new)))
X10Scilife.hsap.new <- calculateQCMetrics(X10Scilife.hsap.new, feature_controls=list(Mt=X10Scilife.hsap.is.mito))
X10Scilife.hsap.nTReads.drop <- which(isOutlier(X10Scilife.hsap.new$mappedReads, nmads=2, type="lower", log=TRUE))
X10Scilife.hsap.mappedPct.drop <- which(X10Scilife.hsap.new$mappedReadsPerc < 65)
X10Scilife.hsap.mito.drop <- which(isOutlier(X10Scilife.hsap.new$log10_total_counts_Mt, nmads=3, type="higher"))
#plot(density(colData(X10Scilife.hsap.new)[-X10Scilife.hsap.mito.drop, "log10_total_counts_Mt"]))
X10Scilife.hsap.final <- X10Scilife.hsap.new[,-c(X10Scilife.hsap.nTReads.drop, X10Scilife.hsap.mappedPct.drop,X10Scilife.hsap.mito.drop)]
X10Scilife.hsap.fin.metadata <- colData(X10Scilife.hsap.final)

#### end ####


#library(data.table)
# Read distribution in each sample in each tech ====
#print("Read distribution in each sample in each tech")
#hsap.read.percentages <- data.frame(MARSseq.hsap.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("MARSseq", nrow(MARSseq.hsap.fin.metadata)))
#hsap.read.percentages <- rbind(hsap.read.percentages, data.frame(CELseq2.hsap.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("CELseq2", nrow(CELseq2.hsap.fin.metadata))))
#hsap.read.percentages <- rbind(hsap.read.percentages, data.frame(QUARTZseq.hsap.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("QUARTZseq", nrow(QUARTZseq.hsap.fin.metadata))))
#hsap.read.percentages <- rbind(hsap.read.percentages, data.frame(Dropseq.hsap.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("Dropseq", nrow(Dropseq.hsap.fin.metadata))))
#hsap.read.percentages <- rbind(hsap.read.percentages, data.frame(SCRBseq.hsap.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("SCRBseq", nrow(SCRBseq.hsap.fin.metadata))))
#hsap.read.percentages <- rbind(hsap.read.percentages, data.frame(SeqwellV1.hsap.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("SeqwellV2", nrow(SeqwellV2.hsap.fin.metadata))))
#hsap.read.percentages <- rbind(hsap.read.percentages, data.frame(SeqwellV2.hsap.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("SeqwellV2", nrow(SeqwellV2.hsap.fin.metadata))))
#hsap.read.percentages <- rbind(hsap.read.percentages, data.frame(Nuclei10X.hsap.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("Nuclei10X", nrow(Nuclei10X.hsap.fin.metadata))))
#hsap.read.percentages <- rbind(hsap.read.percentages, data.frame(ICELL8.hsap.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("ICELL8", nrow(ICELL8.hsap.fin.metadata))))
#hsap.read.percentages <- rbind(hsap.read.percentages, data.frame(ddSEQ.hsap.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("ddSEQ", nrow(ddSEQ.hsap.fin.metadata))))
#hsap.read.percentages <- rbind(hsap.read.percentages, data.frame(ddSEQexp1.hsap.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.hsap.fin.metadata))))
#hsap.read.percentages <- rbind(hsap.read.percentages, data.frame(C1HT.hsap.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("C1HT", nrow(C1HT.hsap.fin.metadata))))
#hsap.read.percentages <- rbind(hsap.read.percentages, data.frame(X10x8x10K.hsap.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("X10x8x10K", nrow(X10x8x10K.hsap.fin.metadata))))
#hsap.read.percentages <- rbind(hsap.read.percentages, data.frame(X10Scilife.hsap.fin.metadata[, c("nTReads", "nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads","Library")], tech = rep("X10Scilife", nrow(X10Scilife.hsap.fin.metadata))))
#hsap.melted.read.percentages <- melt(hsap.read.percentages, id = c("Library","tech"), measure = c("nExonReads", "nIntronReads", "nIntergenicReads", "nUnmappedReads", "nAmbiguityReads", "nMultimapReads"))

#pdf("all_techs/hsap_read_distributions_all.pdf")
#ggplot(data=hsap.melted.read.percentages, aes(Library, value, fill= variable)) + geom_bar(stat = "identity", width = 0.5, position = "stack") + facet_grid(~ tech, scales = "free")+theme (axis.text.x = element_text(angle = 90, hjust = 1))
#ggplot(data=hsap.melted.read.percentages, aes(tech, value, fill= variable)) + geom_bar(stat = "identity", width = 0.5, position = "stack") + facet_grid(~ tech, scales = "free")+theme (axis.text.x = element_text(angle = 90, hjust = 1))
#dev.off()





# Total number of reads per tech per lib ====
print("Total number of reads per tech per lib")
hsap.total.reads <- data.frame(MARSseq.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("MARSseq", nrow(MARSseq.hsap.fin.metadata)))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(CELseq2.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("CELseq2", nrow(CELseq2.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(QUARTZseq.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("QUARTZseq", nrow(QUARTZseq.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(Dropseq.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("Dropseq", nrow(Dropseq.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(SCRBseq.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("SCRBseq", nrow(SCRBseq.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(SeqwellV1.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("SeqwellV1", nrow(SeqwellV1.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(SeqwellV2.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("SeqwellV2", nrow(SeqwellV2.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(Nuclei10X.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("Nuclei10X", nrow(Nuclei10X.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(ICELL8.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("ICELL8", nrow(ICELL8.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(ddSEQ.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("ddSEQ", nrow(ddSEQ.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(ddSEQexp1.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(C1HT.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("C1HT", nrow(C1HT.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(X10x8x10K.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("X10x8x10K", nrow(X10x8x10K.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(X10Scilife.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("X10Scilife", nrow(X10Scilife.hsap.fin.metadata))))
pdf("all_techs/hsap_total_number_Reads_all.pdf")
#ggplot(data=hsap.total.reads, aes(x=Library, y=log(nTReads), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
ggplot(data=hsap.total.reads, aes(x=tech, y=log(nTReads), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
dev.off()




#Total number of genes per tech per lib ====
print("Total number of genes per tech per lib ")
hsap.total.genes <- data.frame(MARSseq.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("MARSseq", nrow(MARSseq.hsap.fin.metadata)))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(CELseq2.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("CELseq2", nrow(CELseq2.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(QUARTZseq.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("QUARTZseq", nrow(QUARTZseq.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(Dropseq.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("Dropseq", nrow(Dropseq.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(SCRBseq.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("SCRBseq", nrow(SCRBseq.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(SeqwellV1.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("SeqwellV1", nrow(SeqwellV1.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(SeqwellV2.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("SeqwellV2", nrow(SeqwellV2.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(Nuclei10X.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("Nuclei10X", nrow(Nuclei10X.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(ICELL8.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("ICELL8", nrow(ICELL8.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(ddSEQ.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("ddSEQ", nrow(ddSEQ.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(ddSEQexp1.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(C1HT.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("C1HT", nrow(C1HT.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(X10x8x10K.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("X10x8x10K", nrow(X10x8x10K.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(X10Scilife.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("X10Scilife", nrow(X10Scilife.hsap.fin.metadata))))
pdf("all_techs/hsap_total_number_Genes_all.pdf")
#ggplot(data=hsap.total.genes, aes(x=Library, y=log(nGenes), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
ggplot(data=hsap.total.genes, aes(x=tech, y=log(nGenes), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
dev.off()


# Total number of UMIs per tech per lib ====
print("Total number of UMIs per tech per lib")
hsap.total.UMIs <- data.frame(MARSseq.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("MARSseq", nrow(MARSseq.hsap.fin.metadata)))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(CELseq2.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("CELseq2", nrow(CELseq2.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(QUARTZseq.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("QUARTZseq", nrow(QUARTZseq.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(Dropseq.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("Dropseq", nrow(Dropseq.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(SCRBseq.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("SCRBseq", nrow(SCRBseq.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(SeqwellV1.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("SeqwellV1", nrow(SeqwellV1.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(SeqwellV2.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("SeqwellV2", nrow(SeqwellV2.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(Nuclei10X.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("Nuclei10X", nrow(Nuclei10X.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(ICELL8.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("ICELL8", nrow(ICELL8.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(ddSEQ.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("ddSEQ", nrow(ddSEQ.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(ddSEQexp1.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(C1HT.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("C1HT", nrow(C1HT.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(X10x8x10K.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("X10x8x10K", nrow(X10x8x10K.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(X10Scilife.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("X10Scilife", nrow(X10Scilife.hsap.fin.metadata))))
pdf("all_techs/hsap_total_number_UMIs_all.pdf")
#ggplot(data=hsap.total.UMIs, aes(x=Library, y=log(nUMIs), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free")
ggplot(data=hsap.total.UMIs, aes(x=tech, y=log(nUMIs), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free")
dev.off()


# nTreads vs nUMI curves ====
print("nTreads vs nUMI curves")
hsap.nread.numi.df <- data.frame(MARSseq.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("MARSseq", nrow(MARSseq.hsap.fin.metadata)))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(CELseq2.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("CELseq2", nrow(CELseq2.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(QUARTZseq.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("QUARTZseq", nrow(QUARTZseq.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(Dropseq.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("Dropseq", nrow(Dropseq.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(SCRBseq.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("SCRBseq", nrow(SCRBseq.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(SeqwellV1.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("SeqwellV1", nrow(SeqwellV1.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(SeqwellV2.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("SeqwellV2", nrow(SeqwellV2.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(Nuclei10X.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("Nuclei10X", nrow(Nuclei10X.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(ICELL8.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("ICELL8", nrow(ICELL8.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(ddSEQ.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("ddSEQ", nrow(ddSEQ.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(ddSEQexp1.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(C1HT.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("C1HT", nrow(C1HT.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(X10x8x10K.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("X10x8x10K", nrow(X10x8x10K.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(X10Scilife.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("X10Scilife", nrow(X10Scilife.hsap.fin.metadata))))

pdf("all_techs/hsap_nReadsVSnUMIs_all.pdf")
ggplot(hsap.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(aes(color=tech))+ xlim(0, 2000000)# + geom_line(aes(color=tech))#
dev.off()



# nTreads vs nGenes curves ====
print("nTreads vs nGenes curves")
hsap.nread.ngene.df <- data.frame(MARSseq.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("MARSseq", nrow(MARSseq.hsap.fin.metadata)))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(CELseq2.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("CELseq2", nrow(CELseq2.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(QUARTZseq.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("QUARTZseq", nrow(QUARTZseq.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(Dropseq.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("Dropseq", nrow(Dropseq.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(SCRBseq.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("SCRBseq", nrow(SCRBseq.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(SeqwellV1.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("SeqwellV1", nrow(SeqwellV1.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(SeqwellV2.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("SeqwellV2", nrow(SeqwellV2.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(Nuclei10X.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("Nuclei10X", nrow(Nuclei10X.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(ICELL8.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("ICELL8", nrow(ICELL8.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(ddSEQ.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("ddSEQ", nrow(ddSEQ.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(ddSEQexp1.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(C1HT.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("C1HT", nrow(C1HT.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(X10x8x10K.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("X10x8x10K", nrow(X10x8x10K.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(X10Scilife.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("X10Scilife", nrow(X10Scilife.hsap.fin.metadata))))
pdf("all_techs/hsap_nReadsVSnGenes_all.pdf")
#ggplot(hsap.nread.ngene.df, aes(x=nGenes, y=nTReads, group=tech)) + geom_point(aes(color=tech)) + geom_smooth(aes(color=tech))+ xlim(0, 18000)# + geom_line(aes(color=tech))#
ggplot(hsap.nread.ngene.df, aes(x=nTReads, y=nGenes, group=tech)) + geom_smooth(aes(color=tech))+ xlim(0, 2000000)# + geom_line(aes(color=tech))#
dev.off()