library(scater)
library(ggplot2)
library(biomaRt)
library(tibble)
library(data.table)
require(minpack.lm)

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


# loading Annotated seurat objects ####
load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/MARSseq_data_seu.obj_res0.5_dim6.RData")
MARSseq.obj <- data
rm(data)
MARSseq.metadata <- MARSseq.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/CELseq2_data_seu.obj_res0.5_dim8.RData")
CELseq2.obj <- data
rm(data)
CELseq2.metadata <- CELseq2.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/QUARTZseq_data_seu.obj_res0.4_dim10.RData")
QUARTZseq.obj <- data
rm(data)
QUARTZseq.metadata <- QUARTZseq.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/Dropseq_data_seu.obj_res0.4_dim8.RData")
Dropseq.obj <- data
rm(data)
Dropseq.metadata <- Dropseq.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/SCRBseq_data_seu.obj_res0.3_dim10.RData")
SCRBseq.obj <- data
rm(data)
SCRBseq.metadata <- SCRBseq.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/SeqwellV2_data_seu.obj_res0.5_dim7.RData")
SeqwellV2.obj <- data
rm(data)
SeqwellV2.metadata <- SeqwellV2.obj@meta.data


load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/X10Scilife_data_seu.obj_res0.1_dim8.RData")
X10Scilife.obj <- data
rm(data)
X10Scilife.metadata <- X10Scilife.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/X10Nuclei_data_seu.obj_res0.4_dim8.RData")
X10Nuclei.obj <- data
rm(data)
X10Nuclei.metadata <- X10Nuclei.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/ICELL8_data_seu.obj_res0.5_dim9.RData")
ICELL8.obj <- data
rm(data)
ICELL8.metadata <- ICELL8.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/ddSEQ_data_seu.obj_res0.5_dim13.RData")
ddSEQ.obj <- data
rm(data)
ddSEQ.metadata <- ddSEQ.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/ddSEQexp1_data_seu.obj_res0.3_dim6.RData")
ddSEQexp1.obj <- data
rm(data)
ddSEQexp1.metadata <- ddSEQexp1.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/C1HTsmall_data_seu.obj_res0.3_dim6.RData")
C1HTsmall.obj <- data
rm(data)
C1HTsmall.metadata <- C1HTsmall.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/C1HTmedium_data_seu.obj_res0.3_dim8.RData")
C1HTmedium.obj <- data
rm(data)
C1HTmedium.metadata <- C1HTmedium.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/CB1_data_seu.obj_res0.5_dim7.RData")
CB1.obj <- data
rm(data)
CB1.metadata <- CB1.obj@meta.data



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
MARSseq.hsap.fin.metadata <- MARSseq.hsap.fin.metadata[rownames(MARSseq.metadata),]
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
CELseq2.hsap.fin.metadata <- CELseq2.hsap.fin.metadata[rownames(CELseq2.metadata),]

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
QUARTZseq.hsap.fin.metadata <- QUARTZseq.hsap.fin.metadata[rownames(QUARTZseq.metadata),]

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
Dropseq.hsap.fin.metadata <- Dropseq.hsap.fin.metadata[rownames(Dropseq.metadata),]

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
SCRBseq.hsap.fin.metadata <- SCRBseq.hsap.fin.metadata[rownames(SCRBseq.metadata),]

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
SeqwellV2.hsap.fin.metadata <- SeqwellV2.hsap.fin.metadata[rownames(SeqwellV2.metadata),]

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
#SeqwellV1.hsap.fin.metadata <- SeqwellV1.hsap.fin.metadata[rownames(SeqwellV1.metadata),]

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
Nuclei10X.hsap.fin.metadata <- Nuclei10X.hsap.fin.metadata[rownames(X10Nuclei.metadata),]

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
ICELL8.hsap.fin.metadata <- ICELL8.hsap.fin.metadata[rownames(ICELL8.metadata),]

                                                                                                       
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
ddSEQ.hsap.fin.metadata <- ddSEQ.hsap.fin.metadata[rownames(ddSEQ.metadata),]

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
ddSEQexp1.hsap.fin.metadata <- ddSEQexp1.hsap.fin.metadata[rownames(ddSEQexp1.metadata),]

#### end ####


# C1HTsmall hsap Preparation ====
print("C1HTsmall is Running...")
load("../SCE_Robjects/C1HTsmall.hsap.full.SCE.Robj")
C1HTsmall.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

C1HTsmall.hsap.count <- as.data.frame(counts(C1HTsmall.hsap))
C1HTsmall.hsap.count.mapped <- mapIDs(C1HTsmall.hsap.count, "hsap")
C1HTsmall.hsap.metadata <- as.data.frame(colData(C1HTsmall.hsap))
C1HTsmall.hsap.metadata[is.na(C1HTsmall.hsap.metadata)] <- 0
C1HTsmall.hsap.tmappedReads <- rowSums(C1HTsmall.hsap.metadata[, c(2:4,6)])
C1HTsmall.hsap.mappedReadsPerc <- (C1HTsmall.hsap.tmappedReads/C1HTsmall.hsap.metadata$nTReads)*100
C1HTsmall.hsap.metadata <- add_column(C1HTsmall.hsap.metadata, mappedReads= C1HTsmall.hsap.tmappedReads, mappedReadsPerc= C1HTsmall.hsap.mappedReadsPerc, .after = 1 )
C1HTsmall.hsap.cells <- which(C1HTsmall.hsap.metadata$Species == "Human")
C1HTsmall.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(C1HTsmall.hsap.count.mapped[,C1HTsmall.hsap.cells])), colData = C1HTsmall.hsap.metadata[C1HTsmall.hsap.cells,c(1:15)])

C1HTsmall.hsap.is.mito <- grep("^MT-", rownames(counts(C1HTsmall.hsap.new)))
C1HTsmall.hsap.new <- calculateQCMetrics(C1HTsmall.hsap.new, feature_controls=list(Mt=C1HTsmall.hsap.is.mito))
C1HTsmall.hsap.nTReads.drop <- which(isOutlier(C1HTsmall.hsap.new$mappedReads, nmads=2, type="lower", log=TRUE))
C1HTsmall.hsap.mappedPct.drop <- which(C1HTsmall.hsap.new$mappedReadsPerc < 65)
C1HTsmall.hsap.mito.drop <- which(isOutlier(C1HTsmall.hsap.new$log10_total_counts_Mt, nmads=2, type="higher"))
#plot(density(colData(C1HTsmall.hsap.new)[-C1HTsmall.hsap.mito.drop, "log10_total_counts_Mt"]))
C1HTsmall.hsap.final <- C1HTsmall.hsap.new[,-c(C1HTsmall.hsap.nTReads.drop, C1HTsmall.hsap.mappedPct.drop,C1HTsmall.hsap.mito.drop)]
C1HTsmall.hsap.fin.metadata <- colData(C1HTsmall.hsap.final)
C1HTsmall.hsap.fin.metadata <- C1HTsmall.hsap.fin.metadata[rownames(C1HTsmall.metadata),]

#### end ####

# C1HTmedium hsap Preparation ====
print("C1HTmedium is Running...")
load("../SCE_Robjects/C1HTmedium.hsap.full.SCE.Robj")
C1HTmedium.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

C1HTmedium.hsap.count <- as.data.frame(counts(C1HTmedium.hsap))
C1HTmedium.hsap.count.mapped <- mapIDs(C1HTmedium.hsap.count, "hsap")
C1HTmedium.hsap.metadata <- as.data.frame(colData(C1HTmedium.hsap))
C1HTmedium.hsap.metadata[is.na(C1HTmedium.hsap.metadata)] <- 0
C1HTmedium.hsap.tmappedReads <- rowSums(C1HTmedium.hsap.metadata[, c(2:4,6)])
C1HTmedium.hsap.mappedReadsPerc <- (C1HTmedium.hsap.tmappedReads/C1HTmedium.hsap.metadata$nTReads)*100
C1HTmedium.hsap.metadata <- add_column(C1HTmedium.hsap.metadata, mappedReads= C1HTmedium.hsap.tmappedReads, mappedReadsPerc= C1HTmedium.hsap.mappedReadsPerc, .after = 1 )
C1HTmedium.hsap.cells <- which(C1HTmedium.hsap.metadata$Species == "Human")
C1HTmedium.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(C1HTmedium.hsap.count.mapped[,C1HTmedium.hsap.cells])), colData = C1HTmedium.hsap.metadata[C1HTmedium.hsap.cells,c(1:15)])

C1HTmedium.hsap.is.mito <- grep("^MT-", rownames(counts(C1HTmedium.hsap.new)))
C1HTmedium.hsap.new <- calculateQCMetrics(C1HTmedium.hsap.new, feature_controls=list(Mt=C1HTmedium.hsap.is.mito))
C1HTmedium.hsap.nTReads.drop <- which(isOutlier(C1HTmedium.hsap.new$mappedReads, nmads=2, type="lower", log=TRUE))
C1HTmedium.hsap.mappedPct.drop <- which(C1HTmedium.hsap.new$mappedReadsPerc < 65)
C1HTmedium.hsap.mito.drop <- which(isOutlier(C1HTmedium.hsap.new$log10_total_counts_Mt, nmads=2, type="higher"))
#plot(density(colData(C1HTmedium.hsap.new)[-C1HTmedium.hsap.mito.drop, "log10_total_counts_Mt"]))
C1HTmedium.hsap.final <- C1HTmedium.hsap.new[,-c(C1HTmedium.hsap.nTReads.drop, C1HTmedium.hsap.mappedPct.drop,C1HTmedium.hsap.mito.drop)]
C1HTmedium.hsap.fin.metadata <- colData(C1HTmedium.hsap.final)
C1HTmedium.hsap.fin.metadata <- C1HTmedium.hsap.fin.metadata[rownames(C1HTmedium.metadata),]

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
X10x8x10K.hsap.fin.metadata <- X10x8x10K.hsap.fin.metadata[rownames(X10x8x10K.metadata),]

#### end ####

#10X2x5Half hsap Preparation ====
print("10X2x5Half is Running...")
load("../SCE_Robjects/10X2x5Half.hsap.full.SCE.Robj")
X10X2x5Half.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

X10X2x5Half.hsap.count <- as.data.frame(counts(X10X2x5Half.hsap))
X10X2x5Half.hsap.count.mapped <- mapIDs(X10X2x5Half.hsap.count, "hsap")
X10X2x5Half.hsap.metadata <- as.data.frame(colData(X10X2x5Half.hsap))
X10X2x5Half.hsap.metadata[is.na(X10X2x5Half.hsap.metadata)] <- 0
X10X2x5Half.hsap.tmappedReads <- rowSums(X10X2x5Half.hsap.metadata[, c(2:4,6)])
X10X2x5Half.hsap.mappedReadsPerc <- (X10X2x5Half.hsap.tmappedReads/X10X2x5Half.hsap.metadata$nTReads)*100
X10X2x5Half.hsap.metadata <- add_column(X10X2x5Half.hsap.metadata, mappedReads= X10X2x5Half.hsap.tmappedReads, mappedReadsPerc= X10X2x5Half.hsap.mappedReadsPerc, .after = 1 )
X10X2x5Half.hsap.cells <- which(X10X2x5Half.hsap.metadata$Species == "Human")
X10X2x5Half.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(X10X2x5Half.hsap.count.mapped[,X10X2x5Half.hsap.cells])), colData = X10X2x5Half.hsap.metadata[X10X2x5Half.hsap.cells,c(1:15)])

X10X2x5Half.hsap.is.mito <- grep("^MT-", rownames(counts(X10X2x5Half.hsap.new)))
X10X2x5Half.hsap.new <- calculateQCMetrics(X10X2x5Half.hsap.new, feature_controls=list(Mt=X10X2x5Half.hsap.is.mito))
X10X2x5Half.hsap.nTReads.drop <- which(isOutlier(X10X2x5Half.hsap.new$mappedReads, nmads=2, type="lower", log=TRUE))
X10X2x5Half.hsap.mappedPct.drop <- which(X10X2x5Half.hsap.new$mappedReadsPerc < 65)
X10X2x5Half.hsap.mito.drop <- which(isOutlier(X10X2x5Half.hsap.new$log10_total_counts_Mt, nmads=2, type="higher"))
#plot(density(colData(X10X2x5Half.hsap.new)[-X10X2x5Half.hsap.mito.drop, "log10_total_counts_Mt"]))
X10X2x5Half.hsap.final <- X10X2x5Half.hsap.new[,-c(X10X2x5Half.hsap.nTReads.drop, X10X2x5Half.hsap.mappedPct.drop,X10X2x5Half.hsap.mito.drop)]
X10X2x5Half.hsap.fin.metadata <- colData(X10X2x5Half.hsap.final)
#X10X2x5Half.hsap.fin.metadata <- X10X2x5Half.hsap.fin.metadata[rownames(X10X2x5Half.metadata),]

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
X10Scilife.hsap.fin.metadata <- X10Scilife.hsap.fin.metadata[rownames(X10Scilife.metadata),]

#### end ####


# CB1 hsap Preparation ====
print("CB1 is Running...")
load("../SCE_Robjects/1CB.hsap.full.SCE.Robj")
CB1.hsap <- full.SCE.hsap
rm(full.SCE.hsap)

CB1.hsap.count <- as.data.frame(counts(CB1.hsap))
CB1.hsap.count.mapped <- mapIDs(CB1.hsap.count, "hsap")
CB1.hsap.metadata <- as.data.frame(colData(CB1.hsap))
CB1.hsap.metadata[is.na(CB1.hsap.metadata)] <- 0
CB1.hsap.tmappedReads <- rowSums(CB1.hsap.metadata[, c(2:4,6)])
CB1.hsap.mappedReadsPerc <- (CB1.hsap.tmappedReads/CB1.hsap.metadata$nTReads)*100
CB1.hsap.metadata <- add_column(CB1.hsap.metadata, mappedReads= CB1.hsap.tmappedReads, mappedReadsPerc= CB1.hsap.mappedReadsPerc, .after = 1 )
CB1.hsap.cells <- which(CB1.hsap.metadata$Species == "Human")
CB1.hsap.new <- SingleCellExperiment(assays = list(counts = as.matrix(CB1.hsap.count.mapped[,CB1.hsap.cells])), colData = CB1.hsap.metadata[CB1.hsap.cells,c(1:15)])

CB1.hsap.is.mito <- grep("^MT-", rownames(counts(CB1.hsap.new)))
CB1.hsap.new <- calculateQCMetrics(CB1.hsap.new, feature_controls=list(Mt=CB1.hsap.is.mito))
CB1.hsap.nTReads.drop <- which(isOutlier(CB1.hsap.new$mappedReads, nmads=2, type="lower", log=TRUE))
CB1.hsap.mappedPct.drop <- which(CB1.hsap.new$mappedReadsPerc < 65)
CB1.hsap.mito.drop <- which(isOutlier(CB1.hsap.new$log10_total_counts_Mt, nmads=2, type="higher"))
#plot(density(colData(CB1.hsap.new)[-CB1.hsap.mito.drop, "log10_total_counts_Mt"]))
CB1.hsap.final <- CB1.hsap.new[,-c(CB1.hsap.nTReads.drop, CB1.hsap.mappedPct.drop,CB1.hsap.mito.drop)]
CB1.hsap.fin.metadata <- colData(CB1.hsap.final)
CB1.hsap.fin.metadata <- CB1.hsap.fin.metadata[rownames(CB1.metadata),]

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

#pdf("all_techs_V2/hsap_read_distributions_all.pdf")
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
#hsap.total.reads <- rbind(hsap.total.reads, data.frame(SeqwellV1.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("SeqwellV1", nrow(SeqwellV1.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(SeqwellV2.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("SeqwellV2", nrow(SeqwellV2.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(Nuclei10X.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("Nuclei10X", nrow(Nuclei10X.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(ICELL8.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("ICELL8", nrow(ICELL8.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(ddSEQ.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("ddSEQ", nrow(ddSEQ.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(ddSEQexp1.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(C1HTsmall.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("C1HTsmall", nrow(C1HTsmall.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(C1HTmedium.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("C1HTmedium", nrow(C1HTmedium.hsap.fin.metadata))))
#hsap.total.reads <- rbind(hsap.total.reads, data.frame(X10x8x10K.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("X10x8x10K", nrow(X10x8x10K.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(X10X2x5Half.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("X10X2x5Half", nrow(X10X2x5Half.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(X10Scilife.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("X10Scilife", nrow(X10Scilife.hsap.fin.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(CB1.hsap.fin.metadata[, c("Library", "nTReads")], tech = rep("CB1", nrow(CB1.hsap.fin.metadata))))
pdf("all_techs_V2/hsap_total_number_Reads_all_V7.pdf")
#ggplot(data=hsap.total.reads, aes(x=Library, y=log(nTReads), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
#ggplot(data=hsap.total.reads, aes(x=tech, y=log(nTReads), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
ggplot(data=hsap.total.reads, aes(x=tech, y=nTReads, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(1000, 20000, 100000, 500000, 1000000)) + facet_grid(. ~ tech, scales = "free") + geom_hline(yintercept=c(1000, 20000, 100000, 500000, 1000000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
#ggplot(data=hsap.total.reads, aes(x=tech, y=nTReads, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(labels= scales::comma, breaks = c(1000, 20000, 100000, 500000, 1000000)) + facet_grid(. ~ tech, scales = "free") + geom_hline(yintercept=c(1000, 20000, 100000, 500000, 1000000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
dev.off()




#Total number of genes per tech per lib ====
print("Total number of genes per tech per lib ")
hsap.total.genes <- data.frame(MARSseq.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("MARSseq", nrow(MARSseq.hsap.fin.metadata)))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(CELseq2.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("CELseq2", nrow(CELseq2.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(QUARTZseq.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("QUARTZseq", nrow(QUARTZseq.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(Dropseq.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("Dropseq", nrow(Dropseq.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(SCRBseq.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("SCRBseq", nrow(SCRBseq.hsap.fin.metadata))))
#hsap.total.genes <- rbind(hsap.total.genes, data.frame(SeqwellV1.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("SeqwellV1", nrow(SeqwellV1.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(SeqwellV2.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("SeqwellV2", nrow(SeqwellV2.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(Nuclei10X.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("Nuclei10X", nrow(Nuclei10X.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(ICELL8.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("ICELL8", nrow(ICELL8.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(ddSEQ.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("ddSEQ", nrow(ddSEQ.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(ddSEQexp1.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(C1HTsmall.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("C1HTsmall", nrow(C1HTsmall.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(C1HTmedium.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("C1HTmedium", nrow(C1HTmedium.hsap.fin.metadata))))
#hsap.total.genes <- rbind(hsap.total.genes, data.frame(X10x8x10K.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("X10x8x10K", nrow(X10x8x10K.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(X10X2x5Half.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("X10X2x5Half", nrow(X10X2x5Half.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(X10Scilife.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("X10Scilife", nrow(X10Scilife.hsap.fin.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(CB1.hsap.fin.metadata[, c("Library", "nGenes")], tech = rep("CB1", nrow(CB1.hsap.fin.metadata))))
pdf("all_techs_V2/hsap_total_number_Genes_all_V7.pdf")
#ggplot(data=hsap.total.genes, aes(x=Library, y=log(nGenes), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
#ggplot(data=hsap.total.genes, aes(x=tech, y=log(nGenes), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
ggplot(data=hsap.total.genes, aes(x=tech, y=nGenes, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(10,150,500, 1500, 5000)) + facet_grid(. ~ tech, scales = "free") + geom_hline(yintercept=c(10,150,500, 1500, 5000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
#ggplot(data=hsap.total.genes, aes(x=tech, y=nGenes, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(labels= scales::comma, breaks = c(10,150,500, 1500, 5000)) + facet_grid(. ~ tech, scales = "free") + geom_hline(yintercept=c(10,150,500, 1500, 5000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
dev.off()


# Total number of UMIs per tech per lib ====
print("Total number of UMIs per tech per lib")
hsap.total.UMIs <- data.frame(MARSseq.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("MARSseq", nrow(MARSseq.hsap.fin.metadata)))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(CELseq2.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("CELseq2", nrow(CELseq2.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(QUARTZseq.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("QUARTZseq", nrow(QUARTZseq.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(Dropseq.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("Dropseq", nrow(Dropseq.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(SCRBseq.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("SCRBseq", nrow(SCRBseq.hsap.fin.metadata))))
#hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(SeqwellV1.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("SeqwellV1", nrow(SeqwellV1.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(SeqwellV2.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("SeqwellV2", nrow(SeqwellV2.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(Nuclei10X.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("Nuclei10X", nrow(Nuclei10X.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(ICELL8.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("ICELL8", nrow(ICELL8.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(ddSEQ.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("ddSEQ", nrow(ddSEQ.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(ddSEQexp1.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(C1HTsmall.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("C1HTsmall", nrow(C1HTsmall.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(C1HTmedium.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("C1HTmedium", nrow(C1HTmedium.hsap.fin.metadata))))
#hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(X10x8x10K.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("X10x8x10K", nrow(X10x8x10K.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(X10X2x5Half.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("X10X2x5Half", nrow(X10X2x5Half.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(X10Scilife.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("X10Scilife", nrow(X10Scilife.hsap.fin.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(CB1.hsap.fin.metadata[, c("Library", "nUMIs")], tech = rep("CB1", nrow(CB1.hsap.fin.metadata))))
pdf("all_techs_V2/hsap_total_number_UMIs_all_V7.pdf")
#ggplot(data=hsap.total.UMIs, aes(x=Library, y=log(nUMIs), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free")
#ggplot(data=hsap.total.UMIs, aes(x=tech, y=log(nUMIs), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free")
ggplot(data=hsap.total.UMIs, aes(x=tech, y=nUMIs, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(200,2000,20000,200000, 1000000)) + facet_grid(. ~ tech, scales = "free")+ geom_hline(yintercept=c(200,2000,20000,200000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
#ggplot(data=hsap.total.UMIs, aes(x=tech, y=nUMIs, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(labels= scales::comma, breaks = c(200,2000,20000,200000, 1000000)) + facet_grid(. ~ tech, scales = "free")+ geom_hline(yintercept=c(200,2000,20000,200000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
dev.off()


# nTreads vs nUMI curves ====
print("nTreads vs nUMI curves")
hsap.nread.numi.df <- data.frame(MARSseq.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("MARSseq", nrow(MARSseq.hsap.fin.metadata)))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(CELseq2.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("CELseq2", nrow(CELseq2.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(QUARTZseq.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("QUARTZseq", nrow(QUARTZseq.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(Dropseq.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("Dropseq", nrow(Dropseq.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(SCRBseq.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("SCRBseq", nrow(SCRBseq.hsap.fin.metadata))))
#hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(SeqwellV1.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("SeqwellV1", nrow(SeqwellV1.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(SeqwellV2.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("SeqwellV2", nrow(SeqwellV2.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(Nuclei10X.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("Nuclei10X", nrow(Nuclei10X.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(ICELL8.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("ICELL8", nrow(ICELL8.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(ddSEQ.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("ddSEQ", nrow(ddSEQ.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(ddSEQexp1.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(C1HTsmall.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("C1HTsmall", nrow(C1HTsmall.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(C1HTmedium.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("C1HTmedium", nrow(C1HTmedium.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(X10x8x10K.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("X10x8x10K", nrow(X10x8x10K.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(X10X2x5Half.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("X10X2x5Half", nrow(X10X2x5Half.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(X10Scilife.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("X10Scilife", nrow(X10Scilife.hsap.fin.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(CB1.hsap.fin.metadata[, c("nTReads", "nUMIs")], tech = rep("CB1", nrow(CB1.hsap.fin.metadata))))

pdf("all_techs_V2/hsap_nReadsVSnUMIs_all_V7.pdf")
#ggplot(hsap.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(aes(color=tech))+ xlim(0, 100000)# + geom_line(aes(color=tech))#
ggplot(hsap.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(method="lm", formula=y~x, fill="grey",aes(color=tech))+ xlim(0, 100000)# + geom_line(aes(color=tech))#
#ggplot(hsap.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(method="lm", formula=y~log(x), fill="grey", aes(color=tech))+ xlim(0, 2000000)# + geom_line(aes(color=tech))#
#ggplot(hsap.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(method='nlsLM',formula=y ~ a*x^b, se=FALSE, method.args = list(start = list(a=1,b=1)), aes(color=tech))+ xlim(0, 1500000)# + geom_line(aes(color=tech))#
#ggplot(hsap.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(method='nls',formula=y ~ SSasympOff(x, A, lrc, c0), se=FALSE, aes(color=tech))+ xlim(0, 1500000)# + geom_line(aes(color=tech))#

dev.off()



# nTreads vs nGenes curves ====
print("nTreads vs nGenes curves")
hsap.nread.ngene.df <- data.frame(MARSseq.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("MARSseq", nrow(MARSseq.hsap.fin.metadata)))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(CELseq2.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("CELseq2", nrow(CELseq2.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(QUARTZseq.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("QUARTZseq", nrow(QUARTZseq.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(Dropseq.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("Dropseq", nrow(Dropseq.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(SCRBseq.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("SCRBseq", nrow(SCRBseq.hsap.fin.metadata))))
#hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(SeqwellV1.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("SeqwellV1", nrow(SeqwellV1.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(SeqwellV2.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("SeqwellV2", nrow(SeqwellV2.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(Nuclei10X.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("Nuclei10X", nrow(Nuclei10X.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(ICELL8.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("ICELL8", nrow(ICELL8.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(ddSEQ.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("ddSEQ", nrow(ddSEQ.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(ddSEQexp1.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(C1HTsmall.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("C1HTsmall", nrow(C1HTsmall.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(C1HTmedium.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("C1HTmedium", nrow(C1HTmedium.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(X10x8x10K.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("X10x8x10K", nrow(X10x8x10K.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(X10X2x5Half.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("X10X2x5Half", nrow(X10X2x5Half.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(X10Scilife.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("X10Scilife", nrow(X10Scilife.hsap.fin.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(CB1.hsap.fin.metadata[, c("nTReads", "nGenes")], tech = rep("CB1", nrow(CB1.hsap.fin.metadata))))
pdf("all_techs_V2/hsap_nReadsVSnGenes_all_V7.pdf")
#ggplot(hsap.nread.ngene.df, aes(x=nGenes, y=nTReads, group=tech)) + geom_point(aes(color=tech)) + geom_smooth(aes(color=tech))+ xlim(0, 100000)# + geom_line(aes(color=tech))#
ggplot(hsap.nread.ngene.df, aes(x=nTReads, y=nGenes, group=tech)) + geom_smooth(method="lm", formula=y~x, fill="grey", aes(color=tech))+ xlim(0, 100000)# + geom_line(aes(color=tech))#
#ggplot(hsap.nread.ngene.df, aes(x=nTReads, y=nGenes, group=tech)) + geom_smooth(method="lm", formula=y~log(x), fill="grey", aes(color=tech))+ xlim(0, 2000000)# + geom_line(aes(color=tech))#
#ggplot(hsap.nread.ngene.df, aes(x=nTReads, y=nGenes, group=tech)) + geom_smooth(method='nlsLM',formula=y ~ a*x^b, se=FALSE, method.args = list(start = list(a=1,b=1)), aes(color=tech))+ xlim(0, 1500000)# + geom_line(aes(color=tech))#
#ggplot(hsap.nread.ngene.df, aes(x=nTReads, y=nGenes, group=tech)) + geom_smooth(method='nls',formula=y ~ SSasympOff(x, A, lrc, c0), se=FALSE, aes(color=tech))+ xlim(0, 1500000)# + geom_line(aes(color=tech))#
dev.off()