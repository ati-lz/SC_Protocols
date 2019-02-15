library(scater)
library(ggplot2)
library(ggExtra)
library(biomaRt)
library(tibble)
library(Seurat)
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

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/CB1_data_seu.obj_res0.5_dim7.RData")
CB1.obj <- data
rm(data)
CB1.metadata <- CB1.obj@meta.data

# End ####

# Taking out HEK cells of each technique ####
MARSseq.HEK <- names(MARSseq.obj@ident)[which(MARSseq.obj@ident == "HEK cells")]
CELseq2.HEK <- names(CELseq2.obj@ident)[which(CELseq2.obj@ident == "HEK cells")]
QUARTZseq.HEK <- names(QUARTZseq.obj@ident)[which(QUARTZseq.obj@ident == "HEK cells")]
Dropseq.HEK <- names(Dropseq.obj@ident)[which(Dropseq.obj@ident == "HEK cells")]
SCRBseq.HEK <- names(SCRBseq.obj@ident)[which(SCRBseq.obj@ident == "HEK cells")]
X10Scilife.HEK <- names(X10Scilife.obj@ident)[which(X10Scilife.obj@ident == "HEK cells")]
X10Nuclei.HEK <- names(X10Nuclei.obj@ident)[which(X10Nuclei.obj@ident == "HEK cells")]
ICELL8.HEK <- names(ICELL8.obj@ident)[which(ICELL8.obj@ident == "HEK cells")]
ddSEQ.HEK <- names(ddSEQ.obj@ident)[which(ddSEQ.obj@ident == "HEK cells")]
ddSEQexp1.HEK <- names(ddSEQexp1.obj@ident)[which(ddSEQexp1.obj@ident == "HEK cells")]
CB1.HEK <- names(CB1.obj@ident)[which(CB1.obj@ident == "HEK cells")]

# End ####

# Taking out Monocytes cells of each technique ####
MARSseq.Monocytes <- names(MARSseq.obj@ident)[which(MARSseq.obj@ident == "CD14+ Monocytes")]
CELseq2.Monocytes <- names(CELseq2.obj@ident)[which(CELseq2.obj@ident == "CD14+ Monocytes")]
QUARTZseq.Monocytes <- names(QUARTZseq.obj@ident)[which(QUARTZseq.obj@ident == "CD14+ Monocytes")]
Dropseq.Monocytes <- names(Dropseq.obj@ident)[which(Dropseq.obj@ident == "CD14+ Monocytes")]
SCRBseq.Monocytes <- names(SCRBseq.obj@ident)[which(SCRBseq.obj@ident == "CD14+ Monocytes")]
X10Scilife.Monocytes <- names(X10Scilife.obj@ident)[which(X10Scilife.obj@ident == "CD14+ and FCGR3A+ Monocytes")]
X10Nuclei.Monocytes <- names(X10Nuclei.obj@ident)[which(X10Nuclei.obj@ident == "CD14+ Monocytes")]
ICELL8.Monocytes <- names(ICELL8.obj@ident)[which(ICELL8.obj@ident == "CD14+ and FCGR3A+ Monocytes")]
ddSEQ.Monocytes <- names(ddSEQ.obj@ident)[which(ddSEQ.obj@ident == "CD14+ and FCGR3A+ Monocytes")]
ddSEQexp1.Monocytes <- names(ddSEQexp1.obj@ident)[which(ddSEQexp1.obj@ident == "CD14+ and FCGR3A+ Monocytes")]
CB1.Monocytes <- names(CB1.obj@ident)[which(CB1.obj@ident == "CD14+ Monocytes")]

# End ####

# Taking out Bcells cells of each technique ####
MARSseq.Bcells <- names(MARSseq.obj@ident)[which(MARSseq.obj@ident == "B cells")]
CELseq2.Bcells <- names(CELseq2.obj@ident)[which(CELseq2.obj@ident == "B cells")]
QUARTZseq.Bcells <- names(QUARTZseq.obj@ident)[which(QUARTZseq.obj@ident == "B cells")]
Dropseq.Bcells <- names(Dropseq.obj@ident)[which(Dropseq.obj@ident == "B cells")]
SCRBseq.Bcells <- names(SCRBseq.obj@ident)[which(SCRBseq.obj@ident == "B cells")]
X10Scilife.Bcells <- names(X10Scilife.obj@ident)[which(X10Scilife.obj@ident == "B cells")]
X10Nuclei.Bcells <- names(X10Nuclei.obj@ident)[which(X10Nuclei.obj@ident == "B cells")]
ICELL8.Bcells <- names(ICELL8.obj@ident)[which(ICELL8.obj@ident == "B cells")]
ddSEQ.Bcells <- names(ddSEQ.obj@ident)[which(ddSEQ.obj@ident == "B cells")]
ddSEQexp1.Bcells <- names(ddSEQexp1.obj@ident)[which(ddSEQexp1.obj@ident == "B cells")]
CB1.Bcells <- names(CB1.obj@ident)[which(CB1.obj@ident == "B cells")]

# End ####

# Loding Downsampleded data ####
load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/MARSseq.hsap.full.SCE.jointDSmat.Robj")
MARSseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
MARSseq.DS.UMI <- MARSseq.DS$UMI
MARSseq.DS.Reads <- MARSseq.DS$Reads
rm(MARSseq.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/CELseq2.hsap.full.SCE.jointDSmat.Robj")
CELseq2.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
CELseq2.DS.UMI <- CELseq2.DS$UMI
CELseq2.DS.Reads <- CELseq2.DS$Reads
rm(CELseq2.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/QUARTZseq.hsap.full.SCE.jointDSmat.Robj")
QUARTZseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
QUARTZseq.DS.UMI <- QUARTZseq.DS$UMI
QUARTZseq.DS.Reads <- QUARTZseq.DS$Reads
rm(QUARTZseq.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Dropseq.hsap.full.SCE.jointDSmat.Robj")
Dropseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
Dropseq.DS.UMI <- Dropseq.DS$UMI
Dropseq.DS.Reads <- Dropseq.DS$Reads
rm(Dropseq.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/SCRBseq.hsap.full.SCE.jointDSmat.Robj")
SCRBseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
SCRBseq.DS.UMI <- SCRBseq.DS$UMI
SCRBseq.DS.Reads <- SCRBseq.DS$Reads
rm(SCRBseq.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/10XScilife.hsap.full.SCE.jointDSmat.Robj")
X10Scilife.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
X10Scilife.DS.UMI <- X10Scilife.DS$UMI
X10Scilife.DS.Reads <- X10Scilife.DS$Reads
rm(X10Scilife.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Nuclei10X.hsap.full.SCE.jointDSmat.Robj")
X10Nuclei.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
X10Nuclei.DS.UMI <- X10Nuclei.DS$UMI
X10Nuclei.DS.Reads <- X10Nuclei.DS$Reads
rm(X10Nuclei.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/ICELL8.hsap.full.SCE.jointDSmat.Robj")
ICELL8.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ICELL8.DS.UMI <- ICELL8.DS$UMI
ICELL8.DS.Reads <- ICELL8.DS$Reads


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/ddSEQexp1.hsap.full.SCE.jointDSmat.Robj")
ddSEQexp1.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ddSEQexp1.DS.UMI <- ddSEQexp1.DS$UMI
ddSEQexp1.DS.Reads <- ddSEQexp1.DS$Reads


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/ddSEQ.hsap.full.SCE.jointDSmat.Robj")
ddSEQ.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ddSEQ.DS.UMI <- ddSEQ.DS$UMI
ddSEQ.DS.Reads <- ddSEQ.DS$Reads

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/C1HT.hsap.full.SCE.jointDSmat.Robj")
C1HT.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
C1HT.DS.UMI <- C1HT.DS$UMI
C1HT.DS.Reads <- C1HT.DS$Reads

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/10X8x10K.hsap.full.SCE.jointDSmat.Robj")
X108x10.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
X108x10.DS.UMI <- X108x10.DS$UMI
X108x10.DS.Reads <- X108x10.DS$Reads

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/1CB.hsap.full.SCE.jointDSmat.Robj")
CB1.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
CB1.DS.UMI <- CB1.DS$UMI
CB1.DS.Reads <- CB1.DS$Reads

# End ####


#Plotting stepwise Downsampling for HEK ####

techniques <- c("MARSseq", "QUARTZseq", "CELseq2", "Dropseq", "SCRBseq", "X10Scilife", "X10Nuclei", "ICELL8", "ddSEQ", "ddSEQexp1", "CB1")#, "X108x10" , "ddSEQ", "ddSEQexp1", "C1HT"
DSth.df <- data.frame()
techs.HEK.20K.list <- list()
for (tech in techniques){
  print(tech)
  tech.DS.UMI <- get(paste(tech,".DS.UMI", sep = ""))
  tech.HEK <- get(paste(tech,".HEK", sep = ""))
  for (DSth in names(tech.DS.UMI)){
    print(paste(tech, DSth, sep = "_"))
    colnames(tech.DS.UMI[[DSth]]) <- gsub(x = colnames(tech.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    comm.cells <- intersect(tech.HEK, colnames(tech.DS.UMI[[DSth]]))
    DS.mat.HEKS <- tech.DS.UMI[[DSth]][, comm.cells]
    if (DSth == "downsampled_20000"){
      DS.mat.HEKS.dup <- DS.mat.HEKS
      DS.mat.HEKS.dup$gene_id <- rownames(DS.mat.HEKS.dup)
      techs.HEK.20K.list[[tech]] <- as.data.frame((DS.mat.HEKS.dup))}
    DS.gene.distribution.HEK <- colSums(DS.mat.HEKS[,]>0)
    DS.UMI.distribution.HEK <- colSums(DS.mat.HEKS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.HEK))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.HEK))
    DS.labels <- rep(DSth, length(DS.gene.distribution.HEK))
    DS.techs <- rep(tech, length(DS.gene.distribution.HEK))
    DS.df <- data.frame(nGenes = DS.gene.distribution.HEK, nUMIs = DS.UMI.distribution.HEK, DSthNum = DSth.number.vec, DStech = DS.techs)
    DSth.df <- rbind(DSth.df, DS.df)
    print(dim(DSth.df))
  }
}

MARSseq.raw.HEK <- techs.HEK.20K.list[["MARSseq"]]
MARSseq.raw.HEK <- mapIDs(MARSseq.raw.HEK, "hsap")
MARSseq.common.cells.HEK <- intersect(colnames(MARSseq.raw.HEK),colnames(MARSseq.obj@scale.data))
MARSseq.raw.HEK <- MARSseq.raw.HEK[rownames(MARSseq.obj@scale.data),MARSseq.common.cells.HEK]
MARSseq.raw.HEK <- as.data.frame(lapply(MARSseq.raw.HEK, as.integer))

CELseq2.raw.HEK <- techs.HEK.20K.list[["CELseq2"]]
CELseq2.raw.HEK <- mapIDs(CELseq2.raw.HEK, "hsap")
CELseq2.common.cells.HEK <- intersect(colnames(CELseq2.raw.HEK),colnames(CELseq2.obj@scale.data))
CELseq2.raw.HEK <- CELseq2.raw.HEK[rownames(CELseq2.obj@scale.data),CELseq2.common.cells.HEK]
CELseq2.raw.HEK <- as.data.frame(lapply(CELseq2.raw.HEK, as.integer))

QUARTZseq.raw.HEK <- techs.HEK.20K.list[["QUARTZseq"]]
QUARTZseq.raw.HEK <- mapIDs(QUARTZseq.raw.HEK, "hsap")
QUARTZseq.common.cells.HEK <- intersect(colnames(QUARTZseq.raw.HEK),colnames(QUARTZseq.obj@scale.data))
QUARTZseq.raw.HEK <- QUARTZseq.raw.HEK[rownames(QUARTZseq.obj@scale.data),QUARTZseq.common.cells.HEK]
QUARTZseq.raw.HEK <- as.data.frame(lapply(QUARTZseq.raw.HEK, as.integer))

SCRBseq.raw.HEK <- techs.HEK.20K.list[["SCRBseq"]]
SCRBseq.raw.HEK <- mapIDs(SCRBseq.raw.HEK, "hsap")
SCRBseq.common.cells.HEK <- intersect(colnames(SCRBseq.raw.HEK),colnames(SCRBseq.obj@scale.data))
SCRBseq.raw.HEK <- SCRBseq.raw.HEK[rownames(SCRBseq.obj@scale.data),SCRBseq.common.cells.HEK]
SCRBseq.raw.HEK <- as.data.frame(lapply(SCRBseq.raw.HEK, as.integer))

Dropseq.raw.HEK <- techs.HEK.20K.list[["Dropseq"]]
Dropseq.raw.HEK <- mapIDs(Dropseq.raw.HEK, "hsap")
Dropseq.common.cells.HEK <- intersect(colnames(Dropseq.raw.HEK),colnames(Dropseq.obj@scale.data))
Dropseq.raw.HEK <- Dropseq.raw.HEK[rownames(Dropseq.obj@scale.data),Dropseq.common.cells.HEK]
Dropseq.raw.HEK <- as.data.frame(lapply(Dropseq.raw.HEK, as.integer))

X10Scilife.raw.HEK <- techs.HEK.20K.list[["X10Scilife"]]
X10Scilife.raw.HEK <- mapIDs(X10Scilife.raw.HEK, "hsap")
X10Scilife.common.cells.HEK <- intersect(colnames(X10Scilife.raw.HEK),colnames(X10Scilife.obj@scale.data))
X10Scilife.raw.HEK <- X10Scilife.raw.HEK[rownames(X10Scilife.obj@scale.data),X10Scilife.common.cells.HEK]
X10Scilife.raw.HEK <- as.data.frame(lapply(X10Scilife.raw.HEK, as.integer))

X10Nuclei.raw.HEK <- techs.HEK.20K.list[["X10Nuclei"]]
X10Nuclei.raw.HEK <- mapIDs(X10Nuclei.raw.HEK, "hsap")
X10Nuclei.common.cells.HEK <- intersect(colnames(X10Nuclei.raw.HEK),colnames(X10Nuclei.obj@scale.data))
X10Nuclei.raw.HEK <- X10Nuclei.raw.HEK[rownames(X10Nuclei.obj@scale.data),X10Nuclei.common.cells.HEK]
X10Nuclei.raw.HEK <- as.data.frame(lapply(X10Nuclei.raw.HEK, as.integer))

ICELL8.raw.HEK <- techs.HEK.20K.list[["ICELL8"]]
ICELL8.raw.HEK <- mapIDs(ICELL8.raw.HEK, "hsap")
ICELL8.common.cells.HEK <- intersect(colnames(ICELL8.raw.HEK),colnames(ICELL8.obj@scale.data))
ICELL8.raw.HEK <- ICELL8.raw.HEK[rownames(ICELL8.obj@scale.data),ICELL8.common.cells.HEK]
ICELL8.raw.HEK <- as.data.frame(lapply(ICELL8.raw.HEK, as.integer))

library(scde)
#Running error models

#MARSseq error model
MARSseq.o.ifm <- scde.error.models(counts = MARSseq.raw.HEK, n.cores = 10, threshold.segmentation = TRUE, save.crossfit.plots = T, save.model.plots = T, verbose = 1)
save(MARSseq.o.ifm)
print("MARSseq model done")

#CELseq2 error model
CELseq2.o.ifm <- scde.error.models(counts = CELseq2.raw.HEK, n.cores = 10, threshold.segmentation = TRUE, save.crossfit.plots = T, save.model.plots = T, verbose = 1)
save(CELseq2.o.ifm)
print("CELseq2 model done")

#QUARTZseq error model
QUARTZseq.o.ifm <- scde.error.models(counts = QUARTZseq.raw.HEK, n.cores = 10, threshold.segmentation = TRUE, save.crossfit.plots = T, save.model.plots = T, verbose = 1)
save(QUARTZseq.o.ifm)
print("QUARTZseq model done")

#Dropseq error model
Dropseq.o.ifm <- scde.error.models(counts = Dropseq.raw.HEK, n.cores = 10, threshold.segmentation = TRUE, save.crossfit.plots = T, save.model.plots = T, verbose = 1)
save(Dropseq.o.ifm)
print("Dropseq model done")


#SCRBseq error model
SCRBseq.o.ifm <- scde.error.models(counts = SCRBseq.raw.HEK, n.cores = 10, threshold.segmentation = TRUE, save.crossfit.plots = T, save.model.plots = T, verbose = 1)
save(SCRBseq.o.ifm)
print("SCRBseq model done")

#X10Scilife error model
X10Scilife.o.ifm <- scde.error.models(counts = X10Scilife.raw.HEK, n.cores = 10, threshold.segmentation = TRUE, save.crossfit.plots = T, save.model.plots = T, verbose = 1)
save(X10Scilife.o.ifm)
print("X10Scilife model done")

#X10Nuclei error model
X10Nuclei.o.ifm <- scde.error.models(counts = X10Nuclei.raw.HEK, n.cores = 10, threshold.segmentation = TRUE, save.crossfit.plots = T, save.model.plots = T, verbose = 1)
save(X10Nuclei.o.ifm)
print("X10Nuclei model done")

#ICELL8 error model
ICELL8.o.ifm <- scde.error.models(counts = ICELL8.raw.HEK, n.cores = 10, threshold.segmentation = TRUE, save.crossfit.plots = T, save.model.plots = T, verbose = 1)
save(ICELL8.o.ifm)
print("ICELL8 model done")




#Running priors

#MARSseq priors
MARSseq.o.prior <- scde.expression.prior(models = MARSseq.o.ifm, counts = MARSseq.raw.HEK, length.out = 400, show.plot = FALSE)
save(MARSseq.o.prior)
print("MARSseq Prior done")

#CELseq2 priors
CELseq2.o.prior <- scde.expression.prior(models = CELseq2.o.ifm, counts = CELseq2.raw.HEK, length.out = 400, show.plot = FALSE)
save(CELseq2.o.prior)
print("CELseq2 Prior done")

#QUARTZseq priors
QUARTZseq.o.prior <- scde.expression.prior(models = QUARTZseq.o.ifm, counts = QUARTZseq.raw.HEK, length.out = 400, show.plot = FALSE)
save(QUARTZseq.o.prior)
print("QUARTZseq Prior done")

#Dropseq priors
Dropseq.o.prior <- scde.expression.prior(models = Dropseq.o.ifm, counts = Dropseq.raw.HEK, length.out = 400, show.plot = FALSE)
save(Dropseq.o.prior)
print("Dropseq Prior done")


#SCRBseq priors
SCRBseq.o.prior <- scde.expression.prior(models = SCRBseq.o.ifm, counts = SCRBseq.raw.HEK, length.out = 400, show.plot = FALSE)
save(SCRBseq.o.prior)
print("SCRBseq Prior done")

#X10Scilife priors
X10Scilife.o.prior <- scde.expression.prior(models = X10Scilife.o.ifm, counts = X10Scilife.raw.HEK, length.out = 400, show.plot = FALSE)
save(X10Scilife.o.prior)
print("X10Scilife Prior done")

#X10Nuclei priors
X10Nuclei.o.prior <- scde.expression.prior(models = X10Nuclei.o.ifm, counts = X10Nuclei.raw.HEK, length.out = 400, show.plot = FALSE)
save(X10Nuclei.o.prior)
print("X10Nuclei Prior done")

#ICELL8 priors
ICELL8.o.prior <- scde.expression.prior(models = ICELL8.o.ifm, counts = ICELL8.raw.HEK, length.out = 400, show.plot = FALSE)
save(ICELL8.o.prior)
print("ICELL8 Prior done")


#get self-fail probabilities (at a given observed count)

#MARSseq failure prob
MARSseq.self.fail <- scde.failure.probability(models = MARSseq.o.ifm, counts = MARSseq.raw.HEK)
save(MARSseq.self.fail)
print("MARSseq failure prob done")

#CELseq2 failure prob
CELseq2.self.fail <- scde.failure.probability(models = CELseq2.o.ifm, counts = CELseq2.raw.HEK)
save(CELseq2.self.fail)
print("CELseq2 failure prob done")

#QUARTZseq failure prob
QUARTZseq.self.fail <- scde.failure.probability(models = QUARTZseq.o.ifm, counts = QUARTZseq.raw.HEK)
save(QUARTZseq.self.fail)
print("QUARTZseq failure prob done")

#Dropseq failure prob
Dropseq.self.fail <- scde.failure.probability(models = Dropseq.o.ifm, counts = Dropseq.raw.HEK)
save(Dropseq.self.fail)
print("Dropseq failure prob done")


#SCRBseq failure prob
SCRBseq.self.fail <- scde.failure.probability(models = SCRBseq.o.ifm, counts = SCRBseq.raw.HEK)
save(SCRBseq.self.fail)
print("SCRBseq failure prob done")

#X10Scilife failure prob
X10Scilife.self.fail <- scde.failure.probability(models = X10Scilife.o.ifm, counts = X10Scilife.raw.HEK)
save(X10Scilife.self.fail)
print("X10Scilife failure prob done")

#X10Nuclei failure prob
X10Nuclei.self.fail <- scde.failure.probability(models = X10Nuclei.o.ifm, counts = X10Nuclei.raw.HEK)
save(X10Nuclei.self.fail)
print("X10Nuclei failure prob done")

#ICELL8 failure prob
ICELL8.self.fail <- scde.failure.probability(models = ICELL8.o.ifm, counts = ICELL8.raw.HEK)
save(ICELL8.self.fail)
print("ICELL8 failure prob done")






load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/gene_cl.ref.RData")
ref.markers <- gene_cl.ref


HEK.ref.markers <- ref.markers[["HEK cells"]]
HEK.markers.common <- Reduce(intersect, list(HEK.ref.markers, rownames(MARSseq.self.fail), rownames(CELseq2.self.fail), rownames(QUARTZseq.self.fail), rownames(SCRBseq.self.fail), rownames(Dropseq.self.fail), rownames(X10Scilife.self.fail), rownames(x10Nuclei.self.fail), rownames(ICELL8.self.fail)))


techniques <- c("MARSseq", "QUARTZseq", "CELseq2", "Dropseq", "SCRBseq",  "ICELL8","X10Scilife", "Nuclei10X")#, "X108x10")
dropout.df <- data.frame()
for (tech in techniques){
  print(tech)
  tech.fail.mat <- get(paste(tech,".self.fail", sep = ""))
  tech.fail.mat.markers <- tech.fail.mat[HEK.markers.common,]
  tech.dropout.prob <- rowSums(tech.fail.mat.markers)
  dropout.df <- rbind(dropout.df, data.frame(dorpout= tech.dropout.prob, tech= rep(tech, length(tech.dropout.prob))))
}
pdf("dropout_boxplots.pdf")
ggplot(data=dropout.df, aes(x=tech, y=dorpout, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
dev.off()



