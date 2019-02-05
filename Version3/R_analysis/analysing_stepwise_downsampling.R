library(scater)
library(ggplot2)
library(biomaRt)
library(tibble)
library(Seurat)
library(data.table)

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


# End ####


#Plotting stepwise Downsampling for HEK ####

techniques <- c("MARSseq", "QUARTZseq", "CELseq2", "Dropseq", "SCRBseq", "X10Scilife", "X10Nuclei", "ICELL8", "ddSEQ", "ddSEQexp1")#, "CB1")#, "X108x10" , "ddSEQ", "ddSEQexp1", "C1HT"
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
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwide_DS_analysis/all_techs_stepwise_DS_HEK_plots_V2.pdf")
ggplot(DSth.df, aes(x=DSthNum, y=nGenes, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))
ggplot(DSth.df, aes(x=DSthNum, y=nUMIs, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))

#ggplot(DSth.df, aes(x=DSthNum, y=nGenes, fill = DSth)) +geom_boxplot() + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(group = 1))

ggplot(data=DSth.df, aes(x=DStech, y=nGenes, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
ggplot(data=DSth.df, aes(x=DStech, y=nUMIs, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
dev.off()

print("Stepwise DS plots for HEK done")

#Plotting stepwise Downsampling for Monocytes ####

techniques <- c("MARSseq", "QUARTZseq", "CELseq2", "Dropseq", "SCRBseq", "X10Scilife", "X10Nuclei", "ICELL8", "ddSEQ", "ddSEQexp1")#, "CB1")#, "X108x10" , "ddSEQ", "ddSEQexp1", "C1HT"
DSth.df <- data.frame()
techs.Monocytes.20K.list <- list()
for (tech in techniques){
  print(tech)
  tech.DS.UMI <- get(paste(tech,".DS.UMI", sep = ""))
  tech.Monocytes <- get(paste(tech,".Monocytes", sep = ""))
  for (DSth in names(tech.DS.UMI)){
    print(paste(tech, DSth, sep = "_"))
    colnames(tech.DS.UMI[[DSth]]) <- gsub(x = colnames(tech.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    comm.cells <- intersect(tech.Monocytes, colnames(tech.DS.UMI[[DSth]]))
    DS.mat.MonocytesS <- tech.DS.UMI[[DSth]][, comm.cells]
    if (DSth == "downsampled_20000"){
      DS.mat.MonocytesS.dup <- DS.mat.MonocytesS
      DS.mat.MonocytesS.dup$gene_id <- rownames(DS.mat.MonocytesS.dup)
      techs.Monocytes.20K.list[[tech]] <- as.data.frame((DS.mat.MonocytesS.dup))}
    DS.gene.distribution.Monocytes <- colSums(DS.mat.MonocytesS[,]>0)
    DS.UMI.distribution.Monocytes <- colSums(DS.mat.MonocytesS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.Monocytes))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.Monocytes))
    DS.labels <- rep(DSth, length(DS.gene.distribution.Monocytes))
    DS.techs <- rep(tech, length(DS.gene.distribution.Monocytes))
    DS.df <- data.frame(nGenes = DS.gene.distribution.Monocytes, nUMIs = DS.UMI.distribution.Monocytes, DSthNum = DSth.number.vec, DStech = DS.techs)
    DSth.df <- rbind(DSth.df, DS.df)
    print(dim(DSth.df))
  }
}
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwide_DS_analysis/all_techs_stepwise_DS_Monocytes_plots_V2.pdf")
ggplot(DSth.df, aes(x=DSthNum, y=nGenes, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))
ggplot(DSth.df, aes(x=DSthNum, y=nUMIs, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))

#ggplot(DSth.df, aes(x=DSthNum, y=nGenes, fill = DSth)) +geom_boxplot() + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(group = 1))

ggplot(data=DSth.df, aes(x=DStech, y=nGenes, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
ggplot(data=DSth.df, aes(x=DStech, y=nUMIs, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
dev.off()

print("Stepwise DS plots for Monocytes done")

#Plotting stepwise Downsampling for Bcell ####

techniques <- c("MARSseq", "QUARTZseq", "CELseq2", "Dropseq", "SCRBseq", "X10Scilife", "X10Nuclei", "ICELL8", "ddSEQ", "ddSEQexp1")#, "CB1")#, "X108x10" , "ddSEQ", "ddSEQexp1", "C1HT"
DSth.df <- data.frame()
techs.Bcells.10K.list <- list()
for (tech in techniques){
  print(tech)
  tech.DS.UMI <- get(paste(tech,".DS.UMI", sep = ""))
  tech.Bcells <- get(paste(tech,".Bcells", sep = ""))
  for (DSth in names(tech.DS.UMI)){
    print(paste(tech, DSth, sep = "_"))
    colnames(tech.DS.UMI[[DSth]]) <- gsub(x = colnames(tech.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    comm.cells <- intersect(tech.Bcells, colnames(tech.DS.UMI[[DSth]]))
    DS.mat.BcellsS <- tech.DS.UMI[[DSth]][, comm.cells]
    if (DSth == "downsampled_10000"){
      DS.mat.BcellsS.dup <- DS.mat.BcellsS
      DS.mat.BcellsS.dup$gene_id <- rownames(DS.mat.BcellsS.dup)
      techs.Bcells.10K.list[[tech]] <- as.data.frame((DS.mat.BcellsS.dup))}
    DS.gene.distribution.Bcells <- colSums(DS.mat.BcellsS[,]>0)
    DS.UMI.distribution.Bcells <- colSums(DS.mat.BcellsS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.Bcells))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.Bcells))
    DS.labels <- rep(DSth, length(DS.gene.distribution.Bcells))
    DS.techs <- rep(tech, length(DS.gene.distribution.Bcells))
    DS.df <- data.frame(nGenes = DS.gene.distribution.Bcells, nUMIs = DS.UMI.distribution.Bcells, DSthNum = DSth.number.vec, DStech = DS.techs)
    DSth.df <- rbind(DSth.df, DS.df)
    print(dim(DSth.df))
  }
}
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwide_DS_analysis/all_techs_stepwise_DS_Bcells_plots_V2.pdf")
ggplot(DSth.df, aes(x=DSthNum, y=nGenes, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))
ggplot(DSth.df, aes(x=DSthNum, y=nUMIs, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))

#ggplot(DSth.df, aes(x=DSthNum, y=nGenes, fill = DSth)) +geom_boxplot() + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(group = 1))

ggplot(data=DSth.df, aes(x=DStech, y=nGenes, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
ggplot(data=DSth.df, aes(x=DStech, y=nUMIs, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
dev.off()

print("Stepwise DS plots for Bcells done")






# The PCA for celltypes HEK
print("PCA analysis started")
library(plyr)
techs.HEK.20K.df=join_all(techs.HEK.20K.list, by = "gene_id", type = 'full')
rownames(techs.HEK.20K.df) <- techs.HEK.20K.df[,"gene_id"]
gene.id.col <- which(colnames(techs.HEK.20K.df) == "gene_id")
techs.HEK.20K.df <- techs.HEK.20K.df[,-gene.id.col]
dim(techs.HEK.20K.df)
techs.HEK.20K.df <- na.omit(techs.HEK.20K.df)
#techs.HEK.20K.log.df <- log(techs.HEK.20K.df + 1)
#pca.res <- prcomp(techs.HEK.20K.log.df, scale. = T)
#pca.df <- as.data.frame(pca.res$rotation[,c("PC1","PC2")])
#tech.vec <- sapply(strsplit(rownames(pca.df), split = "_", fixed = T), `[`,1)
#pca.df <- cbind(pca.df,tech.vec)
#pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwide_DS_analysis/all_techs_stepwise_DS20K_PCA_HEK.pdf")
#ggplot(pca.df, aes(x= PC1, y = PC2, group= tech.vec)) + geom_point(aes(color = tech.vec))
#dev.off()
HEK.seurat <- CreateSeuratObject(techs.HEK.20K.df, min.cells = 5, min.genes = 1, normalization.method = "LogNormalize", scale.factor = 10000)#, total.expr = 1e4, project = "Dicroce")
HEK.seurat <- ScaleData(object= HEK.seurat)
HEK.seurat <- FindVariableGenes(object = HEK.seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.2, x.high.cutoff = 6, y.cutoff = 0.5)
print("HEK HVG:")
print(dim(HEK.seurat@scale.data))
print(length(HEK.seurat@var.genes))
HEK.seurat <- RunPCA(HEK.seurat, pc.genes = HEK.seurat@var.genes, do.print = F)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwide_DS_analysis/all_techs_stepwise_DS20K_PCAseurat_HEK_V2.pdf")
PCAPlot(HEK.seurat, 1,2, group.by = "orig.ident")
dev.off()
print("PCA DS plots for HEK done")


# The PCA for celltypes Monocytes
print("PCA analysis started")
library(plyr)
techs.Monocytes.20K.df=join_all(techs.Monocytes.20K.list, by = "gene_id", type = 'full')
rownames(techs.Monocytes.20K.df) <- techs.Monocytes.20K.df[,"gene_id"]
gene.id.col <- which(colnames(techs.Monocytes.20K.df) == "gene_id")
techs.Monocytes.20K.df <- techs.Monocytes.20K.df[,-gene.id.col]
dim(techs.Monocytes.20K.df)
techs.Monocytes.20K.df <- na.omit(techs.Monocytes.20K.df)
#techs.Monocytes.20K.log.df <- log(techs.Monocytes.20K.df + 1)
#pca.res <- prcomp(techs.Monocytes.20K.log.df, scale. = T)
#pca.df <- as.data.frame(pca.res$rotation[,c("PC1","PC2")])
#tech.vec <- sapply(strsplit(rownames(pca.df), split = "_", fixed = T), `[`,1)
#pca.df <- cbind(pca.df,tech.vec)
#pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwide_DS_analysis/all_techs_stepwise_DS20K_PCA_Monocytes.pdf")
#ggplot(pca.df, aes(x= PC1, y = PC2, group= tech.vec)) + geom_point(aes(color = tech.vec))
#dev.off()
Monocytes.seurat <- CreateSeuratObject(techs.Monocytes.20K.df, min.cells = 5, min.genes = 1, normalization.method = "LogNormalize", scale.factor = 10000)#, total.expr = 1e4, project = "Dicroce")
Monocytes.seurat <- ScaleData(object= Monocytes.seurat)
Monocytes.seurat <- FindVariableGenes(object = Monocytes.seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.2, x.high.cutoff = 6, y.cutoff = 0.5)
print("Monocytes HVG:")
print(dim(Monocytes.seurat@scale.data))
print(length(Monocytes.seurat@var.genes))
Monocytes.seurat <- RunPCA(Monocytes.seurat, pc.genes = Monocytes.seurat@var.genes, do.print = F)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwide_DS_analysis/all_techs_stepwise_DS20K_PCAseurat_Monocytes_V2.pdf")
PCAPlot(Monocytes.seurat, 1,2, group.by = "orig.ident")
dev.off()
print("PCA DS plots for Monocytes done")


# The PCA for celltypes Bcells
print("PCA analysis started")
library(plyr)
techs.Bcells.10K.df=join_all(techs.Bcells.10K.list, by = "gene_id", type = 'full')
rownames(techs.Bcells.10K.df) <- techs.Bcells.10K.df[,"gene_id"]
gene.id.col <- which(colnames(techs.Bcells.10K.df) == "gene_id")
techs.Bcells.10K.df <- techs.Bcells.10K.df[,-gene.id.col]
dim(techs.Bcells.10K.df)
techs.Bcells.10K.df <- na.omit(techs.Bcells.10K.df)
#techs.Bcells.10K.log.df <- log(techs.Bcells.10K.df + 1)
#pca.res <- prcomp(techs.Bcells.10K.log.df, scale. = T)
#pca.df <- as.data.frame(pca.res$rotation[,c("PC1","PC2")])
#tech.vec <- sapply(strsplit(rownames(pca.df), split = "_", fixed = T), `[`,1)
#pca.df <- cbind(pca.df,tech.vec)
#pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwide_DS_analysis/all_techs_stepwise_DS10K_PCA_Bcells.pdf")
#ggplot(pca.df, aes(x= PC1, y = PC2, group= tech.vec)) + geom_point(aes(color = tech.vec))
#dev.off()
Bcells.seurat <- CreateSeuratObject(techs.Bcells.10K.df, min.cells = 5, min.genes = 1, normalization.method = "LogNormalize", scale.factor = 10000)#, total.expr = 1e4, project = "Dicroce")
Bcells.seurat <- ScaleData(object= Bcells.seurat)
Bcells.seurat <- FindVariableGenes(object = Bcells.seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.2, x.high.cutoff = 6, y.cutoff = 0.5)
print("Bcells HVG:")
print(dim(Bcells.seurat@scale.data))
print(length(Bcells.seurat@var.genes))
Bcells.seurat <- RunPCA(Bcells.seurat, pc.genes = Bcells.seurat@var.genes, do.print = F)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwide_DS_analysis/all_techs_stepwise_DS20K_PCAseurat_Bcells_V2.pdf")
PCAPlot(Bcells.seurat, 1,2, group.by = "orig.ident")
dev.off()
print("PCA DS plots for Bcells done")