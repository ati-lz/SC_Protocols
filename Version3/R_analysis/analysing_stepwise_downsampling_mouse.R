library(scater)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(biomaRt)
library(tibble)
library(Seurat)
library(data.table)

# loading Annotated seurat objects ####

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final_rename/Mouse/celseq_seu.obj.RData")
CELseq2.mmus.obj <- celseq
rm(celseq)
CELseq2.mmus.metadata <- CELseq2.mmus.obj@meta.data


load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final_rename/Mouse/marsseq_seu.obj.RData")
MARSseq.mmus.obj <- marsseq
rm(marsseq)
MARSseq.mmus.metadata <- MARSseq.mmus.obj@meta.data


load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final_rename/Mouse/quartzseq_seu.obj.RData")
QUARTZseq.mmus.obj <- quartzseq
rm(quartzseq)
QUARTZseq.mmus.metadata <- QUARTZseq.mmus.obj@meta.data

#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final_rename/Mouse/scrbseq_seu.obj.RData")
#SCRBseq.mmus.obj <- scrbseq
#rm(scrbseq)
#SCRBseq.mmus.metadata <- SCRBseq.mmus.obj@meta.data


load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final_rename/Mouse/smartseq_seu.obj.RData")
SMARTseq2.mmus.obj <- smartseq
rm(smartseq)
SMARTseq2.mmus.metadata <- SMARTseq2.mmus.obj@meta.data


load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final_rename/Mouse/c1ht.s_seu.obj.RData")
C1HTsmall.mmus.obj <- c1ht.s
rm(c1ht.s)
C1HTsmall.mmus.metadata <- C1HTsmall.mmus.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final_rename/Mouse/c1ht.m_seu.obj.RData")
C1HTmedium.mmus.obj <- c1ht.m
rm(c1ht.m)
C1HTmedium.mmus.metadata <- C1HTmedium.mmus.obj@meta.data


load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final_rename/Mouse/chromium_seu.obj.RData")
Chromium.mmus.obj <- chromium
rm(chromium)
Chromium.mmus.metadata <- Chromium.mmus.obj@meta.data


load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final_rename/Mouse/nuclei_seu.obj.RData")
ChromiumNuclei.mmus.obj <- nuclei
rm(nuclei)
ChromiumNuclei.mmus.metadata <- ChromiumNuclei.mmus.obj@meta.data


load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final_rename/Mouse/ddseq_seu.obj.RData")
ddSEQ.mmus.obj <- ddseq
rm(ddseq)
ddSEQ.mmus.metadata <- ddSEQ.mmus.obj@meta.data


load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final_rename/Mouse/dropseq_seu.obj.RData")
Dropseq.mmus.obj <- dropseq
rm(dropseq)
Dropseq.mmus.metadata <- Dropseq.mmus.obj@meta.data


#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final_rename/Mouse/icell8_seu.obj.RData")
#ICELL8.mmus.obj <- icell8
#rm(icell8)
#ICELL8.mmus.metadata <- ICELL8.mmus.obj@meta.data


load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final_rename/Mouse/indrop_seu.obj.RData")
inDrop.mmus.obj <- indrop
rm(indrop)
inDrop.mmus.metadata <- inDrop.mmus.obj@meta.data


# End ####

print("Seurat reading DONE")
# Taking out Secretory cells of each technique ####
CELseq2.Secretory <- rownames(CELseq2.mmus.metadata)[which(CELseq2.mmus.metadata$clean.id == "Secretory")]
MARSseq.Secretory <- rownames(MARSseq.mmus.metadata)[which(MARSseq.mmus.metadata$clean.id == "Secretory")]
QUARTZseq.Secretory <- rownames(QUARTZseq.mmus.metadata)[which(QUARTZseq.mmus.metadata$clean.id == "Secretory")]
#SCRBseq.Secretory <- rownames(SCRBseq.mmus.metadata)[which(SCRBseq.mmus.metadata$clean.id == "Secretory")]
SMARTseq2.Secretory <- rownames(SMARTseq2.mmus.metadata)[which(SMARTseq2.mmus.metadata$clean.id == "Secretory")]
C1HTsmall.Secretory <- rownames(C1HTsmall.mmus.metadata[which(C1HTsmall.mmus.metadata$clean.id == "Secretory"),])
C1HTmedium.Secretory <- rownames(C1HTmedium.mmus.metadata[which(C1HTmedium.mmus.metadata$clean.id == "Secretory"),])
Chromium.Secretory <- rownames(Chromium.mmus.metadata[which(Chromium.mmus.metadata$clean.id == "Secretory"),])
ChromiumNuclei.Secretory <- rownames(ChromiumNuclei.mmus.metadata[which(ChromiumNuclei.mmus.metadata$clean.id == "Secretory"),])
ddSEQ.Secretory <- rownames(ddSEQ.mmus.metadata[which(ddSEQ.mmus.metadata$clean.id == "Secretory"),])
Dropseq.Secretory <- rownames(Dropseq.mmus.metadata[which(Dropseq.mmus.metadata$clean.id == "Secretory"),])
#ICELL8.Secretory <- rownames(ICELL8.mmus.metadata[which(ICELL8.mmus.metadata$clean.id == "Secretory"),])
inDrop.Secretory <- rownames(inDrop.mmus.metadata[which(inDrop.mmus.metadata$clean.id == "Secretory"),])

# End ####

# Transit Amplifying
CELseq2.TA <- rownames(CELseq2.mmus.metadata)[which(CELseq2.mmus.metadata$clean.id == "Transit Amplifying")]
MARSseq.TA <- rownames(MARSseq.mmus.metadata)[which(MARSseq.mmus.metadata$clean.id == "Transit Amplifying")]
QUARTZseq.TA <- rownames(QUARTZseq.mmus.metadata)[which(QUARTZseq.mmus.metadata$clean.id == "Transit Amplifying")]
#SCRBseq.TA <- rownames(SCRBseq.mmus.metadata)[which(SCRBseq.mmus.metadata$clean.id == "Transit Amplifying")]
SMARTseq2.TA <- rownames(SMARTseq2.mmus.metadata)[which(SMARTseq2.mmus.metadata$clean.id == "Transit Amplifying")]
C1HTsmall.TA <- rownames(C1HTsmall.mmus.metadata[which(C1HTsmall.mmus.metadata$clean.id == "Transit Amplifying"),])
C1HTmedium.TA <- rownames(C1HTmedium.mmus.metadata[which(C1HTmedium.mmus.metadata$clean.id == "Transit Amplifying"),])
Chromium.TA <- rownames(Chromium.mmus.metadata[which(Chromium.mmus.metadata$clean.id == "Transit Amplifying"),])
ChromiumNuclei.TA <- rownames(ChromiumNuclei.mmus.metadata[which(ChromiumNuclei.mmus.metadata$clean.id == "Transit Amplifying"),])
ddSEQ.TA <- rownames(ddSEQ.mmus.metadata[which(ddSEQ.mmus.metadata$clean.id == "Transit Amplifying"),])
Dropseq.TA <- rownames(Dropseq.mmus.metadata[which(Dropseq.mmus.metadata$clean.id == "Transit Amplifying"),])
#ICELL8.TA <- rownames(ICELL8.mmus.metadata[which(ICELL8.mmus.metadata$clean.id == "Transit Amplifying"),])
inDrop.TA <- rownames(inDrop.mmus.metadata[which(inDrop.mmus.metadata$clean.id == "Transit Amplifying"),])

# End ####


print("celltype cells DONE")
# Loding Downsampleded data ####
load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Mouse/CELseq2.mmus.full.SCE.jointDSmat.Robj")
CELseq2.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
CELseq2.DS.UMI <- CELseq2.DS$UMI
#CELseq2.DS.Reads <- CELseq2.DS$Reads
rm(CELseq2.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Mouse/MARSseq.mmus.full.SCE.jointDSmat.Robj")
MARSseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
MARSseq.DS.UMI <- MARSseq.DS$UMI
#MARSseq.DS.Reads <- MARSseq.DS$Reads
rm(MARSseq.DS)


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Mouse/QUARTZseq.mmus.full.SCE.jointDSmat.Robj")
QUARTZseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
QUARTZseq.DS.UMI <- QUARTZseq.DS$UMI
#QUARTZseq.DS.Reads <- QUARTZseq.DS$Reads
rm(QUARTZseq.DS)

#load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Mouse/SCRBseq.mmus.full.SCE.jointDSmat.Robj")
#SCRBseq.DS <- output.readcount.umicount.joint.mats
#rm(output.readcount.umicount.joint.mats)
#SCRBseq.DS.UMI <- SCRBseq.DS$UMI
#SCRBseq.DS.Reads <- SCRBseq.DS$Reads
#rm(SCRBseq.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Mouse/SMARTseqFINAL.mmus.full.SCE.jointDSmat.Robj")
SMARTseq2.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
#SMARTseq2.DS.UMI <- SMARTseq2.DS$UMI ### HERE we dont have UMI, so we put Reads instead of UMIs for the name
SMARTseq2.DS.UMI <- SMARTseq2.DS$Reads
rm(SMARTseq2.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Mouse/C1HTsmall.mmus.full.SCE.jointDSmat.Robj")
C1HTsmall.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
C1HTsmall.DS.UMI <- C1HTsmall.DS$UMI
#C1HTsmall.DS.Reads <- C1HTsmall.DS$Reads
rm(C1HTsmall.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Mouse/C1HTmedium.mmus.full.SCE.jointDSmat.Robj")
C1HTmedium.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
C1HTmedium.DS.UMI <- C1HTmedium.DS$UMI
#C1HTmedium.DS.Reads <- C1HTmedium.DS$Reads
rm(C1HTmedium.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Mouse/10X2x5K.mmus.full.SCE.jointDSmat.Robj")
Chromium.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
Chromium.DS.UMI <- Chromium.DS$UMI
#Chromium.DS.Reads <- Chromium.DS$Reads
rm(Chromium.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Mouse/Nuclei10X.mmus.full.SCE.jointDSmat.Robj")
ChromiumNuclei.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ChromiumNuclei.DS.UMI <- ChromiumNuclei.DS$UMI
#ChromiumNuclei.DS.Reads <- ChromiumNuclei.DS$Reads
rm(ChromiumNuclei.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Mouse/ddSEQ.mmus.full.SCE.jointDSmat.Robj")
ddSEQ.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ddSEQ.DS.UMI <- ddSEQ.DS$UMI
#ddSEQ.DS.Reads <- ddSEQ.DS$Reads
rm(ddSEQ.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Mouse/Dropseq.mmus.full.SCE.jointDSmat.Robj")
Dropseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
Dropseq.DS.UMI <- Dropseq.DS$UMI
#Dropseq.DS.Reads <- Dropseq.DS$Reads
rm(Dropseq.DS)

#load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Mouse/ICELL8.mmus.full.SCE.jointDSmat.Robj")
#ICELL8.DS <- output.readcount.umicount.joint.mats
#rm(output.readcount.umicount.joint.mats)
#ICELL8.DS.UMI <- ICELL8.DS$UMI
#ICELL8.DS.Reads <- ICELL8.DS$Reads
#rm(ICELL8.DS)


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Mouse/1CB.mmus.full.SCE.jointDSmat.Robj")
inDrop.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
inDrop.DS.UMI <- inDrop.DS$UMI
#inDrop.DS.Reads <- inDrop.DS$Reads
rm(inDrop.DS)

# End ####

print("Downsampleded reading DONE")
#Plotting stepwise Downsampling for Secretory ####

techniques <- c("CELseq2","MARSseq","QUARTZseq","SMARTseq2", "C1HTsmall","C1HTmedium","Chromium", "ChromiumNuclei","ddSEQ","Dropseq","inDrop") #,"SCRBseq","ICELL8"
DSth.df.Secretory <- data.frame()
techs.Secretory.20K.list <- list()
for (tech in techniques){
  print(tech)
  tech.DS.UMI <- get(paste(tech,".DS.UMI", sep = ""))
  tech.Secretory <- get(paste(tech,".Secretory", sep = ""))
  for (DSth in names(tech.DS.UMI)){
    print(paste(tech, DSth, sep = "_"))
    colnames(tech.DS.UMI[[DSth]]) <- gsub(x = colnames(tech.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    comm.cells <- intersect(tech.Secretory, colnames(tech.DS.UMI[[DSth]]))
    DS.mat.Secretory <- tech.DS.UMI[[DSth]][, comm.cells]
    if (DSth == "downsampled_20000"){
      DS.mat.Secretory.dup <- DS.mat.Secretory
      DS.mat.Secretory.dup$gene_id <- rownames(DS.mat.Secretory.dup)
      techs.Secretory.20K.list[[tech]] <- as.data.frame((DS.mat.Secretory.dup))}
    DS.gene.distribution.Secretory <- colSums(DS.mat.Secretory[,]>0)
    DS.UMI.distribution.Secretory <- colSums(DS.mat.Secretory)
    DS.labels <- rep(DSth, length(DS.gene.distribution.Secretory))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.Secretory))
    DS.labels <- rep(DSth, length(DS.gene.distribution.Secretory))
    DS.techs <- rep(tech, length(DS.gene.distribution.Secretory))
    DS.df <- data.frame(nGenes = DS.gene.distribution.Secretory, nUMIs = DS.UMI.distribution.Secretory, DSthNum = DSth.number.vec, DStech = DS.techs)
    DSth.df.Secretory <- rbind(DSth.df.Secretory, DS.df)
    print(dim(DSth.df.Secretory))
  }
}

save(DSth.df.Secretory, file="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/Stepwise_DS_plotdf_Secretory.RData")
print("Stepwise DS plots for Secretory done")

#Plotting stepwise Downsampling for TA ####

DSth.df.TA <- data.frame()
techs.TA.20K.list <- list()
for (tech in techniques){
  print(tech)
  tech.DS.UMI <- get(paste(tech,".DS.UMI", sep = ""))
  tech.TA <- get(paste(tech,".TA", sep = ""))
  for (DSth in names(tech.DS.UMI)){
    print(paste(tech, DSth, sep = "_"))
    colnames(tech.DS.UMI[[DSth]]) <- gsub(x = colnames(tech.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    comm.cells <- intersect(tech.TA, colnames(tech.DS.UMI[[DSth]]))
    DS.mat.TAS <- tech.DS.UMI[[DSth]][, comm.cells]
    if (DSth == "downsampled_20000"){
      DS.mat.TAS.dup <- DS.mat.TAS
      DS.mat.TAS.dup$gene_id <- rownames(DS.mat.TAS.dup)
      techs.TA.20K.list[[tech]] <- as.data.frame((DS.mat.TAS.dup))}
    DS.gene.distribution.TA <- colSums(DS.mat.TAS[,]>0)
    DS.UMI.distribution.TA <- colSums(DS.mat.TAS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.TA))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.TA))
    DS.labels <- rep(DSth, length(DS.gene.distribution.TA))
    DS.techs <- rep(tech, length(DS.gene.distribution.TA))
    DS.df <- data.frame(nGenes = DS.gene.distribution.TA, nUMIs = DS.UMI.distribution.TA, DSthNum = DSth.number.vec, DStech = DS.techs)
    DSth.df.TA <- rbind(DSth.df.TA, DS.df)
    print(dim(DSth.df.TA))
  }
}

save(DSth.df.TA, file="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/Stepwise_DS_plotdf_TA.RData")

print("Stepwise DS plots for TA done")



# The PCA for celltypes Secretory
print("PCA analysis started")
library(plyr)
techs.Secretory.20K.df=join_all(techs.Secretory.20K.list, by = "gene_id", type = 'full')
rownames(techs.Secretory.20K.df) <- techs.Secretory.20K.df[,"gene_id"]
gene.id.col <- which(colnames(techs.Secretory.20K.df) == "gene_id")
techs.Secretory.20K.df <- techs.Secretory.20K.df[,-gene.id.col]
dim(techs.Secretory.20K.df)
techs.Secretory.20K.df <- na.omit(techs.Secretory.20K.df)
#techs.Secretory.20K.log.df <- log(techs.Secretory.20K.df + 1)
#pca.res <- prcomp(techs.Secretory.20K.log.df, scale. = T)
#pca.df <- as.data.frame(pca.res$rotation[,c("PC1","PC2")])
#tech.vec <- sapply(strsplit(rownames(pca.df), split = "_", fixed = T), `[`,1)
#pca.df <- cbind(pca.df,tech.vec)
#pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCA_Secretory.pdf")
#ggplot(pca.df, aes(x= PC1, y = PC2, group= tech.vec)) + geom_point(aes(color = tech.vec))
#dev.off()
Secretory.seurat <- CreateSeuratObject(techs.Secretory.20K.df, min.cells = 5, min.genes = 1, normalization.method = "LogNormalize", scale.factor = 10000)#, total.expr = 1e4, project = "Dicroce")
Secretory.seurat <- ScaleData(object= Secretory.seurat)
Secretory.seurat <- FindVariableGenes(object = Secretory.seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.2, x.high.cutoff = 6, y.cutoff = 0.5)
print("Secretory HVG:")
print(dim(Secretory.seurat@scale.data))
print(length(Secretory.seurat@var.genes))
Secretory.seurat <- RunPCA(Secretory.seurat, pc.genes = Secretory.seurat@var.genes, do.print = F)
Secretory.colors.nUMI = Secretory.seurat@meta.data[names(Secretory.seurat@ident), "nUMI"]
names(Secretory.colors.nUMI) <- names(Secretory.seurat@ident)
Secretory.data.plot <- as.data.frame(Secretory.seurat@dr$pca@cell.embeddings[,1:8])
Secretory.data.plot <- cbind(Secretory.data.plot, nUMI=Secretory.colors.nUMI, techs=Secretory.seurat@meta.data$orig.ident)

p0 <- PCAPlot(Secretory.seurat, 1,2, group.by = "orig.ident")
p00 <- ggplot(data = Secretory.data.plot, mapping = aes(x = PC1, y = PC2, color=Secretory.colors.nUMI)) + 
  geom_point(size = 1,shape = 16)+ scale_color_gradient( low = "grey", high = "red")

p1 <- ggplot(data = Secretory.data.plot, mapping = aes(x = PC1, y = PC2, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p11 <- ggMarginal(p1, data = Secretory.data.plot, x= Secretory.colors.nUMI, y= Secretory.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

p2 <- ggplot(data = Secretory.data.plot, mapping = aes(x = PC3, y = PC4, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p22 <- ggMarginal(p2, data = Secretory.data.plot, x= Secretory.colors.nUMI, y= Secretory.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

p3 <- ggplot(data = Secretory.data.plot, mapping = aes(x = PC5, y = PC6, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p33 <- ggMarginal(p3, data = Secretory.data.plot, x= Secretory.colors.nUMI, y= Secretory.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

p4 <- ggplot(data = Secretory.data.plot, mapping = aes(x = PC7, y = PC8, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p44 <- ggMarginal(p4, data = Secretory.data.plot, x= Secretory.colors.nUMI, y= Secretory.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

Secretory.plot.list <- list(p0,p00,p1,p11,p2,p22,p3,p33,p4,p44)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCAseurat_Secretory_V6.pdf")
invisible(lapply(Secretory.plot.list, print))
dev.off()
save(Secretory.data.plot, file="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCAseurat_Secretory_dataPlot.RData")
save(Secretory.seurat, file="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCAseurat_Secretory_SeuratObj.RData")

#Secretory Cluster tree
Secretory.seurat@ident <- Secretory.seurat@meta.data$orig.ident
names(Secretory.seurat@ident) <- rownames(Secretory.seurat@meta.data)
Secretory.seurat <- BuildClusterTree(Secretory.seurat, genes.use = Secretory.seurat@var.genes, pcs.use = 1:10)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_Secretory_ClusterTree_pc1to10.pdf")
plotClusterTree(Secretory.seurat)
dev.off()

print("PCA DS plots for Secretory done")


# The PCA for celltypes TA
print("PCA analysis started")
library(plyr)
techs.TA.20K.df=join_all(techs.TA.20K.list, by = "gene_id", type = 'full')
rownames(techs.TA.20K.df) <- techs.TA.20K.df[,"gene_id"]
gene.id.col <- which(colnames(techs.TA.20K.df) == "gene_id")
techs.TA.20K.df <- techs.TA.20K.df[,-gene.id.col]
dim(techs.TA.20K.df)
techs.TA.20K.df <- na.omit(techs.TA.20K.df)
#techs.TA.20K.log.df <- log(techs.TA.20K.df + 1)
#pca.res <- prcomp(techs.TA.20K.log.df, scale. = T)
#pca.df <- as.data.frame(pca.res$rotation[,c("PC1","PC2")])
#tech.vec <- sapply(strsplit(rownames(pca.df), split = "_", fixed = T), `[`,1)
#pca.df <- cbind(pca.df,tech.vec)
#pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCA_TA.pdf")
#ggplot(pca.df, aes(x= PC1, y = PC2, group= tech.vec)) + geom_point(aes(color = tech.vec))
#dev.off()
TA.seurat <- CreateSeuratObject(techs.TA.20K.df, min.cells = 5, min.genes = 1, normalization.method = "LogNormalize", scale.factor = 10000)#, total.expr = 1e4, project = "Dicroce")
TA.seurat <- ScaleData(object= TA.seurat)
TA.seurat <- FindVariableGenes(object = TA.seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.2, x.high.cutoff = 6, y.cutoff = 0.5)
print("TA HVG:")
print(dim(TA.seurat@scale.data))
print(length(TA.seurat@var.genes))
TA.seurat <- RunPCA(TA.seurat, pc.genes = TA.seurat@var.genes, do.print = F)
TA.colors.nUMI = TA.seurat@meta.data[names(TA.seurat@ident), "nUMI"]
names(TA.colors.nUMI) <- names(TA.seurat@ident)
TA.data.plot <- as.data.frame(TA.seurat@dr$pca@cell.embeddings[,1:8])
TA.data.plot <- cbind(TA.data.plot, nUMI=TA.colors.nUMI, techs=TA.seurat@meta.data$orig.ident)

p0 <- PCAPlot(TA.seurat, 1,2, group.by = "orig.ident")
p00 <- ggplot(data = TA.data.plot, mapping = aes(x = PC1, y = PC2, color=TA.colors.nUMI)) + 
  geom_point(size = 1,shape = 16)+ scale_color_gradient( low = "grey", high = "red")

p1 <- ggplot(data = TA.data.plot, mapping = aes(x = PC1, y = PC2, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p11 <- ggMarginal(p1, data = TA.data.plot, x= TA.colors.nUMI, y= TA.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

p2 <- ggplot(data = TA.data.plot, mapping = aes(x = PC3, y = PC4, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p22 <- ggMarginal(p2, data = TA.data.plot, x= TA.colors.nUMI, y= TA.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

p3 <- ggplot(data = TA.data.plot, mapping = aes(x = PC5, y = PC6, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p33 <- ggMarginal(p3, data = TA.data.plot, x= TA.colors.nUMI, y= TA.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

p4 <- ggplot(data = TA.data.plot, mapping = aes(x = PC7, y = PC8, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p44 <- ggMarginal(p4, data = TA.data.plot, x= TA.colors.nUMI, y= TA.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")


TA.plot.list <- list(p0,p00,p1,p11,p2,p22,p3,p33,p4,p44)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCAseurat_TA_V6.pdf")
invisible(lapply(TA.plot.list, print))
dev.off()
save(TA.data.plot, file="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCAseurat_TA_dataPlot.RData")
save(TA.seurat, file="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCAseurat_TA_SeuratObj.RData")

#TA Cluster tree
TA.seurat@ident <- TA.seurat@meta.data$orig.ident
names(TA.seurat@ident) <- rownames(TA.seurat@meta.data)
TA.seurat <- BuildClusterTree(TA.seurat, genes.use = TA.seurat@var.genes, pcs.use = 1:10)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_TA_ClusterTree_pc1to10.pdf")
plotClusterTree(TA.seurat)
dev.off()

print("PCA DS plots for TA done")


#all plottings

pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS_Secretory_plots_V3.pdf")
ggplot(DSth.df.Secretory, aes(x=DSthNum, y=nGenes, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))
ggplot(DSth.df.Secretory, aes(x=DSthNum, y=nUMIs, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))

#ggplot(DSth.df.Secretory, aes(x=DSthNum, y=nGenes, fill = DSth)) +geom_boxplot() + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(group = 1))

ggplot(data=DSth.df.Secretory, aes(x=DStech, y=nGenes, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
ggplot(data=DSth.df.Secretory, aes(x=DStech, y=nUMIs, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
dev.off()



pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS_TA_plots_V3.pdf")
ggplot(DSth.df.TA, aes(x=DSthNum, y=nGenes, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))
ggplot(DSth.df.TA, aes(x=DSthNum, y=nUMIs, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))

#ggplot(DSth.df.TA, aes(x=DSthNum, y=nGenes, fill = DSth)) +geom_boxplot() + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(group = 1))

ggplot(data=DSth.df.TA, aes(x=DStech, y=nGenes, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
ggplot(data=DSth.df.TA, aes(x=DStech, y=nUMIs, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
dev.off()

