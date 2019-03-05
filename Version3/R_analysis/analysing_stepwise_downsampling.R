library(scater)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(biomaRt)
library(tibble)
library(Seurat)
library(data.table)

# loading Annotated seurat objects ####
load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Human/celseq_seu.obj.RData")
CELseq2.hsap.obj <- celseq
rm(celseq)
CELseq2.hsap.metadata <- CELseq2.hsap.obj@meta.data
#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Mouse/celseq_seu.obj.RData")
#CELseq2.mmus.obj <- celseq
#rm(celseq)
#CELseq2.mmus.metadata <- CELseq2.mmus.obj@meta.data


load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Human/marsseq_seu.obj.RData")
MARSseq.hsap.obj <- marsseq
rm(marsseq)
MARSseq.hsap.metadata <- MARSseq.hsap.obj@meta.data
#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Mouse/marsseq_seu.obj.RData")
#MARSseq.mmus.obj <- marsseq
#rm(marsseq)
#MARSseq.mmus.metadata <- MARSseq.mmus.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Human/quartzseq_seu.obj.RData")
QUARTZseq.hsap.obj <- quartzseq
rm(quartzseq)
QUARTZseq.hsap.metadata <- QUARTZseq.hsap.obj@meta.data
#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Mouse/quartzseq_seu.obj.RData")
#QUARTZseq.mmus.obj <- quartzseq
#rm(quartzseq)
#QUARTZseq.mmus.metadata <- QUARTZseq.mmus.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Human/scrbseq_seu.obj.RData")
SCRBseq.hsap.obj <- scrbseq
rm(scrbseq)
SCRBseq.hsap.metadata <- SCRBseq.hsap.obj@meta.data
#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Mouse/scrbseq_seu.obj.RData")
#SCRBseq.mmus.obj <- scrbseq
#rm(scrbseq)
#SCRBseq.mmus.metadata <- SCRBseq.mmus.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Human/smartseq_seu.obj.RData")
SMARTseq2.hsap.obj <- smartseq
rm(smartseq)
SMARTseq2.hsap.metadata <- SMARTseq2.hsap.obj@meta.data
#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Mouse/smartseq_seu.obj.RData")
#SMARTseq2.mmus.obj <- smartseq
#rm(smartseq)
#SMARTseq2.mmus.metadata <- SMARTseq2.mmus.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Human/c1ht.s_seu.obj.RData")
C1HTsmall.hsap.obj <- c1ht.s
rm(c1ht.s)
C1HTsmall.hsap.metadata <- C1HTsmall.hsap.obj@meta.data
#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Mouse/c1ht.s_seu.obj.RData")
#C1HTsmall.mmus.obj <- c1ht.s
#rm(c1ht.s)
#C1HTsmall.mmus.metadata <- C1HTsmall.mmus.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Human/c1ht.m_seu.obj.RData")
C1HTmedium.hsap.obj <- c1ht.m
rm(c1ht.m)
C1HTmedium.hsap.metadata <- C1HTmedium.hsap.obj@meta.data
#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Mouse/c1ht.m_seu.obj.RData")
#C1HTmedium.mmus.obj <- c1ht.m
#rm(c1ht.m)
#C1HTmedium.mmus.metadata <- C1HTmedium.mmus.obj@meta.data


load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Human/chromium_seu.obj.RData")
Chromium.hsap.obj <- chromium
rm(chromium)
Chromium.hsap.metadata <- Chromium.hsap.obj@meta.data
#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Mouse/chromium_seu.obj.RData")
#Chromium.mmus.obj <- chromium
#rm(chromium)
#Chromium.mmus.metadata <- Chromium.mmus.obj@meta.data


load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Human/nuclei_seu.obj.RData")
ChromiumNuclei.hsap.obj <- nuclei
rm(nuclei)
ChromiumNuclei.hsap.metadata <- ChromiumNuclei.hsap.obj@meta.data
#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Mouse/nuclei_seu.obj.RData")
#ChromiumNuclei.mmus.obj <- nuclei
#rm(nuclei)
#ChromiumNuclei.mmus.metadata <- ChromiumNuclei.mmus.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Human/ddseq_seu.obj.RData")
ddSEQ.hsap.obj <- ddseq
rm(ddseq)
ddSEQ.hsap.metadata <- ddSEQ.hsap.obj@meta.data
#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Mouse/ddseq_seu.obj.RData")
#ddSEQ.mmus.obj <- ddseq
#rm(ddseq)
#ddSEQ.mmus.metadata <- ddSEQ.mmus.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Human/dropseq_seu.obj.RData")
Dropseq.hsap.obj <- dropseq
rm(dropseq)
Dropseq.hsap.metadata <- Dropseq.hsap.obj@meta.data
#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Mouse/dropseq_seu.obj.RData")
#Dropseq.mmus.obj <- dropseq
#rm(dropseq)
#Dropseq.mmus.metadata <- Dropseq.mmus.obj@meta.data


load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Human/icell8_seu.obj.RData")
ICELL8.hsap.obj <- icell8
rm(icell8)
ICELL8.hsap.metadata <- ICELL8.hsap.obj@meta.data
#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Mouse/icell8_seu.obj.RData")
#ICELL8.mmus.obj <- icell8
#rm(icell8)
#ICELL8.mmus.metadata <- ICELL8.mmus.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Human/onecb_seu.obj.RData")
inDrop.hsap.obj <- onecb
rm(onecb)
inDrop.hsap.metadata <- inDrop.hsap.obj@meta.data
#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Mouse/indrop_seu.obj.RData")
#inDrop.mmus.obj <- indrop
#rm(indrop)
#inDrop.mmus.metadata <- inDrop.mmus.obj@meta.data

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Human/seqwellV2_seu.obj.RData")
SeqwellV2.hsap.obj <- seqwell
rm(seqwell)
SeqwellV2.hsap.metadata <- SeqwellV2.hsap.obj@meta.data
#load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects_final/Mouse/seqwellV2_seu.obj.RData")
#SeqwellV2.mmus.obj <- seqwell
#rm(seqwell)
#SeqwellV2.mmus.metadata <- SeqwellV2.mmus.obj@meta.data

# End ####

# Taking out HEK cells of each technique ####
CELseq2.HEK <- rownames(CELseq2.hsap.metadata)[which(CELseq2.hsap.metadata$clean.id == "HEK")]
MARSseq.HEK <- rownames(MARSseq.hsap.metadata)[which(MARSseq.hsap.metadata$clean.id == "HEK")]
QUARTZseq.HEK <- rownames(QUARTZseq.hsap.metadata)[which(QUARTZseq.hsap.metadata$clean.id == "HEK")]
SCRBseq.HEK <- rownames(SCRBseq.hsap.metadata)[which(SCRBseq.hsap.metadata$clean.id == "HEK")]
SMARTseq2.HEK <- rownames(SMARTseq2.hsap.metadata)[which(SMARTseq2.hsap.metadata$clean.id == "HEK")]
C1HTsmall.HEK <- rownames(C1HTsmall.hsap.metadata[which(C1HTsmall.hsap.metadata$clean.id == "HEK"),])
C1HTmedium.HEK <- rownames(C1HTmedium.hsap.metadata[which(C1HTmedium.hsap.metadata$clean.id == "HEK"),])
Chromium.HEK <- rownames(Chromium.hsap.metadata[which(Chromium.hsap.metadata$clean.id == "HEK"),])
ChromiumNuclei.HEK <- rownames(ChromiumNuclei.hsap.metadata[which(ChromiumNuclei.hsap.metadata$clean.id == "HEK"),])
ddSEQ.HEK <- rownames(ddSEQ.hsap.metadata[which(ddSEQ.hsap.metadata$clean.id == "HEK"),])
Dropseq.HEK <- rownames(Dropseq.hsap.metadata[which(Dropseq.hsap.metadata$clean.id == "HEK"),])
ICELL8.HEK <- rownames(ICELL8.hsap.metadata[which(ICELL8.hsap.metadata$clean.id == "HEK"),])
inDrop.HEK <- rownames(inDrop.hsap.metadata[which(inDrop.hsap.metadata$clean.id == "HEK"),])
SeqwellV2.HEK <- rownames(SeqwellV2.hsap.metadata[which(SeqwellV2.hsap.metadata$clean.id == "HEK"),])

# End ####

# Taking out Monocytes cells of each technique ####
CELseq2.Monocytes <- rownames(CELseq2.hsap.metadata)[which(CELseq2.hsap.metadata$clean.id == "Monocytes")]
MARSseq.Monocytes <- rownames(MARSseq.hsap.metadata)[which(MARSseq.hsap.metadata$clean.id == "Monocytes")]
QUARTZseq.Monocytes <- rownames(QUARTZseq.hsap.metadata)[which(QUARTZseq.hsap.metadata$clean.id == "Monocytes")]
SCRBseq.Monocytes <- rownames(SCRBseq.hsap.metadata)[which(SCRBseq.hsap.metadata$clean.id == "Monocytes")]
SMARTseq2.Monocytes <- rownames(SMARTseq2.hsap.metadata)[which(SMARTseq2.hsap.metadata$clean.id == "Monocytes")]
C1HTsmall.Monocytes <- rownames(C1HTsmall.hsap.metadata[which(C1HTsmall.hsap.metadata$clean.id == "Monocytes"),])
C1HTmedium.Monocytes <- rownames(C1HTmedium.hsap.metadata[which(C1HTmedium.hsap.metadata$clean.id == "Monocytes"),])
Chromium.Monocytes <- rownames(Chromium.hsap.metadata[which(Chromium.hsap.metadata$clean.id == "Monocytes"),])
ChromiumNuclei.Monocytes <- rownames(ChromiumNuclei.hsap.metadata[which(ChromiumNuclei.hsap.metadata$clean.id == "Monocytes"),])
ddSEQ.Monocytes <- rownames(ddSEQ.hsap.metadata[which(ddSEQ.hsap.metadata$clean.id == "Monocytes"),])
Dropseq.Monocytes <- rownames(Dropseq.hsap.metadata[which(Dropseq.hsap.metadata$clean.id == "Monocytes"),])
ICELL8.Monocytes <- rownames(ICELL8.hsap.metadata[which(ICELL8.hsap.metadata$clean.id == "Monocytes"),])
inDrop.Monocytes <- rownames(inDrop.hsap.metadata[which(inDrop.hsap.metadata$clean.id == "Monocytes"),])
SeqwellV2.Monocytes <- rownames(SeqwellV2.hsap.metadata[which(SeqwellV2.hsap.metadata$clean.id == "Monocytes"),])

# End ####

# Taking out Bcells cells of each technique ####
CELseq2.Bcells <- rownames(CELseq2.hsap.metadata)[which(CELseq2.hsap.metadata$clean.id == "B")]
MARSseq.Bcells <- rownames(MARSseq.hsap.metadata)[which(MARSseq.hsap.metadata$clean.id == "B")]
QUARTZseq.Bcells <- rownames(QUARTZseq.hsap.metadata)[which(QUARTZseq.hsap.metadata$clean.id == "B")]
SCRBseq.Bcells <- rownames(SCRBseq.hsap.metadata)[which(SCRBseq.hsap.metadata$clean.id == "B")]
SMARTseq2.Bcells <- rownames(SMARTseq2.hsap.metadata)[which(SMARTseq2.hsap.metadata$clean.id == "B")]
C1HTsmall.Bcells <- rownames(C1HTsmall.hsap.metadata[which(C1HTsmall.hsap.metadata$clean.id == "B"),])
C1HTmedium.Bcells <- rownames(C1HTmedium.hsap.metadata[which(C1HTmedium.hsap.metadata$clean.id == "B"),])
Chromium.Bcells <- rownames(Chromium.hsap.metadata[which(Chromium.hsap.metadata$clean.id == "B"),])
ChromiumNuclei.Bcells <- rownames(ChromiumNuclei.hsap.metadata[which(ChromiumNuclei.hsap.metadata$clean.id == "B"),])
ddSEQ.Bcells <- rownames(ddSEQ.hsap.metadata[which(ddSEQ.hsap.metadata$clean.id == "B"),])
Dropseq.Bcells <- rownames(Dropseq.hsap.metadata[which(Dropseq.hsap.metadata$clean.id == "B"),])
ICELL8.Bcells <- rownames(ICELL8.hsap.metadata[which(ICELL8.hsap.metadata$clean.id == "B"),])
inDrop.Bcells <- rownames(inDrop.hsap.metadata[which(inDrop.hsap.metadata$clean.id == "B"),])
SeqwellV2.Bcells <- rownames(SeqwellV2.hsap.metadata[which(SeqwellV2.hsap.metadata$clean.id == "B"),])


# End ####

# Loding Downsampleded data ####
load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/CELseq2.hsap.full.SCE.jointDSmat.Robj")
CELseq2.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
CELseq2.DS.UMI <- CELseq2.DS$UMI
CELseq2.DS.Reads <- CELseq2.DS$Reads
rm(CELseq2.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/MARSseq.hsap.full.SCE.jointDSmat.Robj")
MARSseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
MARSseq.DS.UMI <- MARSseq.DS$UMI
MARSseq.DS.Reads <- MARSseq.DS$Reads
rm(MARSseq.DS)


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/QUARTZseq.hsap.full.SCE.jointDSmat.Robj")
QUARTZseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
QUARTZseq.DS.UMI <- QUARTZseq.DS$UMI
QUARTZseq.DS.Reads <- QUARTZseq.DS$Reads
rm(QUARTZseq.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/SCRBseq.hsap.full.SCE.jointDSmat.Robj")
SCRBseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
SCRBseq.DS.UMI <- SCRBseq.DS$UMI
SCRBseq.DS.Reads <- SCRBseq.DS$Reads
rm(SCRBseq.DS)

#load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/SMARTseqFINAL.hsap.full.SCE.jointDSmat.Robj")
#SMARTseq2.DS <- output.readcount.umicount.joint.mats
#rm(output.readcount.umicount.joint.mats)
#SMARTseq2.DS.UMI <- SMARTseq2.DS$UMI
#SMARTseq2.DS.Reads <- SMARTseq2.DS$Reads
#rm(SMARTseq2.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/C1HTsmall.hsap.full.SCE.jointDSmat.Robj")
C1HTsmall.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
C1HTsmall.DS.UMI <- C1HTsmall.DS$UMI
C1HTsmall.DS.Reads <- C1HTsmall.DS$Reads

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/C1HTmedium.hsap.full.SCE.jointDSmat.Robj")
C1HTmedium.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
C1HTmedium.DS.UMI <- C1HTmedium.DS$UMI
C1HTmedium.DS.Reads <- C1HTmedium.DS$Reads

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/10X2x5K.hsap.full.SCE.jointDSmat.Robj")
Chromium.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
Chromium.DS.UMI <- Chromium.DS$UMI
Chromium.DS.Reads <- Chromium.DS$Reads

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Nuclei10X.hsap.full.SCE.jointDSmat.Robj")
ChromiumNuclei.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ChromiumNuclei.DS.UMI <- ChromiumNuclei.DS$UMI
ChromiumNuclei.DS.Reads <- ChromiumNuclei.DS$Reads

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/ddSEQ.hsap.full.SCE.jointDSmat.Robj")
ddSEQ.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ddSEQ.DS.UMI <- ddSEQ.DS$UMI
ddSEQ.DS.Reads <- ddSEQ.DS$Reads

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Dropseq.hsap.full.SCE.jointDSmat.Robj")
Dropseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
Dropseq.DS.UMI <- Dropseq.DS$UMI
Dropseq.DS.Reads <- Dropseq.DS$Reads
rm(Dropseq.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/ICELL8.hsap.full.SCE.jointDSmat.Robj")
ICELL8.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ICELL8.DS.UMI <- ICELL8.DS$UMI
ICELL8.DS.Reads <- ICELL8.DS$Reads


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/1CB.hsap.full.SCE.jointDSmat.Robj")
inDrop.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
inDrop.DS.UMI <- inDrop.DS$UMI
inDrop.DS.Reads <- inDrop.DS$Reads

# End ####


#Plotting stepwise Downsampling for HEK ####

techniques <- c("CELseq2","MARSseq","QUARTZseq","SCRBseq", "C1HTsmall","C1HTmedium","Chromium", "ChromiumNuclei","ddSEQ","Dropseq","ICELL8","inDrop") #"SMARTseq2",
DSth.df.HEK <- data.frame()
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
    DSth.df.HEK <- rbind(DSth.df.HEK, DS.df)
    print(dim(DSth.df.HEK))
  }
}

save(DSth.df.HEK, file="/project/devel/alafzi/SC_Protocols/Version3/stepwise_analysis_final/Stepwise_DS_plotdf_HEK.RData")
print("Stepwise DS plots for HEK done")

#Plotting stepwise Downsampling for Monocytes ####

DSth.df.Monocytes <- data.frame()
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
    DSth.df.Monocytes <- rbind(DSth.df.Monocytes, DS.df)
    print(dim(DSth.df.Monocytes))
  }
}

save(DSth.df.Monocytes, file="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/Stepwise_DS_plotdf_Monocytes.RData")

print("Stepwise DS plots for Monocytes done")

#Plotting stepwise Downsampling for Bcell ####

DSth.df.Bcells <- data.frame()
techs.Bcells.20K.list <- list()
for (tech in techniques){
  print(tech)
  tech.DS.UMI <- get(paste(tech,".DS.UMI", sep = ""))
  tech.Bcells <- get(paste(tech,".Bcells", sep = ""))
  for (DSth in names(tech.DS.UMI)){
    print(paste(tech, DSth, sep = "_"))
    colnames(tech.DS.UMI[[DSth]]) <- gsub(x = colnames(tech.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    comm.cells <- intersect(tech.Bcells, colnames(tech.DS.UMI[[DSth]]))
    DS.mat.BcellsS <- tech.DS.UMI[[DSth]][, comm.cells]
    if (DSth == "downsampled_20000"){
      DS.mat.BcellsS.dup <- DS.mat.BcellsS
      DS.mat.BcellsS.dup$gene_id <- rownames(DS.mat.BcellsS.dup)
      techs.Bcells.20K.list[[tech]] <- as.data.frame((DS.mat.BcellsS.dup))}
    DS.gene.distribution.Bcells <- colSums(DS.mat.BcellsS[,]>0)
    DS.UMI.distribution.Bcells <- colSums(DS.mat.BcellsS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.Bcells))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.Bcells))
    DS.labels <- rep(DSth, length(DS.gene.distribution.Bcells))
    DS.techs <- rep(tech, length(DS.gene.distribution.Bcells))
    DS.df <- data.frame(nGenes = DS.gene.distribution.Bcells, nUMIs = DS.UMI.distribution.Bcells, DSthNum = DSth.number.vec, DStech = DS.techs)
    DSth.df.Bcells <- rbind(DSth.df.Bcells, DS.df)
    print(dim(DSth.df.Bcells))
  }
}

save(DSth.df.Bcells, file="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/Stepwise_DS_plotdf_Bcells.RData")

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
#pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCA_HEK.pdf")
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
HEK.colors.nUMI = HEK.seurat@meta.data[names(HEK.seurat@ident), "nUMI"]
names(HEK.colors.nUMI) <- names(HEK.seurat@ident)
HEK.data.plot <- as.data.frame(HEK.seurat@dr$pca@cell.embeddings[,1:8])
HEK.data.plot <- cbind(HEK.data.plot, nUMI=HEK.colors.nUMI, techs=HEK.seurat@meta.data$orig.ident)

p0 <- PCAPlot(HEK.seurat, 1,2, group.by = "orig.ident")
p00 <- ggplot(data = HEK.data.plot, mapping = aes(x = PC1, y = PC2, color=HEK.colors.nUMI)) + 
  geom_point(size = 1,shape = 16)+ scale_color_gradient( low = "grey", high = "red")

p1 <- ggplot(data = HEK.data.plot, mapping = aes(x = PC1, y = PC2, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p11 <- ggMarginal(p1, data = HEK.data.plot, x= HEK.colors.nUMI, y= HEK.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

p2 <- ggplot(data = HEK.data.plot, mapping = aes(x = PC3, y = PC4, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p22 <- ggMarginal(p2, data = HEK.data.plot, x= HEK.colors.nUMI, y= HEK.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

p3 <- ggplot(data = HEK.data.plot, mapping = aes(x = PC5, y = PC6, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p33 <- ggMarginal(p3, data = HEK.data.plot, x= HEK.colors.nUMI, y= HEK.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

p4 <- ggplot(data = HEK.data.plot, mapping = aes(x = PC7, y = PC8, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p44 <- ggMarginal(p4, data = HEK.data.plot, x= HEK.colors.nUMI, y= HEK.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

HEK.plot.list <- list(p0,p00,p1,p11,p2,p22,p3,p33,p4,p44)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCAseurat_HEK_V6.pdf")
invisible(lapply(HEK.plot.list, print))
dev.off()
save(HEK.data.plot, file="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCAseurat_HEK_dataPlot.RData")
save(HEK.seurat, file="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCAseurat_HEK_SeuratObj.RData")

#HEK Cluster tree
HEK.seurat@ident <- HEK.seurat@meta.data$orig.ident
names(HEK.seurat@ident) <- rownames(HEK.seurat@meta.data)
HEK.seurat <- BuildClusterTree(HEK.seurat, genes.use = HEK.seurat@var.genes, pcs.use = 1:10)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_HEK_ClusterTree_pc1to10.pdf")
plotClusterTree(HEK.seurat)
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
#pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCA_Monocytes.pdf")
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
Monocytes.colors.nUMI = Monocytes.seurat@meta.data[names(Monocytes.seurat@ident), "nUMI"]
names(Monocytes.colors.nUMI) <- names(Monocytes.seurat@ident)
Monocytes.data.plot <- as.data.frame(Monocytes.seurat@dr$pca@cell.embeddings[,1:8])
Monocytes.data.plot <- cbind(Monocytes.data.plot, nUMI=Monocytes.colors.nUMI, techs=Monocytes.seurat@meta.data$orig.ident)

p0 <- PCAPlot(Monocytes.seurat, 1,2, group.by = "orig.ident")
p00 <- ggplot(data = Monocytes.data.plot, mapping = aes(x = PC1, y = PC2, color=Monocytes.colors.nUMI)) + 
  geom_point(size = 1,shape = 16)+ scale_color_gradient( low = "grey", high = "red")

p1 <- ggplot(data = Monocytes.data.plot, mapping = aes(x = PC1, y = PC2, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p11 <- ggMarginal(p1, data = Monocytes.data.plot, x= Monocytes.colors.nUMI, y= Monocytes.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

p2 <- ggplot(data = Monocytes.data.plot, mapping = aes(x = PC3, y = PC4, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p22 <- ggMarginal(p2, data = Monocytes.data.plot, x= Monocytes.colors.nUMI, y= Monocytes.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

p3 <- ggplot(data = Monocytes.data.plot, mapping = aes(x = PC5, y = PC6, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p33 <- ggMarginal(p3, data = Monocytes.data.plot, x= Monocytes.colors.nUMI, y= Monocytes.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

p4 <- ggplot(data = Monocytes.data.plot, mapping = aes(x = PC7, y = PC8, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p44 <- ggMarginal(p4, data = Monocytes.data.plot, x= Monocytes.colors.nUMI, y= Monocytes.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")


Monocytes.plot.list <- list(p0,p00,p1,p11,p2,p22,p3,p33,p4,p44)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCAseurat_Monocytes_V6.pdf")
invisible(lapply(Monocytes.plot.list, print))
dev.off()
save(Monocytes.data.plot, file="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCAseurat_Monocytes_dataPlot.RData")
save(Monocytes.seurat, file="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCAseurat_Monocytes_SeuratObj.RData")

#Monocytes Cluster tree
Monocytes.seurat@ident <- Monocytes.seurat@meta.data$orig.ident
names(Monocytes.seurat@ident) <- rownames(Monocytes.seurat@meta.data)
Monocytes.seurat <- BuildClusterTree(Monocytes.seurat, genes.use = Monocytes.seurat@var.genes, pcs.use = 1:10)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_Monocytes_ClusterTree_pc1to10.pdf")
plotClusterTree(Monocytes.seurat)
dev.off()

print("PCA DS plots for Monocytes done")


# The PCA for celltypes Bcells
print("PCA analysis started")
library(plyr)
techs.Bcells.20K.df=join_all(techs.Bcells.20K.list, by = "gene_id", type = 'full')
rownames(techs.Bcells.20K.df) <- techs.Bcells.20K.df[,"gene_id"]
gene.id.col <- which(colnames(techs.Bcells.20K.df) == "gene_id")
techs.Bcells.20K.df <- techs.Bcells.20K.df[,-gene.id.col]
dim(techs.Bcells.20K.df)
techs.Bcells.20K.df <- na.omit(techs.Bcells.20K.df)
#techs.Bcells.20K.log.df <- log(techs.Bcells.20K.df + 1)
#pca.res <- prcomp(techs.Bcells.20K.log.df, scale. = T)
#pca.df <- as.data.frame(pca.res$rotation[,c("PC1","PC2")])
#tech.vec <- sapply(strsplit(rownames(pca.df), split = "_", fixed = T), `[`,1)
#pca.df <- cbind(pca.df,tech.vec)
#pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCA_Bcells.pdf")
#ggplot(pca.df, aes(x= PC1, y = PC2, group= tech.vec)) + geom_point(aes(color = tech.vec))
#dev.off()
Bcells.seurat <- CreateSeuratObject(techs.Bcells.20K.df, min.cells = 5, min.genes = 1, normalization.method = "LogNormalize", scale.factor = 10000)#, total.expr = 1e4, project = "Dicroce")
Bcells.seurat <- ScaleData(object= Bcells.seurat)
Bcells.seurat <- FindVariableGenes(object = Bcells.seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.2, x.high.cutoff = 6, y.cutoff = 0.5)
print("Bcells HVG:")
print(dim(Bcells.seurat@scale.data))
print(length(Bcells.seurat@var.genes))
Bcells.seurat <- RunPCA(Bcells.seurat, pc.genes = Bcells.seurat@var.genes, do.print = F)
Bcells.colors.nUMI = Bcells.seurat@meta.data[names(Bcells.seurat@ident), "nUMI"]
names(Bcells.colors.nUMI) <- names(Bcells.seurat@ident)
Bcells.data.plot <- as.data.frame(Bcells.seurat@dr$pca@cell.embeddings[,1:8])
Bcells.data.plot <- cbind(Bcells.data.plot, nUMI=Bcells.colors.nUMI, techs=Bcells.seurat@meta.data$orig.ident)

p0 <- PCAPlot(Bcells.seurat, 1,2, group.by = "orig.ident")
p00 <- ggplot(data = Bcells.data.plot, mapping = aes(x = PC1, y = PC2, color=Bcells.colors.nUMI)) + 
  geom_point(size = 1,shape = 16)+ scale_color_gradient( low = "grey", high = "red")

p1 <- ggplot(data = Bcells.data.plot, mapping = aes(x = PC1, y = PC2, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p11 <- ggMarginal(p1, data = Bcells.data.plot, x= Bcells.colors.nUMI, y= Bcells.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

p2 <- ggplot(data = Bcells.data.plot, mapping = aes(x = PC3, y = PC4, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p22 <-ggMarginal(p2, data = Bcells.data.plot, x= Bcells.colors.nUMI, y= Bcells.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

p3 <- ggplot(data = Bcells.data.plot, mapping = aes(x = PC5, y = PC6, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p33 <- ggMarginal(p3, data = Bcells.data.plot, x= Bcells.colors.nUMI, y= Bcells.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

p4 <- ggplot(data = Bcells.data.plot, mapping = aes(x = PC7, y = PC8, color=techs)) + 
  geom_point(size = 1,shape = 16)+ theme(legend.position="none")
p44 <- ggMarginal(p4, data = Bcells.data.plot, x= Bcells.colors.nUMI, y= Bcells.colors.nUMI, type = "density", margins = "both", size = 4, color = "pink", fill = "pink")

Bcells.plot.list <- list(p0,p00,p1,p11,p2,p22,p3,p33,p4,p44)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCAseurat_Bcells_V6.pdf")
invisible(lapply(Bcells.plot.list, print))
dev.off()
save(Bcells.data.plot, file="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCAseurat_Bcells_dataPlot.RData")
save(Bcells.seurat, file="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_PCAseurat_Bcells_SeuratObj.RData")

#Bcells Cluster tree
Bcells.seurat@ident <- Bcells.seurat@meta.data$orig.ident
names(Bcells.seurat@ident) <- rownames(Bcells.seurat@meta.data)
Bcells.seurat <- BuildClusterTree(Bcells.seurat, genes.use = Bcells.seurat@var.genes, pcs.use = 1:10)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS20K_Bcells_ClusterTree_pc1to10.pdf")
plotClusterTree(Bcells.seurat)
dev.off()

print("PCA DS plots for Bcells done")


#all plottings

pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS_HEK_plots_V3.pdf")
ggplot(DSth.df.HEK, aes(x=DSthNum, y=nGenes, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))
ggplot(DSth.df.HEK, aes(x=DSthNum, y=nUMIs, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))

#ggplot(DSth.df.HEK, aes(x=DSthNum, y=nGenes, fill = DSth)) +geom_boxplot() + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(group = 1))

ggplot(data=DSth.df.HEK, aes(x=DStech, y=nGenes, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
ggplot(data=DSth.df.HEK, aes(x=DStech, y=nUMIs, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
dev.off()



pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS_Monocytes_plots_V3.pdf")
ggplot(DSth.df.Monocytes, aes(x=DSthNum, y=nGenes, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))
ggplot(DSth.df.Monocytes, aes(x=DSthNum, y=nUMIs, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))

#ggplot(DSth.df.Monocytes, aes(x=DSthNum, y=nGenes, fill = DSth)) +geom_boxplot() + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(group = 1))

ggplot(data=DSth.df.Monocytes, aes(x=DStech, y=nGenes, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
ggplot(data=DSth.df.Monocytes, aes(x=DStech, y=nUMIs, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
dev.off()


pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwise_analysis_final/all_techs_stepwise_DS_Bcells_plots_V3.pdf")
ggplot(DSth.df.Bcells, aes(x=DSthNum, y=nGenes, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))
ggplot(DSth.df.Bcells, aes(x=DSthNum, y=nUMIs, group =DStech))  + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(color=DStech))

#ggplot(DSth.df.Bcells, aes(x=DSthNum, y=nGenes, fill = DSth)) +geom_boxplot() + geom_smooth(method = "lm", formula = y ~ log(x), se = T, aes(group = 1))

ggplot(data=DSth.df.Bcells, aes(x=DStech, y=nGenes, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
ggplot(data=DSth.df.Bcells, aes(x=DStech, y=nUMIs, fill=DStech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +facet_wrap(. ~ DSthNum, scales = "free") 
dev.off()
