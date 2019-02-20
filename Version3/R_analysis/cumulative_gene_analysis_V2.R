library(scater)
library(ggplot2)
library(ggExtra)
library(gridExtra)
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

load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/CB1_data_seu.obj_res0.5_dim7.RData")
CB1.obj <- data
rm(data)
CB1.metadata <- CB1.obj@meta.data

# End ####

# Taking out HEK cells of each technique ####
MARSseq.HEK <- names(MARSseq.obj@ident)[which(MARSseq.obj@ident == "HEK cells")]
if(length(MARSseq.HEK) > 50){ MARSseq.HEK <- sample(MARSseq.HEK,50)}
CELseq2.HEK <- names(CELseq2.obj@ident)[which(CELseq2.obj@ident == "HEK cells")]
if(length(CELseq2.HEK) > 50){ CELseq2.HEK <- sample(CELseq2.HEK,50)}
QUARTZseq.HEK <- names(QUARTZseq.obj@ident)[which(QUARTZseq.obj@ident == "HEK cells")]
if(length(QUARTZseq.HEK) > 50){ QUARTZseq.HEK <- sample(QUARTZseq.HEK,50)}
Dropseq.HEK <- names(Dropseq.obj@ident)[which(Dropseq.obj@ident == "HEK cells")]
if(length(Dropseq.HEK) > 50){ Dropseq.HEK <- sample(Dropseq.HEK,50)}
SCRBseq.HEK <- names(SCRBseq.obj@ident)[which(SCRBseq.obj@ident == "HEK cells")]
if(length(SCRBseq.HEK) > 50){ SCRBseq.HEK <- sample(SCRBseq.HEK,50)}
X10Scilife.HEK <- names(X10Scilife.obj@ident)[which(X10Scilife.obj@ident == "HEK cells")]
if(length(X10Scilife.HEK) > 50){ X10Scilife.HEK <- sample(X10Scilife.HEK,50)}
X10Nuclei.HEK <- names(X10Nuclei.obj@ident)[which(X10Nuclei.obj@ident == "HEK cells")]
if(length(X10Nuclei.HEK) > 50){ X10Nuclei.HEK <- sample(X10Nuclei.HEK,50)}
ICELL8.HEK <- names(ICELL8.obj@ident)[which(ICELL8.obj@ident == "HEK cells")]
if(length(ICELL8.HEK) > 50){ ICELL8.HEK <- sample(ICELL8.HEK,50)}
ddSEQ.HEK <- names(ddSEQ.obj@ident)[which(ddSEQ.obj@ident == "HEK cells1")] # NOTICE THIS **********
if(length(ddSEQ.HEK) > 50){ ddSEQ.HEK <- sample(ddSEQ.HEK,50)}
ddSEQexp1.HEK <- names(ddSEQexp1.obj@ident)[which(ddSEQexp1.obj@ident == "HEK cells")]
if(length(ddSEQexp1.HEK) > 50){ ddSEQexp1.HEK <- sample(ddSEQexp1.HEK,50)}
CB1.HEK <- names(CB1.obj@ident)[which(CB1.obj@ident == "HEK cells")]
if(length(CB1.HEK) > 50){ CB1.HEK <- sample(CB1.HEK,50)}

# End ####

# Taking out Monocytes cells of each technique ####
MARSseq.Monocytes <- names(MARSseq.obj@ident)[which(MARSseq.obj@ident == "CD14+ Monocytes")]
if(length(MARSseq.Monocytes) > 50){ MARSseq.Monocytes <- sample(MARSseq.Monocytes,50)}
CELseq2.Monocytes <- names(CELseq2.obj@ident)[which(CELseq2.obj@ident == "CD14+ Monocytes")]
if(length(CELseq2.Monocytes) > 50){ CELseq2.Monocytes <- sample(CELseq2.Monocytes,50)}
QUARTZseq.Monocytes <- names(QUARTZseq.obj@ident)[which(QUARTZseq.obj@ident == "CD14+ Monocytes")]
if(length(QUARTZseq.Monocytes) > 50){ QUARTZseq.Monocytes <- sample(QUARTZseq.Monocytes,50)}
Dropseq.Monocytes <- names(Dropseq.obj@ident)[which(Dropseq.obj@ident == "CD14+ Monocytes")]
if(length(Dropseq.Monocytes) > 50){ Dropseq.Monocytes <- sample(Dropseq.Monocytes,50)}
SCRBseq.Monocytes <- names(SCRBseq.obj@ident)[which(SCRBseq.obj@ident == "CD14+ Monocytes")]
if(length(SCRBseq.Monocytes) > 50){ SCRBseq.Monocytes <- sample(SCRBseq.Monocytes,50)}
X10Scilife.Monocytes <- names(X10Scilife.obj@ident)[which(X10Scilife.obj@ident == "CD14+ and FCGR3A+ Monocytes")]
if(length(X10Scilife.Monocytes) > 50){ X10Scilife.Monocytes <- sample(X10Scilife.Monocytes,50)}
X10Nuclei.Monocytes <- names(X10Nuclei.obj@ident)[which(X10Nuclei.obj@ident == "CD14+ Monocytes")]
if(length(X10Nuclei.Monocytes) > 50){ X10Nuclei.Monocytes <- sample(X10Nuclei.Monocytes,50)}
ICELL8.Monocytes <- names(ICELL8.obj@ident)[which(ICELL8.obj@ident == "CD14+ and FCGR3A+ Monocytes")]
if(length(ICELL8.Monocytes) > 50){ ICELL8.Monocytes <- sample(ICELL8.Monocytes,50)}
ddSEQ.Monocytes <- names(ddSEQ.obj@ident)[which(ddSEQ.obj@ident == "CD14+ and FCGR3A+ Monocytes")]
if(length(ddSEQ.Monocytes) > 50){ ddSEQ.Monocytes <- sample(ddSEQ.Monocytes,50)}
ddSEQexp1.Monocytes <- names(ddSEQexp1.obj@ident)[which(ddSEQexp1.obj@ident == "CD14+ and FCGR3A+ Monocytes")]
if(length(ddSEQexp1.Monocytes) > 50){ ddSEQexp1.Monocytes <- sample(ddSEQexp1.Monocytes,50)}
CB1.Monocytes <- names(CB1.obj@ident)[which(CB1.obj@ident == "CD14+ Monocytes")]
if(length(CB1.Monocytes) > 50){ CB1.Monocytes <- sample(CB1.Monocytes,50)}
# End ####

# Taking out Bcells cells of each technique ####
MARSseq.Bcells <- names(MARSseq.obj@ident)[which(MARSseq.obj@ident == "B cells")]
if(length(MARSseq.Bcells) > 50){ MARSseq.Bcells <- sample(MARSseq.Bcells,50)}
CELseq2.Bcells <- names(CELseq2.obj@ident)[which(CELseq2.obj@ident == "B cells")]
if(length(CELseq2.Bcells) > 50){ CELseq2.Bcells <- sample(CELseq2.Bcells,50)}
QUARTZseq.Bcells <- names(QUARTZseq.obj@ident)[which(QUARTZseq.obj@ident == "B cells")]
if(length(QUARTZseq.Bcells) > 50){ QUARTZseq.Bcells <- sample(QUARTZseq.Bcells,50)}
Dropseq.Bcells <- names(Dropseq.obj@ident)[which(Dropseq.obj@ident == "B cells")]
if(length(Dropseq.Bcells) > 50){ Dropseq.Bcells <- sample(Dropseq.Bcells,50)}
SCRBseq.Bcells <- names(SCRBseq.obj@ident)[which(SCRBseq.obj@ident == "B cells")]
if(length(SCRBseq.Bcells) > 50){ SCRBseq.Bcells <- sample(SCRBseq.Bcells,50)}
X10Scilife.Bcells <- names(X10Scilife.obj@ident)[which(X10Scilife.obj@ident == "B cells")]
if(length(X10Scilife.Bcells) > 50){ X10Scilife.Bcells <- sample(X10Scilife.Bcells,50)}
X10Nuclei.Bcells <- names(X10Nuclei.obj@ident)[which(X10Nuclei.obj@ident == "B cells")]
if(length(X10Nuclei.Bcells) > 50){ X10Nuclei.Bcells <- sample(X10Nuclei.Bcells,50)}
ICELL8.Bcells <- names(ICELL8.obj@ident)[which(ICELL8.obj@ident == "B cells")]
if(length(ICELL8.Bcells) > 50){ ICELL8.Bcells <- sample(ICELL8.Bcells,50)}
ddSEQ.Bcells <- names(ddSEQ.obj@ident)[which(ddSEQ.obj@ident == "B cells")]
if(length(ddSEQ.Bcells) > 50){ ddSEQ.Bcells <- sample(ddSEQ.Bcells,50)}
ddSEQexp1.Bcells <- names(ddSEQexp1.obj@ident)[which(ddSEQexp1.obj@ident == "B cells")]
if(length(ddSEQexp1.Bcells) > 50){ ddSEQexp1.Bcells <- sample(ddSEQexp1.Bcells,50)}
CB1.Bcells <- names(CB1.obj@ident)[which(CB1.obj@ident == "B cells")]
if(length(CB1.Bcells) > 50){ CB1.Bcells <- sample(CB1.Bcells,50)}

# End ####

# End ####

# Loding Downsampleded data ####
load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/MARSseq.hsap.full.SCE.jointDSmat.Robj")
MARSseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
MARSseq.DS.UMI <- MARSseq.DS$UMI
rm(MARSseq.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/CELseq2.hsap.full.SCE.jointDSmat.Robj")
CELseq2.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
CELseq2.DS.UMI <- CELseq2.DS$UMI
rm(CELseq2.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/QUARTZseq.hsap.full.SCE.jointDSmat.Robj")
QUARTZseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
QUARTZseq.DS.UMI <- QUARTZseq.DS$UMI
rm(QUARTZseq.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Dropseq.hsap.full.SCE.jointDSmat.Robj")
Dropseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
Dropseq.DS.UMI <- Dropseq.DS$UMI
rm(Dropseq.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/SCRBseq.hsap.full.SCE.jointDSmat.Robj")
SCRBseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
SCRBseq.DS.UMI <- SCRBseq.DS$UMI
rm(SCRBseq.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/10XScilife.hsap.full.SCE.jointDSmat.Robj")
X10Scilife.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
X10Scilife.DS.UMI <- X10Scilife.DS$UMI
rm(X10Scilife.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Nuclei10X.hsap.full.SCE.jointDSmat.Robj")
X10Nuclei.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
X10Nuclei.DS.UMI <- X10Nuclei.DS$UMI
rm(X10Nuclei.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/ICELL8.hsap.full.SCE.jointDSmat.Robj")
ICELL8.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ICELL8.DS.UMI <- ICELL8.DS$UMI
rm(ICELL8.DS)


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/ddSEQexp1.hsap.full.SCE.jointDSmat.Robj")
ddSEQexp1.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ddSEQexp1.DS.UMI <- ddSEQexp1.DS$UMI
rm(ddSEQexp1.DS)


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/ddSEQ.hsap.full.SCE.jointDSmat.Robj")
ddSEQ.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ddSEQ.DS.UMI <- ddSEQ.DS$UMI
rm(ddSEQ.DS)


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/C1HT.hsap.full.SCE.jointDSmat.Robj")
C1HT.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
C1HT.DS.UMI <- C1HT.DS$UMI
rm(C1HT.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/10X8x10K.hsap.full.SCE.jointDSmat.Robj")
X108x10.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
X108x10.DS.UMI <- X108x10.DS$UMI
rm(X108x10.DS)

load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/1CB.hsap.full.SCE.jointDSmat.Robj")
CB1.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
CB1.DS.UMI <- CB1.DS$UMI
rm(CB1.DS)

# End ####


#Plotting stepwise Downsampling for HEK ####

techniques <- c("MARSseq", "QUARTZseq", "CELseq2", "Dropseq", "SCRBseq", "X10Scilife", "X10Nuclei", "ICELL8", "ddSEQexp1", "CB1")#, "X108x10" , "ddSEQ", "ddSEQexp1", "C1HT"
DSth.df <- data.frame()
techs.HEK.20K.list <- list()
for (tech in techniques){
  print(tech)
  tech.DS.UMI <- get(paste(tech,".DS.UMI", sep = ""))
  tech.HEK <- get(paste(tech,".HEK", sep = ""))
  DSth <- "downsampled_20000"
  print(paste(tech, DSth, sep = "_"))
  tech.DS.UMI.20K <- tech.DS.UMI$downsampled_20000
  colnames(tech.DS.UMI.20K) <- gsub(x = colnames(tech.DS.UMI.20K), pattern = "\\.", replacement = "_")
  comm.cells <- intersect(tech.HEK, colnames(tech.DS.UMI.20K))
  DS.mat.HEKS <- tech.DS.UMI.20K[, comm.cells]
    
  #DS.mat.HEKS.dup <- DS.mat.HEKS
  #DS.mat.HEKS.dup$gene_id <- rownames(DS.mat.HEKS.dup)
  techs.HEK.20K.list[[tech]] <- as.data.frame((DS.mat.HEKS))
  rm(tech.DS.UMI)
}

#Plotting stepwise Downsampling for Monocytes ####

techniques <- c("MARSseq", "QUARTZseq", "CELseq2", "Dropseq", "SCRBseq", "X10Scilife", "X10Nuclei", "ICELL8", "ddSEQexp1", "CB1")#, "X108x10" , "ddSEQ", "ddSEQexp1", "C1HT"
DSth.df <- data.frame()
techs.Monocytes.20K.list <- list()
for (tech in techniques){
  print(tech)
  tech.DS.UMI <- get(paste(tech,".DS.UMI", sep = ""))
  tech.Monocytes <- get(paste(tech,".Monocytes", sep = ""))
  DSth <- "downsampled_20000"
  print(paste(tech, DSth, sep = "_"))
  tech.DS.UMI.20K <- tech.DS.UMI$downsampled_20000
  colnames(tech.DS.UMI.20K) <- gsub(x = colnames(tech.DS.UMI.20K), pattern = "\\.", replacement = "_")
  comm.cells <- intersect(tech.Monocytes, colnames(tech.DS.UMI.20K))
  DS.mat.Monocytes <- tech.DS.UMI.20K[, comm.cells]
    
  #DS.mat.Monocytes.dup <- DS.mat.Monocytes
  #DS.mat.Monocytes.dup$gene_id <- rownames(DS.mat.Monocytes.dup)
  techs.Monocytes.20K.list[[tech]] <- as.data.frame((DS.mat.Monocytes))
  rm(tech.DS.UMI)
}
#Plotting stepwise Downsampling for Bcell ####

techniques <- c("MARSseq", "QUARTZseq", "CELseq2", "Dropseq", "SCRBseq", "X10Scilife", "X10Nuclei", "ICELL8", "ddSEQexp1", "CB1")#, "X108x10" , "ddSEQ", "ddSEQexp1", "C1HT"
DSth.df <- data.frame()
techs.Bcells.10K.list <- list()
for (tech in techniques){
  print(tech)
  tech.DS.UMI <- get(paste(tech,".DS.UMI", sep = ""))
  tech.Bcells <- get(paste(tech,".Bcells", sep = ""))
  DSth <- "downsampled_20000"
  print(paste(tech, DSth, sep = "_"))
  tech.DS.UMI.10K <- tech.DS.UMI$downsampled_20000
  colnames(tech.DS.UMI.10K) <- gsub(x = colnames(tech.DS.UMI.10K), pattern = "\\.", replacement = "_")
  comm.cells <- intersect(tech.Bcells, colnames(tech.DS.UMI.10K))
  DS.mat.Bcells <- tech.DS.UMI.10K[, comm.cells]
    
  #DS.mat.Bcells.dup <- DS.mat.Bcells
  #DS.mat.Bcells.dup$gene_id <- rownames(DS.mat.Bcells.dup)
  techs.Bcells.10K.list[[tech]] <- as.data.frame((DS.mat.Bcells))
  rm(tech.DS.UMI)
}


# calculating the cumulatives
print("calculating the cumulatives For HEK")
HEK.plot.df <- data.frame()
for (tech in names(techs.HEK.20K.list)){
    print(tech)
    tech.data <- techs.HEK.20K.list[[tech]]
    tech.cumul.gene.numbers <- rep(NA, ncol(tech.data))
    for (cell in 1:ncol(tech.data)){
      print(cell)
      sample.size = cell
      sample.gene.numbers <- rep(NA,50)
      for (i in 1:50){
        selected.cells <- sample(colnames(tech.data), sample.size)
        sample.gene.names <- rep(NA, 70000)
        for (cell2 in selected.cells){
          vec.start.point <- (length(sample.gene.names[which(is.na(sample.gene.names) == FALSE)])) +1
          detected.genes <- rownames(tech.data[which(tech.data[,cell2] > 0),])
          new.genes <- setdiff(detected.genes,sample.gene.names)
          vec.end.point <- vec.start.point + length(new.genes) -1
          sample.gene.names[vec.start.point:vec.end.point] <- new.genes
          #detected.genes <- rownames(tech.data[which(tech.data[,cell2] > 0),])
          #sample.gene.names <- c(sample.gene.names,detected.genes)
          }
        sample.gene.names <- na.omit(sample.gene.names)
        num.uniq.genes <- length(unique(sample.gene.names))
        sample.gene.numbers[i] <- num.uniq.genes
      }
      sample.gene.numbers <- na.omit(sample.gene.numbers)
      print(paste(cell," = ", mean(sample.gene.numbers), sep = ""))
      tech.cumul.gene.numbers[cell] <- mean(sample.gene.numbers)
    }
    HEK.plot.df <- rbind(HEK.plot.df, data.frame(Cumul= tech.cumul.gene.numbers, tech= rep(tech, length(tech.cumul.gene.numbers)), cell.num= seq(1:length(tech.cumul.gene.numbers))))
}

save(HEK.plot.df, "/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene/Cumulative_gene_dist_DS20K_HEK_dataPlot.RData")

pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene/Cumulative_gene_dist_DS20K_HEK_V2.pdf")
ggplot(data=HEK.plot.df, aes(x=cell.num, y=Cumul, group=tech)) + geom_point(aes(color=tech)) + geom_smooth(aes(color=tech))
dev.off()
print("HEK cumulative done!")


for (tech in names(techs.Monocytes.20K.list)){
    print(tech)
    tech.data <- techs.Monocytes.20K.list[[tech]]
    tech.cumul.gene.numbers <- rep(NA, ncol(tech.data))
    for (cell in 1:ncol(tech.data)){
      print(cell)
      sample.size = cell
      sample.gene.numbers <- rep(NA,50)
      for (i in 1:50){
        selected.cells <- sample(colnames(tech.data), sample.size)
        sample.gene.names <- rep(NA, 70000)
        for (cell2 in selected.cells){
          vec.start.point <- (length(sample.gene.names[which(is.na(sample.gene.names) == FALSE)])) +1
          detected.genes <- rownames(tech.data[which(tech.data[,cell2] > 0),])
          new.genes <- setdiff(detected.genes,sample.gene.names)
          vec.end.point <- vec.start.point + length(new.genes) -1
          sample.gene.names[vec.start.point:vec.end.point] <- new.genes
          #detected.genes <- rownames(tech.data[which(tech.data[,cell2] > 0),])
          #sample.gene.names <- c(sample.gene.names,detected.genes)
          }
        sample.gene.names <- na.omit(sample.gene.names)
        num.uniq.genes <- length(unique(sample.gene.names))
        sample.gene.numbers[i] <- num.uniq.genes
      }
      sample.gene.numbers <- na.omit(sample.gene.numbers)
      print(paste(cell," = ", mean(sample.gene.numbers), sep = ""))
      tech.cumul.gene.numbers[cell] <- mean(sample.gene.numbers)
    }
    Monocytes.plot.df <- rbind(Monocytes.plot.df, data.frame(Cumul= tech.cumul.gene.numbers, tech= rep(tech, length(tech.cumul.gene.numbers)), cell.num= seq(1:length(tech.cumul.gene.numbers))))
}

save(Monocytes.plot.df, "/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene/Cumulative_gene_dist_DS20K_Monocytes_dataPlot.RData")

pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene/Cumulative_gene_dist_DS20K_Monocytes_V2.pdf")
ggplot(data=Monocytes.plot.df, aes(x=cell.num, y=Cumul, group=tech)) + geom_point(aes(color=tech)) + geom_smooth(aes(color=tech))
dev.off()
print("Monocytes cumulative done!")


Bcells.plot.df <- data.frame()
for (tech in names(techs.Bcells.10K.list)){
    print(tech)
    tech.data <- techs.Bcells.10K.list[[tech]]
    tech.cumul.gene.numbers <- rep(NA, ncol(tech.data))
    for (cell in 1:ncol(tech.data)){
      print(cell)
      sample.size = cell
      sample.gene.numbers <- rep(NA,50)
      for (i in 1:50){
        selected.cells <- sample(colnames(tech.data), sample.size)
        sample.gene.names <- rep(NA, 70000)
        for (cell2 in selected.cells){
          vec.start.point <- (length(sample.gene.names[which(is.na(sample.gene.names) == FALSE)])) +1
          detected.genes <- rownames(tech.data[which(tech.data[,cell2] > 0),])
          new.genes <- setdiff(detected.genes,sample.gene.names)
          vec.end.point <- vec.start.point + length(new.genes) -1
          sample.gene.names[vec.start.point:vec.end.point] <- new.genes
          #detected.genes <- rownames(tech.data[which(tech.data[,cell2] > 0),])
          #sample.gene.names <- c(sample.gene.names,detected.genes)
          }
        sample.gene.names <- na.omit(sample.gene.names)
        num.uniq.genes <- length(unique(sample.gene.names))
        sample.gene.numbers[i] <- num.uniq.genes
      }
      sample.gene.numbers <- na.omit(sample.gene.numbers)
      print(paste(cell," = ", mean(sample.gene.numbers), sep = ""))
      tech.cumul.gene.numbers[cell] <- mean(sample.gene.numbers)
    }
    Bcells.plot.df <- rbind(Bcells.plot.df, data.frame(Cumul= tech.cumul.gene.numbers, tech= rep(tech, length(tech.cumul.gene.numbers)), cell.num= seq(1:length(tech.cumul.gene.numbers))))
}

save(Bcells.plot.df, "/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene/Cumulative_gene_dist_DS10K_Bcells_dataPlot.RData")

pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene/Cumulative_gene_dist_DS10K_Bcells.pdf")
ggplot(data=Bcells.plot.df, aes(x=cell.num, y=Cumul, group=tech)) + geom_point(aes(color=tech)) + geom_smooth(aes(color=tech))
dev.off()
print("Bcells cumulative done!")
