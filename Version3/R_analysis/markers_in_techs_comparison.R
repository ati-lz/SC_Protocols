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

# End ####

# Taking out Monocytes cells of each technique ####
MARSseq.Monocytes <- names(MARSseq.obj@ident)[which(MARSseq.obj@ident == "CD14+ Monocytes")]
#if(length(MARSseq.Monocytes) > 50){ MARSseq.Monocytes <- sample(MARSseq.Monocytes,50)}
CELseq2.Monocytes <- names(CELseq2.obj@ident)[which(CELseq2.obj@ident == "CD14+ Monocytes")]
#if(length(CELseq2.Monocytes) > 50){ CELseq2.Monocytes <- sample(CELseq2.Monocytes,50)}
QUARTZseq.Monocytes <- names(QUARTZseq.obj@ident)[which(QUARTZseq.obj@ident == "CD14+ Monocytes")]
#if(length(QUARTZseq.Monocytes) > 50){ QUARTZseq.Monocytes <- sample(QUARTZseq.Monocytes,50)}
Dropseq.Monocytes <- names(Dropseq.obj@ident)[which(Dropseq.obj@ident == "CD14+ Monocytes")]
#if(length(Dropseq.Monocytes) > 50){ Dropseq.Monocytes <- sample(Dropseq.Monocytes,50)}
SCRBseq.Monocytes <- names(SCRBseq.obj@ident)[which(SCRBseq.obj@ident == "CD14+ Monocytes")]
#if(length(SCRBseq.Monocytes) > 50){ SCRBseq.Monocytes <- sample(SCRBseq.Monocytes,50)}
X10Scilife.Monocytes <- names(X10Scilife.obj@ident)[which(X10Scilife.obj@ident == "CD14+ and FCGR3A+ Monocytes")]
#if(length(X10Scilife.Monocytes) > 50){ X10Scilife.Monocytes <- sample(X10Scilife.Monocytes,50)}
X10Nuclei.Monocytes <- names(X10Nuclei.obj@ident)[which(X10Nuclei.obj@ident == "CD14+ Monocytes")]
#if(length(X10Nuclei.Monocytes) > 50){ X10Nuclei.Monocytes <- sample(X10Nuclei.Monocytes,50)}
ICELL8.Monocytes <- names(ICELL8.obj@ident)[which(ICELL8.obj@ident == "CD14+ and FCGR3A+ Monocytes")]
#if(length(ICELL8.Monocytes) > 50){ ICELL8.Monocytes <- sample(ICELL8.Monocytes,50)}

# End ####

# Taking out Bcells cells of each technique ####
MARSseq.Bcells <- names(MARSseq.obj@ident)[which(MARSseq.obj@ident == "B cells")]
#if(length(MARSseq.Bcells) > 50){ MARSseq.Bcells <- sample(MARSseq.Bcells,50)}
CELseq2.Bcells <- names(CELseq2.obj@ident)[which(CELseq2.obj@ident == "B cells")]
#if(length(CELseq2.Bcells) > 50){ CELseq2.Bcells <- sample(CELseq2.Bcells,50)}
QUARTZseq.Bcells <- names(QUARTZseq.obj@ident)[which(QUARTZseq.obj@ident == "B cells")]
#if(length(QUARTZseq.Bcells) > 50){ QUARTZseq.Bcells <- sample(QUARTZseq.Bcells,50)}
Dropseq.Bcells <- names(Dropseq.obj@ident)[which(Dropseq.obj@ident == "B cells")]
#if(length(Dropseq.Bcells) > 50){ Dropseq.Bcells <- sample(Dropseq.Bcells,50)}
SCRBseq.Bcells <- names(SCRBseq.obj@ident)[which(SCRBseq.obj@ident == "B cells")]
#if(length(SCRBseq.Bcells) > 50){ SCRBseq.Bcells <- sample(SCRBseq.Bcells,50)}
X10Scilife.Bcells <- names(X10Scilife.obj@ident)[which(X10Scilife.obj@ident == "B cells")]
#if(length(X10Scilife.Bcells) > 50){ X10Scilife.Bcells <- sample(X10Scilife.Bcells,50)}
X10Nuclei.Bcells <- names(X10Nuclei.obj@ident)[which(X10Nuclei.obj@ident == "B cells")]
#if(length(X10Nuclei.Bcells) > 50){ X10Nuclei.Bcells <- sample(X10Nuclei.Bcells,50)}
ICELL8.Bcells <- names(ICELL8.obj@ident)[which(ICELL8.obj@ident == "B cells")]
#if(length(ICELL8.Bcells) > 50){ ICELL8.Bcells <- sample(ICELL8.Bcells,50)}

# End ####


# Loding Downsampleded data ####
load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/MARSseq.hsap.full.SCE.jointDSmat.Robj")
MARSseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
MARSseq.DS.UMI <- MARSseq.DS$UMI


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/CELseq2.hsap.full.SCE.jointDSmat.Robj")
CELseq2.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
CELseq2.DS.UMI <- CELseq2.DS$UMI


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/QUARTZseq.hsap.full.SCE.jointDSmat.Robj")
QUARTZseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
QUARTZseq.DS.UMI <- QUARTZseq.DS$UMI


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Dropseq.hsap.full.SCE.jointDSmat.Robj")
Dropseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
Dropseq.DS.UMI <- Dropseq.DS$UMI


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/SCRBseq.hsap.full.SCE.jointDSmat.Robj")
SCRBseq.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
SCRBseq.DS.UMI <- SCRBseq.DS$UMI


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/10XScilife.hsap.full.SCE.jointDSmat.Robj")
X10Scilife.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
X10Scilife.DS.UMI <- X10Scilife.DS$UMI


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/Nuclei10X.hsap.full.SCE.jointDSmat.Robj")
X10Nuclei.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
X10Nuclei.DS.UMI <- X10Nuclei.DS$UMI


load("/project/devel/alafzi/SC_Protocols/Downsampling/DS_secondRound/ICELL8.hsap.full.SCE.jointDSmat.Robj")
ICELL8.DS <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ICELL8.DS.UMI <- ICELL8.DS$UMI


# End ####

# preparing  B cells downsampleded matrix and plot
MARSseq.10K.UMI.mat <- mapIDs(MARSseq.DS.UMI$downsampled_10000, "hsap")
MARSseq.Bcells.common <- intersect(colnames(MARSseq.10K.UMI.mat),MARSseq.Bcells)

CELseq2.10K.UMI.mat <- mapIDs(CELseq2.DS.UMI$downsampled_10000, "hsap")
colnames(CELseq2.10K.UMI.mat) <- gsub(x = colnames(CELseq2.10K.UMI.mat), pattern = "\\.", replacement = "_")
CELseq2.Bcells.common <- intersect(colnames(CELseq2.10K.UMI.mat),CELseq2.Bcells)

QUARTZseq.10K.UMI.mat <- mapIDs(QUARTZseq.DS.UMI$downsampled_10000, "hsap")
QUARTZseq.Bcells.common <- intersect(colnames(QUARTZseq.10K.UMI.mat),QUARTZseq.Bcells)

Dropseq.10K.UMI.mat <- mapIDs(Dropseq.DS.UMI$downsampled_10000, "hsap")
Dropseq.Bcells.common <- intersect(colnames(Dropseq.10K.UMI.mat),Dropseq.Bcells)

SCRBseq.10K.UMI.mat <- mapIDs(SCRBseq.DS.UMI$downsampled_10000, "hsap")
SCRBseq.Bcells.common <- intersect(colnames(SCRBseq.10K.UMI.mat),SCRBseq.Bcells)

X10Scilife.10K.UMI.mat <- mapIDs(X10Scilife.DS.UMI$downsampled_10000, "hsap")
X10Scilife.Bcells.common <- intersect(colnames(X10Scilife.10K.UMI.mat),X10Scilife.Bcells)

X10Nuclei.10K.UMI.mat <- mapIDs(X10Nuclei.DS.UMI$downsampled_10000, "hsap")
X10Nuclei.Bcells.common <- intersect(colnames(X10Nuclei.10K.UMI.mat),X10Nuclei.Bcells)

ICELL8.10K.UMI.mat <- mapIDs(ICELL8.DS.UMI$downsampled_10000, "hsap")
ICELL8.Bcells.common <- intersect(colnames(ICELL8.10K.UMI.mat),ICELL8.Bcells)

#Separating cluster specific markers ####
load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/gene_cl.ref.RData")
ref.markers <- gene_cl.ref
Bcells.common <- Reduce(intersect, list(ref.markers[["B cells"]], rownames(MARSseq.10K.UMI.mat), rownames(CELseq2.10K.UMI.mat), rownames(QUARTZseq.10K.UMI.mat), rownames(Dropseq.10K.UMI.mat), rownames(SCRBseq.10K.UMI.mat) , rownames(X10Scilife.10K.UMI.mat), rownames(X10Nuclei.10K.UMI.mat), rownames(ICELL8.10K.UMI.mat)))
print(paste("length common Bcell markers = ", length(Bcells.common), sep=""))
#Monocytes.common <- Reduce(intersect, list(ref.markers[["CD14+ Monocytes"]], rownames(MARSseq.10K.UMI.mat), rownames(CELseq2.10K.UMI.mat), rownames(QUARTZseq.10K.UMI.mat), rownames(SCRBseq.10K.UMI.mat)))

MARSseq.10K.UMI.Bcell.mat <- MARSseq.10K.UMI.mat[Bcells.common, MARSseq.Bcells.common]
CELseq2.10K.UMI.Bcell.mat <- CELseq2.10K.UMI.mat[Bcells.common, CELseq2.Bcells.common]
SCRBseq.10K.UMI.Bcell.mat <- SCRBseq.10K.UMI.mat[Bcells.common, SCRBseq.Bcells.common]
QUARTZseq.10K.UMI.Bcell.mat <- QUARTZseq.10K.UMI.mat[Bcells.common, QUARTZseq.Bcells.common]
Dropseq.10K.UMI.Bcell.mat <- Dropseq.10K.UMI.mat[Bcells.common, Dropseq.Bcells.common]
X10Scilife.10K.UMI.Bcell.mat <- X10Scilife.10K.UMI.mat[Bcells.common, X10Scilife.Bcells.common]
X10Nuclei.10K.UMI.Bcell.mat <- X10Nuclei.10K.UMI.mat[Bcells.common, X10Nuclei.Bcells.common]
ICELL8.10K.UMI.Bcell.mat <- ICELL8.10K.UMI.mat[Bcells.common, ICELL8.Bcells.common]

library(dplyr)
library(ggplot2)
library(gplots)

Bcell.heatmap.df <- log(as.matrix(bind_cols(MARSseq.10K.UMI.Bcell.mat, CELseq2.10K.UMI.Bcell.mat, SCRBseq.10K.UMI.Bcell.mat, QUARTZseq.10K.UMI.Bcell.mat, Dropseq.10K.UMI.Bcell.mat, X10Scilife.10K.UMI.Bcell.mat, X10Nuclei.10K.UMI.Bcell.mat, ICELL8.10K.UMI.Bcell.mat)) + 1)
rownames(Bcell.heatmap.df) <- Bcells.common
col.separators = c(ncol(MARSseq.10K.UMI.Bcell.mat), ncol(MARSseq.10K.UMI.Bcell.mat) + ncol(CELseq2.10K.UMI.Bcell.mat),
 ncol(MARSseq.10K.UMI.Bcell.mat) + ncol(CELseq2.10K.UMI.Bcell.mat) + ncol(SCRBseq.10K.UMI.Bcell.mat),
 ncol(MARSseq.10K.UMI.Bcell.mat) + ncol(CELseq2.10K.UMI.Bcell.mat) + ncol(SCRBseq.10K.UMI.Bcell.mat) + ncol(QUARTZseq.10K.UMI.Bcell.mat),
 ncol(MARSseq.10K.UMI.Bcell.mat) + ncol(CELseq2.10K.UMI.Bcell.mat) + ncol(SCRBseq.10K.UMI.Bcell.mat) + ncol(QUARTZseq.10K.UMI.Bcell.mat) + ncol(Dropseq.10K.UMI.Bcell.mat),
 ncol(MARSseq.10K.UMI.Bcell.mat) + ncol(CELseq2.10K.UMI.Bcell.mat) + ncol(SCRBseq.10K.UMI.Bcell.mat) + ncol(QUARTZseq.10K.UMI.Bcell.mat) + ncol(Dropseq.10K.UMI.Bcell.mat) + ncol(X10Scilife.10K.UMI.Bcell.mat),
 ncol(MARSseq.10K.UMI.Bcell.mat) + ncol(CELseq2.10K.UMI.Bcell.mat) + ncol(SCRBseq.10K.UMI.Bcell.mat) + ncol(QUARTZseq.10K.UMI.Bcell.mat) + ncol(Dropseq.10K.UMI.Bcell.mat) + ncol(X10Scilife.10K.UMI.Bcell.mat) + ncol(X10Nuclei.10K.UMI.Bcell.mat))

col.sep.color = c(rep("slateblue1", ncol(MARSseq.10K.UMI.Bcell.mat)), rep("orchid1",ncol(CELseq2.10K.UMI.Bcell.mat)), rep("red", ncol(SCRBseq.10K.UMI.Bcell.mat)), rep("olivedrab3", ncol(QUARTZseq.10K.UMI.Bcell.mat)), rep("purple", ncol(Dropseq.10K.UMI.Bcell.mat)), rep("grey", ncol(X10Scilife.10K.UMI.Bcell.mat)), rep("orange", ncol(X10Nuclei.10K.UMI.Bcell.mat)), rep("yellow", ncol(ICELL8.10K.UMI.Bcell.mat)))
my_palette <- colorRampPalette(c("steelblue2","yellow", "orangered2"))(n = 299)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwide_DS_analysis/Markers_comparison_Bcells.pdf")
heatmap.2(Bcell.heatmap.df, Rowv = F, Colv = F, trace = "none", 
          colsep = col.separators,sepcolor = "black",sepwidth = c(0.8,0.8),
          ColSideColors= col.sep.color, labCol = F, col = my_palette, cexRow = 0.3, main = "B Cells markers downsampled to 10K")
legend(0.7,1.1,legend=c("MARSseq","CELseq2","SCRBseq", "QUARTZseq", "Dropseq", "X10Scilife", "X10Nuclei", "ICELL8"),
       fill=c("slateblue1", "orchid1", "red", "olivedrab3", "purple", "grey", "orange", "yellow"),border=FALSE, bty="n", y.intersp = 0.7, cex=0.7, xpd = T)

dev.off()

vec.Bcell.matrices <- list(log(as.matrix(MARSseq.10K.UMI.Bcell.mat) +1), log(as.matrix(CELseq2.10K.UMI.Bcell.mat)+1), log(as.matrix(SCRBseq.10K.UMI.Bcell.mat)+1), log(as.matrix(QUARTZseq.10K.UMI.Bcell.mat)+1), log(as.matrix(Dropseq.10K.UMI.Bcell.mat)+1), log(as.matrix(X10Scilife.10K.UMI.Bcell.mat)+), log(as.matrix(X10Nuclei.10K.UMI.Bcell.mat)+1), log(as.matrix(ICELL8.10K.UMI.Bcell.mat)+1))
mean.tech.matrices.Bcell <- unlist(lapply(vec.Bcell.matrices, function(x) mean(x)))
names(mean.tech.matrices.Bcell) <- c("MARSseq", "CELseq2", "SCRBseq", "QUARTZseq","Dropseq", "X10Scilife", "X10Nuclei", "ICELL8")



# preparing  Monocytes downsampleded matrix and plot
MARSseq.10K.UMI.mat <- mapIDs(MARSseq.DS.UMI$downsampled_10000, "hsap")
MARSseq.Monocytes.common <- intersect(colnames(MARSseq.10K.UMI.mat),MARSseq.Monocytes)

CELseq2.10K.UMI.mat <- mapIDs(CELseq2.DS.UMI$downsampled_10000, "hsap")
colnames(CELseq2.10K.UMI.mat) <- gsub(x = colnames(CELseq2.10K.UMI.mat), pattern = "\\.", replacement = "_")
CELseq2.Monocytes.common <- intersect(colnames(CELseq2.10K.UMI.mat),CELseq2.Monocytes)

QUARTZseq.10K.UMI.mat <- mapIDs(QUARTZseq.DS.UMI$downsampled_10000, "hsap")
QUARTZseq.Monocytes.common <- intersect(colnames(QUARTZseq.10K.UMI.mat),QUARTZseq.Monocytes)

Dropseq.10K.UMI.mat <- mapIDs(Dropseq.DS.UMI$downsampled_10000, "hsap")
Dropseq.Monocytes.common <- intersect(colnames(Dropseq.10K.UMI.mat),Dropseq.Monocytes)

SCRBseq.10K.UMI.mat <- mapIDs(SCRBseq.DS.UMI$downsampled_10000, "hsap")
SCRBseq.Monocytes.common <- intersect(colnames(SCRBseq.10K.UMI.mat),SCRBseq.Monocytes)

X10Scilife.10K.UMI.mat <- mapIDs(X10Scilife.DS.UMI$downsampled_10000, "hsap")
X10Scilife.Monocytes.common <- intersect(colnames(X10Scilife.10K.UMI.mat),X10Scilife.Monocytes)

X10Nuclei.10K.UMI.mat <- mapIDs(X10Nuclei.DS.UMI$downsampled_10000, "hsap")
X10Nuclei.Monocytes.common <- intersect(colnames(X10Nuclei.10K.UMI.mat),X10Nuclei.Monocytes)

ICELL8.10K.UMI.mat <- mapIDs(ICELL8.DS.UMI$downsampled_10000, "hsap")
ICELL8.Monocytes.common <- intersect(colnames(ICELL8.10K.UMI.mat),ICELL8.Monocytes)

#Separating cluster specific markers ####
load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/gene_cl.ref.RData")
ref.markers <- gene_cl.ref
Monocytes.common <- Reduce(intersect, list(ref.markers[["CD14+ Monocytes"]], rownames(MARSseq.10K.UMI.mat), rownames(CELseq2.10K.UMI.mat), rownames(QUARTZseq.10K.UMI.mat), rownames(Dropseq.10K.UMI.mat), rownames(SCRBseq.10K.UMI.mat) , rownames(X10Scilife.10K.UMI.mat), rownames(X10Nuclei.10K.UMI.mat), rownames(ICELL8.10K.UMI.mat)))
print(paste("length common CD14+ Monocytes markers = ", length(Monocytes.common), sep=""))
#Monocytes.common <- Reduce(intersect, list(ref.markers[["CD14+ Monocytes"]], rownames(MARSseq.10K.UMI.mat), rownames(CELseq2.10K.UMI.mat), rownames(QUARTZseq.10K.UMI.mat), rownames(SCRBseq.10K.UMI.mat)))

MARSseq.10K.UMI.Monocytes.mat <- MARSseq.10K.UMI.mat[Monocytes.common, MARSseq.Monocytes.common]
CELseq2.10K.UMI.Monocytes.mat <- CELseq2.10K.UMI.mat[Monocytes.common, CELseq2.Monocytes.common]
SCRBseq.10K.UMI.Monocytes.mat <- SCRBseq.10K.UMI.mat[Monocytes.common, SCRBseq.Monocytes.common]
QUARTZseq.10K.UMI.Monocytes.mat <- QUARTZseq.10K.UMI.mat[Monocytes.common, QUARTZseq.Monocytes.common]
Dropseq.10K.UMI.Monocytes.mat <- Dropseq.10K.UMI.mat[Monocytes.common, Dropseq.Monocytes.common]
X10Scilife.10K.UMI.Monocytes.mat <- X10Scilife.10K.UMI.mat[Monocytes.common, X10Scilife.Monocytes.common]
X10Nuclei.10K.UMI.Monocytes.mat <- X10Nuclei.10K.UMI.mat[Monocytes.common, X10Nuclei.Monocytes.common]
ICELL8.10K.UMI.Monocytes.mat <- ICELL8.10K.UMI.mat[Monocytes.common, ICELL8.Monocytes.common]

library(dplyr)
library(ggplot2)
library(gplots)

Monocytes.heatmap.df <- log(as.matrix(bind_cols(MARSseq.10K.UMI.Monocytes.mat, CELseq2.10K.UMI.Monocytes.mat, SCRBseq.10K.UMI.Monocytes.mat, QUARTZseq.10K.UMI.Monocytes.mat, Dropseq.10K.UMI.Monocytes.mat, X10Scilife.10K.UMI.Monocytes.mat, X10Nuclei.10K.UMI.Monocytes.mat, ICELL8.10K.UMI.Monocytes.mat)) + 1)
rownames(Monocytes.heatmap.df) <- Monocytes.common
col.separators = c(ncol(MARSseq.10K.UMI.Monocytes.mat), ncol(MARSseq.10K.UMI.Monocytes.mat) + ncol(CELseq2.10K.UMI.Monocytes.mat),
 ncol(MARSseq.10K.UMI.Monocytes.mat) + ncol(CELseq2.10K.UMI.Monocytes.mat) + ncol(SCRBseq.10K.UMI.Monocytes.mat),
 ncol(MARSseq.10K.UMI.Monocytes.mat) + ncol(CELseq2.10K.UMI.Monocytes.mat) + ncol(SCRBseq.10K.UMI.Monocytes.mat) + ncol(QUARTZseq.10K.UMI.Monocytes.mat),
 ncol(MARSseq.10K.UMI.Monocytes.mat) + ncol(CELseq2.10K.UMI.Monocytes.mat) + ncol(SCRBseq.10K.UMI.Monocytes.mat) + ncol(QUARTZseq.10K.UMI.Monocytes.mat) + ncol(Dropseq.10K.UMI.Monocytes.mat),
 ncol(MARSseq.10K.UMI.Monocytes.mat) + ncol(CELseq2.10K.UMI.Monocytes.mat) + ncol(SCRBseq.10K.UMI.Monocytes.mat) + ncol(QUARTZseq.10K.UMI.Monocytes.mat) + ncol(Dropseq.10K.UMI.Monocytes.mat) + ncol(X10Scilife.10K.UMI.Monocytes.mat),
 ncol(MARSseq.10K.UMI.Monocytes.mat) + ncol(CELseq2.10K.UMI.Monocytes.mat) + ncol(SCRBseq.10K.UMI.Monocytes.mat) + ncol(QUARTZseq.10K.UMI.Monocytes.mat) + ncol(Dropseq.10K.UMI.Monocytes.mat) + ncol(X10Scilife.10K.UMI.Monocytes.mat) + ncol(X10Nuclei.10K.UMI.Monocytes.mat))

col.sep.color = c(rep("slateblue1", ncol(MARSseq.10K.UMI.Monocytes.mat)), rep("orchid1",ncol(CELseq2.10K.UMI.Monocytes.mat)), rep("red", ncol(SCRBseq.10K.UMI.Monocytes.mat)), rep("olivedrab3", ncol(QUARTZseq.10K.UMI.Monocytes.mat)), rep("purple", ncol(Dropseq.10K.UMI.Monocytes.mat)), rep("grey", ncol(X10Scilife.10K.UMI.Monocytes.mat)), rep("orange", ncol(X10Nuclei.10K.UMI.Monocytes.mat)), rep("yellow", ncol(ICELL8.10K.UMI.Monocytes.mat)))
my_palette <- colorRampPalette(c("steelblue2","yellow", "orangered2"))(n = 299)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwide_DS_analysis/Markers_comparison_Monocytes.pdf")
heatmap.2(Monocytes.heatmap.df, Rowv = F, Colv = F, trace = "none", 
          colsep = col.separators,sepcolor = "black",sepwidth = c(0.8,0.8),
          ColSideColors= col.sep.color, labCol = F, col = my_palette, cexRow = 0.3, main = "Monocytes markers downsampled to 10K")
legend(0.7,1.1,legend=c("MARSseq","CELseq2","SCRBseq", "QUARTZseq", "Dropseq", "X10Scilife", "X10Nuclei", "ICELL8"),
       fill=c("slateblue1", "orchid1", "red", "olivedrab3", "purple", "grey", "orange", "yellow"),border=FALSE, bty="n", y.intersp = 0.7, cex=0.7, xpd = T)

dev.off()

vec.Monocytes.matrices <- list(log(as.matrix(MARSseq.10K.UMI.Monocytes.mat)+1), log(as.matrix(CELseq2.10K.UMI.Monocytes.mat)+1), log(as.matrix(SCRBseq.10K.UMI.Monocytes.mat)+1), log(as.matrix(QUARTZseq.10K.UMI.Monocytes.mat)+1), log(as.matrix(Dropseq.10K.UMI.Monocytes.mat)+1), log(as.matrix(X10Scilife.10K.UMI.Monocytes.mat)+1), log(as.matrix(X10Nuclei.10K.UMI.Monocytes.mat)+1), log(as.matrix(ICELL8.10K.UMI.Monocytes.mat)+1))
mean.tech.matrices.Monocytes <- unlist(lapply(vec.Monocytes.matrices, function(x) mean(x)))
names(mean.tech.matrices.Monocytes) <- c("MARSseq", "CELseq2", "SCRBseq", "QUARTZseq","Dropseq", "X10Scilife", "X10Nuclei", "ICELL8")



# preparing  Hek cells downsampleded matrix and plot
MARSseq.20K.UMI.mat <- mapIDs(MARSseq.DS.UMI$downsampled_20000, "hsap")
MARSseq.HEK.common <- intersect(colnames(MARSseq.20K.UMI.mat),MARSseq.HEK)

CELseq2.20K.UMI.mat <- mapIDs(CELseq2.DS.UMI$downsampled_20000, "hsap")
colnames(CELseq2.20K.UMI.mat) <- gsub(x = colnames(CELseq2.20K.UMI.mat), pattern = "\\.", replacement = "_")
CELseq2.HEK.common <- intersect(colnames(CELseq2.20K.UMI.mat),CELseq2.HEK)

QUARTZseq.20K.UMI.mat <- mapIDs(QUARTZseq.DS.UMI$downsampled_20000, "hsap")
QUARTZseq.HEK.common <- intersect(colnames(QUARTZseq.20K.UMI.mat),QUARTZseq.HEK)

Dropseq.20K.UMI.mat <- mapIDs(Dropseq.DS.UMI$downsampled_20000, "hsap")
Dropseq.HEK.common <- intersect(colnames(Dropseq.20K.UMI.mat),Dropseq.HEK)

SCRBseq.20K.UMI.mat <- mapIDs(SCRBseq.DS.UMI$downsampled_20000, "hsap")
SCRBseq.HEK.common <- intersect(colnames(SCRBseq.20K.UMI.mat),SCRBseq.HEK)

X10Scilife.20K.UMI.mat <- mapIDs(X10Scilife.DS.UMI$downsampled_20000, "hsap")
X10Scilife.HEK.common <- intersect(colnames(X10Scilife.20K.UMI.mat),X10Scilife.HEK)

X10Nuclei.20K.UMI.mat <- mapIDs(X10Nuclei.DS.UMI$downsampled_20000, "hsap")
X10Nuclei.HEK.common <- intersect(colnames(X10Nuclei.20K.UMI.mat),X10Nuclei.HEK)

ICELL8.20K.UMI.mat <- mapIDs(ICELL8.DS.UMI$downsampled_20000, "hsap")
ICELL8.HEK.common <- intersect(colnames(ICELL8.20K.UMI.mat),ICELL8.HEK)

#Separating cluster specific markers ####
load("/project/devel/alafzi/SC_Protocols/Version3/Seurat_objects/gene_cl.ref.RData")
ref.markers <- gene_cl.ref
HEK.common <- Reduce(intersect, list(ref.markers[["HEK cells"]], rownames(MARSseq.20K.UMI.mat), rownames(CELseq2.20K.UMI.mat), rownames(QUARTZseq.20K.UMI.mat), rownames(Dropseq.20K.UMI.mat), rownames(SCRBseq.20K.UMI.mat) , rownames(X10Scilife.20K.UMI.mat), rownames(X10Nuclei.20K.UMI.mat), rownames(ICELL8.20K.UMI.mat)))
print(paste("length common HEK markers = ", length(HEK.common), sep=""))
#HEK.common <- Reduce(intersect, list(ref.markers[["CD14+ HEK"]], rownames(MARSseq.20K.UMI.mat), rownames(CELseq2.20K.UMI.mat), rownames(QUARTZseq.20K.UMI.mat), rownames(SCRBseq.20K.UMI.mat)))

MARSseq.20K.UMI.HEK.mat <- MARSseq.20K.UMI.mat[HEK.common, MARSseq.HEK.common]
CELseq2.20K.UMI.HEK.mat <- CELseq2.20K.UMI.mat[HEK.common, CELseq2.HEK.common]
SCRBseq.20K.UMI.HEK.mat <- SCRBseq.20K.UMI.mat[HEK.common, SCRBseq.HEK.common]
QUARTZseq.20K.UMI.HEK.mat <- QUARTZseq.20K.UMI.mat[HEK.common, QUARTZseq.HEK.common]
Dropseq.20K.UMI.HEK.mat <- Dropseq.20K.UMI.mat[HEK.common, Dropseq.HEK.common]
X10Scilife.20K.UMI.HEK.mat <- X10Scilife.20K.UMI.mat[HEK.common, X10Scilife.HEK.common]
X10Nuclei.20K.UMI.HEK.mat <- X10Nuclei.20K.UMI.mat[HEK.common, X10Nuclei.HEK.common]
ICELL8.20K.UMI.HEK.mat <- ICELL8.20K.UMI.mat[HEK.common, ICELL8.HEK.common]

library(dplyr)
library(ggplot2)
library(gplots)

HEK.heatmap.df <- log(as.matrix(bind_cols(MARSseq.20K.UMI.HEK.mat, CELseq2.20K.UMI.HEK.mat, SCRBseq.20K.UMI.HEK.mat, QUARTZseq.20K.UMI.HEK.mat, Dropseq.20K.UMI.HEK.mat, X10Scilife.20K.UMI.HEK.mat, X10Nuclei.20K.UMI.HEK.mat, ICELL8.20K.UMI.HEK.mat)) + 1)
rownames(HEK.heatmap.df) <- HEK.common
col.separators = c(ncol(MARSseq.20K.UMI.HEK.mat), ncol(MARSseq.20K.UMI.HEK.mat) + ncol(CELseq2.20K.UMI.HEK.mat),
 ncol(MARSseq.20K.UMI.HEK.mat) + ncol(CELseq2.20K.UMI.HEK.mat) + ncol(SCRBseq.20K.UMI.HEK.mat),
 ncol(MARSseq.20K.UMI.HEK.mat) + ncol(CELseq2.20K.UMI.HEK.mat) + ncol(SCRBseq.20K.UMI.HEK.mat) + ncol(QUARTZseq.20K.UMI.HEK.mat),
 ncol(MARSseq.20K.UMI.HEK.mat) + ncol(CELseq2.20K.UMI.HEK.mat) + ncol(SCRBseq.20K.UMI.HEK.mat) + ncol(QUARTZseq.20K.UMI.HEK.mat) + ncol(Dropseq.20K.UMI.HEK.mat),
 ncol(MARSseq.20K.UMI.HEK.mat) + ncol(CELseq2.20K.UMI.HEK.mat) + ncol(SCRBseq.20K.UMI.HEK.mat) + ncol(QUARTZseq.20K.UMI.HEK.mat) + ncol(Dropseq.20K.UMI.HEK.mat) + ncol(X10Scilife.20K.UMI.HEK.mat),
 ncol(MARSseq.20K.UMI.HEK.mat) + ncol(CELseq2.20K.UMI.HEK.mat) + ncol(SCRBseq.20K.UMI.HEK.mat) + ncol(QUARTZseq.20K.UMI.HEK.mat) + ncol(Dropseq.20K.UMI.HEK.mat) + ncol(X10Scilife.20K.UMI.HEK.mat) + ncol(X10Nuclei.20K.UMI.HEK.mat))

col.sep.color = c(rep("slateblue1", ncol(MARSseq.20K.UMI.HEK.mat)), rep("orchid1",ncol(CELseq2.20K.UMI.HEK.mat)), rep("red", ncol(SCRBseq.20K.UMI.HEK.mat)), rep("olivedrab3", ncol(QUARTZseq.20K.UMI.HEK.mat)), rep("purple", ncol(Dropseq.20K.UMI.HEK.mat)), rep("grey", ncol(X10Scilife.20K.UMI.HEK.mat)), rep("orange", ncol(X10Nuclei.20K.UMI.HEK.mat)), rep("yellow", ncol(ICELL8.20K.UMI.HEK.mat)))
my_palette <- colorRampPalette(c("steelblue2","yellow", "orangered2"))(n = 299)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwide_DS_analysis/Markers_comparison_HEK.pdf")
heatmap.2(HEK.heatmap.df, Rowv = F, Colv = F, trace = "none", 
          colsep = col.separators,sepcolor = "black",sepwidth = c(0.8,0.8),
          ColSideColors= col.sep.color, labCol = F, col = my_palette, cexRow = 0.3, main = "HEK markers downsampled to 20K")
legend(0.7,1.1,legend=c("MARSseq","CELseq2","SCRBseq", "QUARTZseq", "Dropseq", "X10Scilife", "X10Nuclei", "ICELL8"),
       fill=c("slateblue1", "orchid1", "red", "olivedrab3", "purple", "grey", "orange", "yellow"),border=FALSE, bty="n", y.intersp = 0.7, cex=0.7, xpd = T)

dev.off()


vec.HEK.matrices <- list(log(as.matrix(MARSseq.20K.UMI.HEK.mat)+1), log(as.matrix(CELseq2.20K.UMI.HEK.mat)+1), log(as.matrix(SCRBseq.20K.UMI.HEK.mat)+1), log(as.matrix(QUARTZseq.20K.UMI.HEK.mat)+1), log(as.matrix(Dropseq.20K.UMI.HEK.mat)+1), log(as.matrix(X10Scilife.20K.UMI.HEK.mat)+1), log(as.matrix(X10Nuclei.20K.UMI.HEK.mat)+1), log(as.matrix(ICELL8.20K.UMI.HEK.mat)+1))
mean.tech.matrices.HEK <- unlist(lapply(vec.HEK.matrices, function(x) mean(x)))
names(mean.tech.matrices.HEK) <- c("MARSseq", "CELseq2", "SCRBseq", "QUARTZseq","Dropseq", "X10Scilife", "X10Nuclei", "ICELL8")

#The overall mean hetamap per tech per celltype
whole.mean.df <- cbind(Monocytes=mean.tech.matrices.Monocytes, Bcells=mean.tech.matrices.Bcell, HEK=mean.tech.matrices.HEK)
col.sep.color = c("slateblue1", "orchid1")#, "red", "olivedrab3")
my_palette <- colorRampPalette(c("yellow", "orangered2"))(n = 10)
pdf("/project/devel/alafzi/SC_Protocols/Version3/R_analysis/stepwide_DS_analysis/Markers_comparison_MeanPerTechPerCelltype.pdf")
heatmap.2(whole.mean.df, Rowv = F, Colv = F, trace = "none", col = my_palette, colsep = c(1,2),sepcolor = "black",sepwidth = c(0.02,0.02), cexRow = 1.5, cexCol = 1.5, srtCol= 0, adjCol= 0.5, margins = c(5,10)) #ColSideColors= col.sep.color,
dev.off()