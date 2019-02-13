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
Nuclei10X.metadata <- X10Nuclei.obj@meta.data

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
X10X2x5Half.metadata <- colData(X10X2x5Half.hsap.final)
#X10X2x5Half.hsap.fin.metadata <- X10X2x5Half.hsap.fin.metadata[rownames(X10X2x5Half.metadata),]

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
hsap.total.reads <- data.frame(MARSseq.metadata[, c("Library", "nTReads")], tech = rep("MARSseq", nrow(MARSseq.metadata)))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(CELseq2.metadata[, c("Library", "nTReads")], tech = rep("CELseq2", nrow(CELseq2.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(QUARTZseq.metadata[, c("Library", "nTReads")], tech = rep("QUARTZseq", nrow(QUARTZseq.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(Dropseq.metadata[, c("Library", "nTReads")], tech = rep("Dropseq", nrow(Dropseq.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(SCRBseq.metadata[, c("Library", "nTReads")], tech = rep("SCRBseq", nrow(SCRBseq.metadata))))
#hsap.total.reads <- rbind(hsap.total.reads, data.frame(SeqwellV1.metadata[, c("Library", "nTReads")], tech = rep("SeqwellV1", nrow(SeqwellV1.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(SeqwellV2.metadata[, c("Library", "nTReads")], tech = rep("SeqwellV2", nrow(SeqwellV2.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(Nuclei10X.metadata[, c("Library", "nTReads")], tech = rep("Nuclei10X", nrow(Nuclei10X.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(ICELL8.metadata[, c("Library", "nTReads")], tech = rep("ICELL8", nrow(ICELL8.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(ddSEQ.metadata[, c("Library", "nTReads")], tech = rep("ddSEQ", nrow(ddSEQ.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(ddSEQexp1.metadata[, c("Library", "nTReads")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(C1HTsmall.metadata[, c("Library", "nTReads")], tech = rep("C1HTsmall", nrow(C1HTsmall.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(C1HTmedium.metadata[, c("Library", "nTReads")], tech = rep("C1HTmedium", nrow(C1HTmedium.metadata))))
#hsap.total.reads <- rbind(hsap.total.reads, data.frame(X10x8x10K.metadata[, c("Library", "nTReads")], tech = rep("X10x8x10K", nrow(X10x8x10K.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(X10X2x5Half.metadata[, c("Library", "nTReads")], tech = rep("X10X2x5Half", nrow(X10X2x5Half.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(X10Scilife.metadata[, c("Library", "nTReads")], tech = rep("X10Scilife", nrow(X10Scilife.metadata))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(CB1.metadata[, c("Library", "nTReads")], tech = rep("CB1", nrow(CB1.metadata))))
pdf("all_techs_V2/hsap_total_number_Reads_all_V7.pdf")
#ggplot(data=hsap.total.reads, aes(x=Library, y=log(nTReads), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
#ggplot(data=hsap.total.reads, aes(x=tech, y=log(nTReads), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
ggplot(data=hsap.total.reads, aes(x=tech, y=nTReads, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(1000, 20000, 100000, 500000, 1000000)) + facet_grid(. ~ tech, scales = "free") + geom_hline(yintercept=c(1000, 20000, 100000, 500000, 1000000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
#ggplot(data=hsap.total.reads, aes(x=tech, y=nTReads, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(labels= scales::comma, breaks = c(1000, 20000, 100000, 500000, 1000000)) + facet_grid(. ~ tech, scales = "free") + geom_hline(yintercept=c(1000, 20000, 100000, 500000, 1000000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
dev.off()




#Total number of genes per tech per lib ====
print("Total number of genes per tech per lib ")
hsap.total.genes <- data.frame(MARSseq.metadata[, c("Library", "nGenes")], tech = rep("MARSseq", nrow(MARSseq.metadata)))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(CELseq2.metadata[, c("Library", "nGenes")], tech = rep("CELseq2", nrow(CELseq2.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(QUARTZseq.metadata[, c("Library", "nGenes")], tech = rep("QUARTZseq", nrow(QUARTZseq.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(Dropseq.metadata[, c("Library", "nGenes")], tech = rep("Dropseq", nrow(Dropseq.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(SCRBseq.metadata[, c("Library", "nGenes")], tech = rep("SCRBseq", nrow(SCRBseq.metadata))))
#hsap.total.genes <- rbind(hsap.total.genes, data.frame(SeqwellV1.metadata[, c("Library", "nGenes")], tech = rep("SeqwellV1", nrow(SeqwellV1.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(SeqwellV2.metadata[, c("Library", "nGenes")], tech = rep("SeqwellV2", nrow(SeqwellV2.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(Nuclei10X.metadata[, c("Library", "nGenes")], tech = rep("Nuclei10X", nrow(Nuclei10X.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(ICELL8.metadata[, c("Library", "nGenes")], tech = rep("ICELL8", nrow(ICELL8.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(ddSEQ.metadata[, c("Library", "nGenes")], tech = rep("ddSEQ", nrow(ddSEQ.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(ddSEQexp1.metadata[, c("Library", "nGenes")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(C1HTsmall.metadata[, c("Library", "nGenes")], tech = rep("C1HTsmall", nrow(C1HTsmall.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(C1HTmedium.metadata[, c("Library", "nGenes")], tech = rep("C1HTmedium", nrow(C1HTmedium.metadata))))
#hsap.total.genes <- rbind(hsap.total.genes, data.frame(X10x8x10K.metadata[, c("Library", "nGenes")], tech = rep("X10x8x10K", nrow(X10x8x10K.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(X10X2x5Half.metadata[, c("Library", "nGenes")], tech = rep("X10X2x5Half", nrow(X10X2x5Half.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(X10Scilife.metadata[, c("Library", "nGenes")], tech = rep("X10Scilife", nrow(X10Scilife.metadata))))
hsap.total.genes <- rbind(hsap.total.genes, data.frame(CB1.metadata[, c("Library", "nGenes")], tech = rep("CB1", nrow(CB1.metadata))))
pdf("all_techs_V2/hsap_total_number_Genes_all_V7.pdf")
#ggplot(data=hsap.total.genes, aes(x=Library, y=log(nGenes), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
#ggplot(data=hsap.total.genes, aes(x=tech, y=log(nGenes), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
ggplot(data=hsap.total.genes, aes(x=tech, y=nGenes, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(10,150,500, 1500, 5000)) + facet_grid(. ~ tech, scales = "free") + geom_hline(yintercept=c(10,150,500, 1500, 5000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
#ggplot(data=hsap.total.genes, aes(x=tech, y=nGenes, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(labels= scales::comma, breaks = c(10,150,500, 1500, 5000)) + facet_grid(. ~ tech, scales = "free") + geom_hline(yintercept=c(10,150,500, 1500, 5000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
dev.off()


# Total number of UMIs per tech per lib ====
print("Total number of UMIs per tech per lib")
hsap.total.UMIs <- data.frame(MARSseq.metadata[, c("Library", "nUMIs")], tech = rep("MARSseq", nrow(MARSseq.metadata)))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(CELseq2.metadata[, c("Library", "nUMIs")], tech = rep("CELseq2", nrow(CELseq2.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(QUARTZseq.metadata[, c("Library", "nUMIs")], tech = rep("QUARTZseq", nrow(QUARTZseq.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(Dropseq.metadata[, c("Library", "nUMIs")], tech = rep("Dropseq", nrow(Dropseq.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(SCRBseq.metadata[, c("Library", "nUMIs")], tech = rep("SCRBseq", nrow(SCRBseq.metadata))))
#hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(SeqwellV1.metadata[, c("Library", "nUMIs")], tech = rep("SeqwellV1", nrow(SeqwellV1.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(SeqwellV2.metadata[, c("Library", "nUMIs")], tech = rep("SeqwellV2", nrow(SeqwellV2.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(Nuclei10X.metadata[, c("Library", "nUMIs")], tech = rep("Nuclei10X", nrow(Nuclei10X.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(ICELL8.metadata[, c("Library", "nUMIs")], tech = rep("ICELL8", nrow(ICELL8.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(ddSEQ.metadata[, c("Library", "nUMIs")], tech = rep("ddSEQ", nrow(ddSEQ.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(ddSEQexp1.metadata[, c("Library", "nUMIs")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(C1HTsmall.metadata[, c("Library", "nUMIs")], tech = rep("C1HTsmall", nrow(C1HTsmall.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(C1HTmedium.metadata[, c("Library", "nUMIs")], tech = rep("C1HTmedium", nrow(C1HTmedium.metadata))))
#hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(X10x8x10K.metadata[, c("Library", "nUMIs")], tech = rep("X10x8x10K", nrow(X10x8x10K.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(X10X2x5Half.metadata[, c("Library", "nUMIs")], tech = rep("X10X2x5Half", nrow(X10X2x5Half.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(X10Scilife.metadata[, c("Library", "nUMIs")], tech = rep("X10Scilife", nrow(X10Scilife.metadata))))
hsap.total.UMIs <- rbind(hsap.total.UMIs, data.frame(CB1.metadata[, c("Library", "nUMIs")], tech = rep("CB1", nrow(CB1.metadata))))
pdf("all_techs_V2/hsap_total_number_UMIs_all_V7.pdf")
#ggplot(data=hsap.total.UMIs, aes(x=Library, y=log(nUMIs), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free")
#ggplot(data=hsap.total.UMIs, aes(x=tech, y=log(nUMIs), fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free")
ggplot(data=hsap.total.UMIs, aes(x=tech, y=nUMIs, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(200,2000,20000,200000, 1000000)) + facet_grid(. ~ tech, scales = "free")+ geom_hline(yintercept=c(200,2000,20000,200000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
#ggplot(data=hsap.total.UMIs, aes(x=tech, y=nUMIs, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) + scale_y_continuous(labels= scales::comma, breaks = c(200,2000,20000,200000, 1000000)) + facet_grid(. ~ tech, scales = "free")+ geom_hline(yintercept=c(200,2000,20000,200000), color = "orange", linetype = "dashed") # , limits = c(5000, NA)
dev.off()


# nTreads vs nUMI curves ====
print("nTreads vs nUMI curves")
hsap.nread.numi.df <- data.frame(MARSseq.metadata[, c("nTReads", "nUMIs")], tech = rep("MARSseq", nrow(MARSseq.metadata)))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(CELseq2.metadata[, c("nTReads", "nUMIs")], tech = rep("CELseq2", nrow(CELseq2.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(QUARTZseq.metadata[, c("nTReads", "nUMIs")], tech = rep("QUARTZseq", nrow(QUARTZseq.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(Dropseq.metadata[, c("nTReads", "nUMIs")], tech = rep("Dropseq", nrow(Dropseq.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(SCRBseq.metadata[, c("nTReads", "nUMIs")], tech = rep("SCRBseq", nrow(SCRBseq.metadata))))
#hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(SeqwellV1.metadata[, c("nTReads", "nUMIs")], tech = rep("SeqwellV1", nrow(SeqwellV1.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(SeqwellV2.metadata[, c("nTReads", "nUMIs")], tech = rep("SeqwellV2", nrow(SeqwellV2.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(Nuclei10X.metadata[, c("nTReads", "nUMIs")], tech = rep("Nuclei10X", nrow(Nuclei10X.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(ICELL8.metadata[, c("nTReads", "nUMIs")], tech = rep("ICELL8", nrow(ICELL8.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(ddSEQ.metadata[, c("nTReads", "nUMIs")], tech = rep("ddSEQ", nrow(ddSEQ.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(ddSEQexp1.metadata[, c("nTReads", "nUMIs")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(C1HTsmall.metadata[, c("nTReads", "nUMIs")], tech = rep("C1HTsmall", nrow(C1HTsmall.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(C1HTmedium.metadata[, c("nTReads", "nUMIs")], tech = rep("C1HTmedium", nrow(C1HTmedium.metadata))))
#hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(X10x8x10K.metadata[, c("nTReads", "nUMIs")], tech = rep("X10x8x10K", nrow(X10x8x10K.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(X10X2x5Half.metadata[, c("nTReads", "nUMIs")], tech = rep("X10X2x5Half", nrow(X10X2x5Half.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(X10Scilife.metadata[, c("nTReads", "nUMIs")], tech = rep("X10Scilife", nrow(X10Scilife.metadata))))
hsap.nread.numi.df <- rbind(hsap.nread.numi.df, data.frame(CB1.metadata[, c("nTReads", "nUMIs")], tech = rep("CB1", nrow(CB1.metadata))))

pdf("all_techs_V2/hsap_nReadsVSnUMIs_all_V7.pdf")
#ggplot(hsap.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(aes(color=tech))+ xlim(0, 100000)# + geom_line(aes(color=tech))#
ggplot(hsap.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(method="lm", formula=y~x, fill="grey",aes(color=tech))+ xlim(0, 100000)# + geom_line(aes(color=tech))#
#ggplot(hsap.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(method="lm", formula=y~log(x), fill="grey", aes(color=tech))+ xlim(0, 2000000)# + geom_line(aes(color=tech))#
#ggplot(hsap.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(method='nlsLM',formula=y ~ a*x^b, se=FALSE, method.args = list(start = list(a=1,b=1)), aes(color=tech))+ xlim(0, 1500000)# + geom_line(aes(color=tech))#
#ggplot(hsap.nread.numi.df, aes(x=nTReads, y=nUMIs, group=tech)) + geom_smooth(method='nls',formula=y ~ SSasympOff(x, A, lrc, c0), se=FALSE, aes(color=tech))+ xlim(0, 1500000)# + geom_line(aes(color=tech))#

dev.off()



# nTreads vs nGenes curves ====
print("nTreads vs nGenes curves")
hsap.nread.ngene.df <- data.frame(MARSseq.metadata[, c("nTReads", "nGenes")], tech = rep("MARSseq", nrow(MARSseq.metadata)))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(CELseq2.metadata[, c("nTReads", "nGenes")], tech = rep("CELseq2", nrow(CELseq2.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(QUARTZseq.metadata[, c("nTReads", "nGenes")], tech = rep("QUARTZseq", nrow(QUARTZseq.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(Dropseq.metadata[, c("nTReads", "nGenes")], tech = rep("Dropseq", nrow(Dropseq.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(SCRBseq.metadata[, c("nTReads", "nGenes")], tech = rep("SCRBseq", nrow(SCRBseq.metadata))))
#hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(SeqwellV1.metadata[, c("nTReads", "nGenes")], tech = rep("SeqwellV1", nrow(SeqwellV1.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(SeqwellV2.metadata[, c("nTReads", "nGenes")], tech = rep("SeqwellV2", nrow(SeqwellV2.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(Nuclei10X.metadata[, c("nTReads", "nGenes")], tech = rep("Nuclei10X", nrow(Nuclei10X.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(ICELL8.metadata[, c("nTReads", "nGenes")], tech = rep("ICELL8", nrow(ICELL8.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(ddSEQ.metadata[, c("nTReads", "nGenes")], tech = rep("ddSEQ", nrow(ddSEQ.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(ddSEQexp1.metadata[, c("nTReads", "nGenes")], tech = rep("ddSEQexp1", nrow(ddSEQexp1.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(C1HTsmall.metadata[, c("nTReads", "nGenes")], tech = rep("C1HTsmall", nrow(C1HTsmall.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(C1HTmedium.metadata[, c("nTReads", "nGenes")], tech = rep("C1HTmedium", nrow(C1HTmedium.metadata))))
#hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(X10x8x10K.metadata[, c("nTReads", "nGenes")], tech = rep("X10x8x10K", nrow(X10x8x10K.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(X10X2x5Half.metadata[, c("nTReads", "nGenes")], tech = rep("X10X2x5Half", nrow(X10X2x5Half.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(X10Scilife.metadata[, c("nTReads", "nGenes")], tech = rep("X10Scilife", nrow(X10Scilife.metadata))))
hsap.nread.ngene.df <- rbind(hsap.nread.ngene.df, data.frame(CB1.metadata[, c("nTReads", "nGenes")], tech = rep("CB1", nrow(CB1.metadata))))
pdf("all_techs_V2/hsap_nReadsVSnGenes_all_V7.pdf")
#ggplot(hsap.nread.ngene.df, aes(x=nGenes, y=nTReads, group=tech)) + geom_point(aes(color=tech)) + geom_smooth(aes(color=tech))+ xlim(0, 100000)# + geom_line(aes(color=tech))#
ggplot(hsap.nread.ngene.df, aes(x=nTReads, y=nGenes, group=tech)) + geom_smooth(method="lm", formula=y~x, fill="grey", aes(color=tech))+ xlim(0, 100000)# + geom_line(aes(color=tech))#
#ggplot(hsap.nread.ngene.df, aes(x=nTReads, y=nGenes, group=tech)) + geom_smooth(method="lm", formula=y~log(x), fill="grey", aes(color=tech))+ xlim(0, 2000000)# + geom_line(aes(color=tech))#
#ggplot(hsap.nread.ngene.df, aes(x=nTReads, y=nGenes, group=tech)) + geom_smooth(method='nlsLM',formula=y ~ a*x^b, se=FALSE, method.args = list(start = list(a=1,b=1)), aes(color=tech))+ xlim(0, 1500000)# + geom_line(aes(color=tech))#
#ggplot(hsap.nread.ngene.df, aes(x=nTReads, y=nGenes, group=tech)) + geom_smooth(method='nls',formula=y ~ SSasympOff(x, A, lrc, c0), se=FALSE, aes(color=tech))+ xlim(0, 1500000)# + geom_line(aes(color=tech))#
dev.off()