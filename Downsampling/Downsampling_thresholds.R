load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/MARS-Seq/data_seu.obj_res0.5_dim6.RData")
MARSseq.obj <- data
rm(data)
MARSseq.metadata <- MARSseq.obj@meta.data
TSNEPlot(object = MARSseq.obj)
MARSseq.obj@ident <- MARSseq.metadata[names(MARSseq.obj@ident), "nTReads"]

load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/CEL-Seq2/data_seu.obj_res0.5_dim8.RData")
CELseq2.obj <- data
rm(data)
CELseq2.metadata <- CELseq2.obj@meta.data

load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Quartz-seq/data_seu.obj_res0.4_dim10.RData")
QUARTZseq.obj <- data
rm(data)
QUARTZseq.metadata <- QUARTZseq.obj@meta.data

load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Drop-Seq/data_seu.obj_res0.4_dim8.RData")
Dropseq.obj <- data
rm(data)
Dropseq.metadata <- Dropseq.obj@meta.data

load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/SCRB-Seq/data_seu.obj_res0.3_dim10.RData")
SCRBseq.obj <- data
rm(data)
SCRBseq.metadata <- SCRBseq.obj@meta.data

load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/10X_Scilife/data_seu.obj_res0.1_dim8.RData")
X10Scilife.obj <- data
rm(data)
X10Scilife.metadata <- X10Scilife.obj@meta.data

load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/10X-Nuclei/data_seu.obj_res0.4_dim8.RData")
X10Nuclei.obj <- data
rm(data)
X10Nuclei.metadata <- X10Nuclei.obj@meta.data

load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/ICELL8/data_seu.obj_res0.5_dim9.RData")
ICELL8.obj <- data
rm(data)
ICELL8.metadata <- ICELL8.obj@meta.data

# HEKS 
MARSseq.HEK <- names(MARSseq.obj@ident)[which(MARSseq.obj@ident == "HEK cells")]
CELseq2.HEK <- names(CELseq2.obj@ident)[which(CELseq2.obj@ident == "HEK cells")]
QUARTZseq.HEK <- names(QUARTZseq.obj@ident)[which(QUARTZseq.obj@ident == "HEK cells")]
Dropseq.HEK <- names(Dropseq.obj@ident)[which(Dropseq.obj@ident == "HEK cells")]
SCRBseq.HEK <- names(SCRBseq.obj@ident)[which(SCRBseq.obj@ident == "HEK cells")]
X10Scilife.HEK <- names(X10Scilife.obj@ident)[which(X10Scilife.obj@ident == "HEK cells")]
X10Nuclei.HEK <- names(X10Nuclei.obj@ident)[which(X10Nuclei.obj@ident == "HEK cells")]
ICELL8.HEK <- names(ICELL8.obj@ident)[which(ICELL8.obj@ident == "HEK cells")]

MARSseq.numbs.HEK <- paste(length(MARSseq.HEK), "/", dim(MARSseq.metadata)[1], sep = "")
CELseq2.numbs.HEK <- paste(length(CELseq2.HEK), "/", dim(CELseq2.metadata)[1], sep = "")
QUARTZseq.numbs.HEK <- paste(length(QUARTZseq.HEK), "/", dim(QUARTZseq.metadata)[1], sep = "")
Dropseq.numbs.HEK <- paste(length(Dropseq.HEK), "/", dim(Dropseq.metadata)[1], sep = "")
SCRBseq.numbs.HEK <- paste(length(SCRBseq.HEK), "/", dim(SCRBseq.metadata)[1], sep = "")
X10Scilife.numbs.HEK <- paste(length(X10Scilife.HEK), "/", dim(X10Scilife.metadata)[1], sep = "")
X10Nuclei.numbs.HEK <- paste(length(X10Nuclei.HEK), "/", dim(X10Nuclei.metadata)[1], sep = "")
ICELL8.numbs.HEK <- paste(length(ICELL8.HEK), "/", dim(ICELL8.metadata)[1], sep = "")
numb.labels.HEK <- c(MARSseq.numbs.HEK,CELseq2.numbs.HEK,QUARTZseq.numbs.HEK,Dropseq.numbs.HEK,SCRBseq.numbs.HEK, X10Scilife.numbs.HEK,X10Nuclei.numbs.HEK, ICELL8.numbs.HEK)

hsap.total.reads <- data.frame(nreads= MARSseq.metadata[MARSseq.HEK, c("nTReads")], tech = rep("MARSseq", length(MARSseq.HEK)))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= CELseq2.metadata[CELseq2.HEK, c("nTReads")], tech = rep("CELseq2", length(CELseq2.HEK))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= QUARTZseq.metadata[QUARTZseq.HEK, c("nTReads")], tech = rep("QUARTZseq", length(QUARTZseq.HEK))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= Dropseq.metadata[Dropseq.HEK, c("nTReads")], tech = rep("Dropseq", length(Dropseq.HEK))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= SCRBseq.metadata[SCRBseq.HEK, c("nTReads")], tech = rep("SCRBseq", length(SCRBseq.HEK))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= X10Scilife.metadata[X10Scilife.HEK, c("nTReads")], tech = rep("X10Scilife", length(X10Scilife.HEK))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= X10Nuclei.metadata[X10Nuclei.HEK, c("nTReads")], tech = rep("X10Nuclei", length(X10Nuclei.HEK))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= ICELL8.metadata[ICELL8.HEK, c("nTReads")], tech = rep("ICELL8", length(ICELL8.HEK))))

a <- aggregate(nreads ~ tech , hsap.total.reads, function(i) round(mean(i)))
pdf("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/HEKS_DS_thresholds.pdf")
#ggplot(data=hsap.total.reads, aes(x=tech, y=nreads, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
ggplot(data=hsap.total.reads, aes(x=tech, y=nreads, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(10000,20000,50000,100000,500000,1000000))+ facet_grid(. ~ tech, scales = "free") + geom_hline(yintercept=c(20000,50000,150000,1000000), color = "orange", linetype = "dashed") +
  geom_text(data = a, aes(label = numb.labels.HEK), position = position_dodge(width=0.9),vjust = -5)
###    stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.9)
#stat_box_data <- function(y, upper_limit = max(log(hsap.total.reads$nreads)) * 1.15) {
#  return(data.frame(y = 0.95 * upper_limit,label = paste('count =', length(y), '\n')))}
dev.off()



# MonocytesS 
MARSseq.Monocytes <- names(MARSseq.obj@ident)[which(MARSseq.obj@ident == "CD14+ Monocytes")]
CELseq2.Monocytes <- names(CELseq2.obj@ident)[which(CELseq2.obj@ident == "CD14+ Monocytes")]
QUARTZseq.Monocytes <- names(QUARTZseq.obj@ident)[which(QUARTZseq.obj@ident == "CD14+ Monocytes")]
Dropseq.Monocytes <- names(Dropseq.obj@ident)[which(Dropseq.obj@ident == "CD14+ Monocytes")]
SCRBseq.Monocytes <- names(SCRBseq.obj@ident)[which(SCRBseq.obj@ident == "CD14+ Monocytes")]
X10Scilife.Monocytes <- names(X10Scilife.obj@ident)[which(X10Scilife.obj@ident == "CD14+ and FCGR3A+ Monocytes")]
X10Nuclei.Monocytes <- names(X10Nuclei.obj@ident)[which(X10Nuclei.obj@ident == "CD14+ Monocytes")]
ICELL8.Monocytes <- names(ICELL8.obj@ident)[which(ICELL8.obj@ident == "CD14+ and FCGR3A+ Monocytes")]

MARSseq.numbs.Monocytes <- paste(length(MARSseq.Monocytes), "/", dim(MARSseq.metadata)[1], sep = "")
CELseq2.numbs.Monocytes <- paste(length(CELseq2.Monocytes), "/", dim(CELseq2.metadata)[1], sep = "")
QUARTZseq.numbs.Monocytes <- paste(length(QUARTZseq.Monocytes), "/", dim(QUARTZseq.metadata)[1], sep = "")
Dropseq.numbs.Monocytes <- paste(length(Dropseq.Monocytes), "/", dim(Dropseq.metadata)[1], sep = "")
SCRBseq.numbs.Monocytes <- paste(length(SCRBseq.Monocytes), "/", dim(SCRBseq.metadata)[1], sep = "")
X10Scilife.numbs.Monocytes <- paste(length(X10Scilife.Monocytes), "/", dim(X10Scilife.metadata)[1], sep = "")
X10Nuclei.numbs.Monocytes <- paste(length(X10Nuclei.Monocytes), "/", dim(X10Nuclei.metadata)[1], sep = "")
ICELL8.numbs.Monocytes <- paste(length(ICELL8.Monocytes), "/", dim(ICELL8.metadata)[1], sep = "")
numb.labels.Monocytes <- c(MARSseq.numbs.Monocytes,CELseq2.numbs.Monocytes,QUARTZseq.numbs.Monocytes,Dropseq.numbs.Monocytes,SCRBseq.numbs.Monocytes, X10Scilife.numbs.Monocytes,X10Nuclei.numbs.Monocytes, ICELL8.numbs.Monocytes)

hsap.total.reads <- data.frame(nreads= MARSseq.metadata[MARSseq.Monocytes, c("nTReads")], tech = rep("MARSseq", length(MARSseq.Monocytes)))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= CELseq2.metadata[CELseq2.Monocytes, c("nTReads")], tech = rep("CELseq2", length(CELseq2.Monocytes))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= QUARTZseq.metadata[QUARTZseq.Monocytes, c("nTReads")], tech = rep("QUARTZseq", length(QUARTZseq.Monocytes))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= Dropseq.metadata[Dropseq.Monocytes, c("nTReads")], tech = rep("Dropseq", length(Dropseq.Monocytes))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= SCRBseq.metadata[SCRBseq.Monocytes, c("nTReads")], tech = rep("SCRBseq", length(SCRBseq.Monocytes))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= X10Scilife.metadata[X10Scilife.Monocytes, c("nTReads")], tech = rep("X10Scilife", length(X10Scilife.Monocytes))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= X10Nuclei.metadata[X10Nuclei.Monocytes, c("nTReads")], tech = rep("X10Nuclei", length(X10Nuclei.Monocytes))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= ICELL8.metadata[ICELL8.Monocytes, c("nTReads")], tech = rep("ICELL8", length(ICELL8.Monocytes))))

a <- aggregate(nreads ~ tech , hsap.total.reads, function(i) round(mean(i)))
pdf("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/MonocytesS_DS_thresholds.pdf")
#ggplot(data=hsap.total.reads, aes(x=tech, y=nreads, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
ggplot(data=hsap.total.reads, aes(x=tech, y=nreads, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(10000,15000,20000,50000,100000,500000,1000000))+ facet_grid(. ~ tech, scales = "free") + geom_hline(yintercept=c(10000,15000,20000,50000,150000,1000000), color = "orange", linetype = "dashed") +
  geom_text(data = a, aes(label = numb.labels.Monocytes), position = position_dodge(width=0.9),vjust = -5)
###    stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.9)
#stat_box_data <- function(y, upper_limit = max(log(hsap.total.reads$nreads)) * 1.15) {
#  return(data.frame(y = 0.95 * upper_limit,label = paste('count =', length(y), '\n')))}
dev.off()



# BcellsS 
MARSseq.Bcells <- names(MARSseq.obj@ident)[which(MARSseq.obj@ident == "B cells")]
CELseq2.Bcells <- names(CELseq2.obj@ident)[which(CELseq2.obj@ident == "B cells")]
QUARTZseq.Bcells <- names(QUARTZseq.obj@ident)[which(QUARTZseq.obj@ident == "B cells")]
Dropseq.Bcells <- names(Dropseq.obj@ident)[which(Dropseq.obj@ident == "B cells")]
SCRBseq.Bcells <- names(SCRBseq.obj@ident)[which(SCRBseq.obj@ident == "B cells")]
X10Scilife.Bcells <- names(X10Scilife.obj@ident)[which(X10Scilife.obj@ident == "B cells")]
X10Nuclei.Bcells <- names(X10Nuclei.obj@ident)[which(X10Nuclei.obj@ident == "B cells")]
ICELL8.Bcells <- names(ICELL8.obj@ident)[which(ICELL8.obj@ident == "B cells")]

MARSseq.numbs.Bcells <- paste(length(MARSseq.Bcells), "/", dim(MARSseq.metadata)[1], sep = "")
CELseq2.numbs.Bcells <- paste(length(CELseq2.Bcells), "/", dim(CELseq2.metadata)[1], sep = "")
QUARTZseq.numbs.Bcells <- paste(length(QUARTZseq.Bcells), "/", dim(QUARTZseq.metadata)[1], sep = "")
Dropseq.numbs.Bcells <- paste(length(Dropseq.Bcells), "/", dim(Dropseq.metadata)[1], sep = "")
SCRBseq.numbs.Bcells <- paste(length(SCRBseq.Bcells), "/", dim(SCRBseq.metadata)[1], sep = "")
X10Scilife.numbs.Bcells <- paste(length(X10Scilife.Bcells), "/", dim(X10Scilife.metadata)[1], sep = "")
X10Nuclei.numbs.Bcells <- paste(length(X10Nuclei.Bcells), "/", dim(X10Nuclei.metadata)[1], sep = "")
ICELL8.numbs.Bcells <- paste(length(ICELL8.Bcells), "/", dim(ICELL8.metadata)[1], sep = "")
numb.labels.Bcells <- c(MARSseq.numbs.Bcells,CELseq2.numbs.Bcells,QUARTZseq.numbs.Bcells,Dropseq.numbs.Bcells,SCRBseq.numbs.Bcells, X10Scilife.numbs.Bcells,X10Nuclei.numbs.Bcells, ICELL8.numbs.Bcells)

hsap.total.reads <- data.frame(nreads= MARSseq.metadata[MARSseq.Bcells, c("nTReads")], tech = rep("MARSseq", length(MARSseq.Bcells)))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= CELseq2.metadata[CELseq2.Bcells, c("nTReads")], tech = rep("CELseq2", length(CELseq2.Bcells))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= QUARTZseq.metadata[QUARTZseq.Bcells, c("nTReads")], tech = rep("QUARTZseq", length(QUARTZseq.Bcells))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= Dropseq.metadata[Dropseq.Bcells, c("nTReads")], tech = rep("Dropseq", length(Dropseq.Bcells))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= SCRBseq.metadata[SCRBseq.Bcells, c("nTReads")], tech = rep("SCRBseq", length(SCRBseq.Bcells))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= X10Scilife.metadata[X10Scilife.Bcells, c("nTReads")], tech = rep("X10Scilife", length(X10Scilife.Bcells))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= X10Nuclei.metadata[X10Nuclei.Bcells, c("nTReads")], tech = rep("X10Nuclei", length(X10Nuclei.Bcells))))
hsap.total.reads <- rbind(hsap.total.reads, data.frame(nreads= ICELL8.metadata[ICELL8.Bcells, c("nTReads")], tech = rep("ICELL8", length(ICELL8.Bcells))))

a <- aggregate(nreads ~ tech , hsap.total.reads, function(i) round(mean(i)))
pdf("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/BcellsS_DS_thresholds.pdf")
#ggplot(data=hsap.total.reads, aes(x=tech, y=nreads, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + facet_grid(. ~ tech, scales = "free") 
ggplot(data=hsap.total.reads, aes(x=tech, y=nreads, fill=tech)) + geom_boxplot() +theme (axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",panel.spacing.x = unit(0, "null")) +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(10000,15000,20000,50000,100000,500000,1000000))+ facet_grid(. ~ tech, scales = "free") + geom_hline(yintercept=c(10000,15000,20000,50000,150000,1000000), color = "orange", linetype = "dashed") +
  geom_text(data = a, aes(label = numb.labels.Bcells), position = position_dodge(width=0.9),vjust = -5)
###    stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.9)
#stat_box_data <- function(y, upper_limit = max(log(hsap.total.reads$nreads)) * 1.15) {
#  return(data.frame(y = 0.95 * upper_limit,label = paste('count =', length(y), '\n')))}
dev.off()


