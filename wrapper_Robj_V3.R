
"
Wrap all the output from previous steps of the pipline and produce a single SCE object for each sc technology

Usage:
Robject_wrapper.R --hsapExp <DATA_IN> --mmusExp <DATA_IN> --hnReads <DATA_IN> --hnUMI <DATA_IN> --hnGene <DATA_IN> --hnFeatures <DATA_IN> --mnReads <DATA_IN> --mnUMI <DATA_IN> --mnGene <DATA_IN> --mnFeatures <DATA_IN> --species <DATA_IN> --vcfs <DATA_IN> --output_SCEobj <FILE_OUT> [--technology <sc_technology> ]

Options:
-hsapExp --hsapExp <DATA_IN>                list of expression object from Human mapping for all pools
-mmusExp --mmusExp <DATA_IN>                list of expression object from Mouse mapping for all pools
-hnReads --hnReads <DATA_IN>                  list of Number of read files for all pools mapped to human only
-hnUMI --hnUMI <DATA_IN>                      list of Number of UMIs files for all pools mapped to human only
-hnGene --hnGene <DATA_IN>                    list of Number of genes files for all pools mapped to human only
-hnFeatures --hnFeatures <DATA_IN>            list of Number of features files for all pools mapped to human only
-mnReads --mnReads <DATA_IN>                  list of Number of read files for all pools mapped to mouse only
-mnUMI --mnUMI <DATA_IN>                      list of Number of UMIs files for all pools mapped to mouse only
-mnGene --mnGene <DATA_IN>                    list of Number of genes files for all pools mapped to mouse only
-mnFeatures --mnFeatures <DATA_IN>            list of Number of features files for all pools mapped to mouse only
-species --species <DATA_IN>                list of species infor from species deconvolution step from mixed mapping for all pools
-vcfs --vcfs <DATA_IN>                      list of VCF files for all pools
-Output --output_SCEobj <FILE_OUT>          The output location for two R SCE object that has expression data mapped to human ref and mouse Ref reparately
--technology <tech_id>                      The scRNASeq technology
-h --help                                   show this
-v --version                                print version and stop

" -> doc

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(cardelino))

main <- function(hsapExp,mmusExp,hnReads,hnUMI,hnGene,hnFeatures,mnReads,mnUMI,mnGene,mnFeatures, species, 
                 vcfs, output_SCEobj, technology) {
  
  #ng_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/GeneNumber_list.txt", what="", sep="\n")
  #nUMI_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/UMINumber_list.txt", what="", sep="\n")
  #nread_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/ReadNumber_list.txt", what="", sep="\n")
  #species_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/species_list.txt", what="", sep="\n")
  #hsapExp_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/Hsap_expression_list.txt", what="", sep=" ")
  #mmusExp_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/Mmus_expression_list.txt", what="", sep="\n")
  #vcf_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/VCF_list.txt", what="", sep="\n")
  
  #hsapExp = "/Volumes/Ati-Archive/HCA/donor_id/Hsap_expression_list.txt"
  #mmusExp = "/Volumes/Ati-Archive/HCA/donor_id/Mmus_expression_list.txt"
  #nReads = "/Volumes/Ati-Archive/HCA/donor_id/ReadNumber_list.txt"
  #nUMI = "/Volumes/Ati-Archive/HCA/donor_id/UMINumber_list.txt"
  #nGene = "/Volumes/Ati-Archive/HCA/donor_id/GeneNumber_list.txt"
  #nFeatures = "/Volumes/Ati-Archive/HCA/donor_id/Features_list.txt"
  #species = "/Volumes/Ati-Archive/HCA/donor_id/species_list.txt"
  #vcfs = "/Volumes/Ati-Archive/HCA/donor_id/VCF_list.txt"
  
  print("Wrapping starts")
  hsapExp_list <- scan(hsapExp, what="", sep=" ")
  mmusExp_list <- scan(mmusExp, what="", sep=" ")
  hsap_nread_list <- scan(hnReads, what="", sep=" ")
  hsap_nUMI_list <- scan(hnUMI, what="", sep=" ")
  hsap_nGene_list <- scan(hnGene, what="", sep=" ")
  hsap_nFeatures_list <- scan(hnFeatures, what="", sep=" ")
  mmus_nread_list <- scan(mnReads, what="", sep=" ")
  mmus_nUMI_list <- scan(mnUMI, what="", sep=" ")
  mmus_nGene_list <- scan(mnGene, what="", sep=" ")
  mmus_nFeatures_list <- scan(mnFeatures, what="", sep=" ")
  species_list <- scan(species, what="", sep=" ")
  vcf_list <- scan(vcfs, what="", sep=" ")
  number_of_samples <- length(hsapExp_list)
  
  #SCEobj.list <- vector("list", number_of_samples)
  hsap.ExpsMat.list <- list()
  hsap.metadata.list <- list()
  mmus.ExpsMat.list <- list()
  mmus.metadata.list <- list()
  vcfs.list <- list()
  for (sample in 1:number_of_samples){
    sample.hng.file <- read.table(hsap_nGene_list[sample], header = T)
    sample.hnUMI.file <- read.table(hsap_nUMI_list[sample], header = T)
    sample.hnread.file <- read.table(hsap_nread_list[sample], header = T)
    sample.hfeatures.file <- read.table(hsap_nFeatures_list[sample], header = T)
    if (NA %in% sample.hfeatures.file$AssignmentType){
      sample.hfeatures.file$AssignmentType <- factor(sample.hfeatures.file$AssignmentType, levels = levels(addNA(sample.hfeatures.file$AssignmentType)), labels = c(levels(sample.hfeatures.file$AssignmentType), "Multimapping"), exclude = NULL)}

    sample.mng.file <- read.table(mmus_nGene_list[sample], header = T)
    sample.mnUMI.file <- read.table(mmus_nUMI_list[sample], header = T)
    sample.mnread.file <- read.table(mmus_nread_list[sample], header = T)
    sample.mfeatures.file <- read.table(mmus_nFeatures_list[sample], header = T)
    if (NA %in% sample.mfeatures.file$AssignmentType){
      sample.mfeatures.file$AssignmentType <- factor(sample.mfeatures.file$AssignmentType, levels = levels(addNA(sample.mfeatures.file$AssignmentType)), labels = c(levels(sample.mfeatures.file$AssignmentType), "Multimapping"), exclude = NULL)}
    
    sample.species.file <- read.table(species_list[sample])
    sample.hsapExp.obj <- readRDS(hsapExp_list[sample])
    sample.mmusExp.obj <- readRDS(mmusExp_list[sample])
    #sample.vcf <- read.vcfR(vcf_list[sample])
    sample.vcf <- read_vcf(vcf_list[sample], genome = "GRCh38")
    sample.ID.pre <- unlist(strsplit(unlist(strsplit(hsap_nread_list[sample], "/"))[2],"_")) #for local 6
    sample.ID <- paste(sample.ID.pre[-length(sample.ID.pre)], collapse = "_")
    
    #sample.exon.reads <- as.data.frame(sample.features.file[which(sample.features.file$AssignmentType == "exon"), "NumberOfReads"], row.names = create_cell_IDs(sample.features.file[which(sample.features.file$AssignmentType == "exon"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.exon.reads) <- "nExonReads" 
    #sample.intron.reads <- as.data.frame(sample.features.file[which(sample.features.file$AssignmentType == "intron"), "NumberOfReads"], row.names = create_cell_IDs(sample.features.file[which(sample.features.file$AssignmentType == "intron"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.exon.reads) <- "nIntronReads" 
    #sample.intergenic.reads <- as.data.frame(sample.features.file[which(sample.features.file$AssignmentType == "Intergenic"), "NumberOfReads"], row.names = create_cell_IDs(sample.features.file[which(sample.features.file$AssignmentType == "Intergenic"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.exon.reads) <- "nIntergenicReads" 
    #sample.Unmapped.reads <- as.data.frame(sample.features.file[which(sample.features.file$AssignmentType == "Unmapped"), "NumberOfReads"], row.names = create_cell_IDs(sample.features.file[which(sample.features.file$AssignmentType == "Unmapped"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.exon.reads) <- "nUnmappedReads" 
    #sample.Ambiguity.reads <- as.data.frame(sample.features.file[which(sample.features.file$AssignmentType == "Ambiguity"), "NumberOfReads"], row.names = create_cell_IDs(sample.features.file[which(sample.features.file$AssignmentType == "Ambiguity"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.exon.reads) <- "nAmbiguityReads" 
    #sample.ng <- as.data.frame(sample.ng.file[which(sample.ng.file$featureType == "intron.exon"), "Count"], row.names = create_cell_IDs(sample.ng.file[which(sample.ng.file$featureType == "intron.exon"), "SampleID"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.ng) <- "nGenes"
    #sample.nUMI <- as.data.frame(sample.nUMI.file[which(sample.nUMI.file$featureType == "intron.exon"), "Count"], row.names = create_cell_IDs(sample.nUMI.file[which(sample.nUMI.file$featureType == "intron.exon"), "SampleID"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.nUMI) <- "nUMIs"
    #sample.nread <- as.data.frame(sample.nread.file$Total, row.names = create_cell_IDs(sample.nread.file$XC, id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.nread) <- "nReads"
    #sample.species <- as.data.frame(sample.species.file[,2], row.names = create_cell_IDs(sample.species.file, id.type = "standard"), col.names = c("species")); colnames(sample.species) <- "species"
    #sample.species <- as.data.frame(sample.species[rownames(sample.nread),],row.names = rownames(sample.nread)); colnames(sample.species) <- "Species"
    
    sample.hsap.exon.reads <- data.frame(nExonReads=sample.hfeatures.file[which(sample.hfeatures.file$AssignmentType == "exon"), "NumberOfReads"], cellID=create_cell_IDs(sample.hfeatures.file[which(sample.hfeatures.file$AssignmentType == "exon"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID))#; colnames(sample.exon.reads) <- "nExonReads" 
    sample.hsap.intron.reads <- data.frame(nIntronReads=sample.hfeatures.file[which(sample.hfeatures.file$AssignmentType == "intron"), "NumberOfReads"], cellID = create_cell_IDs(sample.hfeatures.file[which(sample.hfeatures.file$AssignmentType == "intron"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)) 
    sample.hsap.intergenic.reads <- data.frame(nIntergenicReads=sample.hfeatures.file[which(sample.hfeatures.file$AssignmentType == "Intergenic"), "NumberOfReads"], cellID = create_cell_IDs(sample.hfeatures.file[which(sample.hfeatures.file$AssignmentType == "Intergenic"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)) 
    sample.hsap.Unmapped.reads <- data.frame(nUnmappedReads=sample.hfeatures.file[which(sample.hfeatures.file$AssignmentType == "Unmapped"), "NumberOfReads"], cellID = create_cell_IDs(sample.hfeatures.file[which(sample.hfeatures.file$AssignmentType == "Unmapped"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID))
    sample.hsap.Ambiguity.reads <- data.frame(nAmbiguityReads=sample.hfeatures.file[which(sample.hfeatures.file$AssignmentType == "Ambiguity"), "NumberOfReads"], cellID = create_cell_IDs(sample.hfeatures.file[which(sample.hfeatures.file$AssignmentType == "Ambiguity"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)) 
    sample.hsap.Multimap.reads <- data.frame(nMultimapReads=sample.hfeatures.file[which(sample.hfeatures.file$AssignmentType == "Multimapping"), "NumberOfReads"], cellID = create_cell_IDs(sample.hfeatures.file[which(sample.hfeatures.file$AssignmentType == "Multimapping"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)) 
    sample.hsap.ng <- data.frame(nGenes=sample.hng.file[which(sample.hng.file$featureType == "intron.exon"), "Count"], cellID = create_cell_IDs(sample.hng.file[which(sample.hng.file$featureType == "intron.exon"), "SampleID"], id.type = "cell_Barcode",tech = technology, lib = sample.ID))
    sample.hsap.nUMI <- data.frame(nUMIs=sample.hnUMI.file[which(sample.hnUMI.file$featureType == "intron.exon"), "Count"], cellID = create_cell_IDs(sample.hnUMI.file[which(sample.hnUMI.file$featureType == "intron.exon"), "SampleID"], id.type = "cell_Barcode",tech = technology, lib = sample.ID))
    sample.hsap.nread <- data.frame(nTReads=sample.hnread.file$Total, cellID = create_cell_IDs(sample.hnread.file$XC, id.type = "cell_Barcode",tech = technology, lib = sample.ID))

    sample.mmus.exon.reads <- data.frame(nExonReads=sample.mfeatures.file[which(sample.mfeatures.file$AssignmentType == "exon"), "NumberOfReads"], cellID=create_cell_IDs(sample.mfeatures.file[which(sample.mfeatures.file$AssignmentType == "exon"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID))#; colnames(sample.exon.reads) <- "nExonReads" 
    sample.mmus.intron.reads <- data.frame(nIntronReads=sample.mfeatures.file[which(sample.mfeatures.file$AssignmentType == "intron"), "NumberOfReads"], cellID = create_cell_IDs(sample.mfeatures.file[which(sample.mfeatures.file$AssignmentType == "intron"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)) 
    sample.mmus.intergenic.reads <- data.frame(nIntergenicReads=sample.mfeatures.file[which(sample.mfeatures.file$AssignmentType == "Intergenic"), "NumberOfReads"], cellID = create_cell_IDs(sample.mfeatures.file[which(sample.mfeatures.file$AssignmentType == "Intergenic"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)) 
    sample.mmus.Unmapped.reads <- data.frame(nUnmappedReads=sample.mfeatures.file[which(sample.mfeatures.file$AssignmentType == "Unmapped"), "NumberOfReads"], cellID = create_cell_IDs(sample.mfeatures.file[which(sample.mfeatures.file$AssignmentType == "Unmapped"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID))
    sample.mmus.Ambiguity.reads <- data.frame(nAmbiguityReads=sample.mfeatures.file[which(sample.mfeatures.file$AssignmentType == "Ambiguity"), "NumberOfReads"], cellID = create_cell_IDs(sample.mfeatures.file[which(sample.mfeatures.file$AssignmentType == "Ambiguity"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)) 
    sample.mmus.Multimap.reads <- data.frame(nMultimapReads=sample.mfeatures.file[which(sample.mfeatures.file$AssignmentType == "Multimapping"), "NumberOfReads"], cellID = create_cell_IDs(sample.mfeatures.file[which(sample.mfeatures.file$AssignmentType == "Multimapping"), "XC"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)) 
    sample.mmus.ng <- data.frame(nGenes=sample.mng.file[which(sample.mng.file$featureType == "intron.exon"), "Count"], cellID = create_cell_IDs(sample.mng.file[which(sample.mng.file$featureType == "intron.exon"), "SampleID"], id.type = "cell_Barcode",tech = technology, lib = sample.ID))
    sample.mmus.nUMI <- data.frame(nUMIs=sample.mnUMI.file[which(sample.mnUMI.file$featureType == "intron.exon"), "Count"], cellID = create_cell_IDs(sample.mnUMI.file[which(sample.mnUMI.file$featureType == "intron.exon"), "SampleID"], id.type = "cell_Barcode",tech = technology, lib = sample.ID))
    sample.mmus.nread <- data.frame(nTReads=sample.mnread.file$Total, cellID = create_cell_IDs(sample.mnread.file$XC, id.type = "cell_Barcode",tech = technology, lib = sample.ID))
    
    sample.species <- data.frame(Species=sample.species.file[,2], cellID = create_cell_IDs(sample.species.file, id.type = "standard",tech = technology, lib = sample.ID))
    sample.species <- sample.species[rownames(sample.nread),]
    
        
    ids <- donor_id(sample.vcf, n_donor= 5, n_vars_threshold = 5)
    #table(ids$assigned$donor_id)
    #head(ids$assigned)
    #donor.assigned <- ids$assigned$donor_id
    #names(donor.assigned) <- create_cell_IDs(ids$assigned$cell, id.type = "vcf_id", tech = technology, lib = sample.ID)
    nvars.per.cell <- ids$assigned$n_vars
    names(nvars.per.cell) <- create_cell_IDs(ids$assigned$cell, id.type = "vcf_id", tech = technology, lib = sample.ID)
    s.ProbMat <- ids$prob
    rownames(s.ProbMat) <- create_cell_IDs(rownames(s.ProbMat), id.type = "vcf_id", tech = technology, lib = sample.ID)
    donor.assigned= rep("Not_human", nrow(sample.nread))
    names(donor.assigned) <- sample.nread$cellID
    for (cellBC in rownames(s.ProbMat)){
      donorID = "unassigned"
      for (donor in colnames(s.ProbMat)){
        if (s.ProbMat[cellBC, donor] > 0.5){
          donorID = donor}
        donor.assigned[cellBC] <- donorID}}
    sample.donor = data.frame(Donor=donor.assigned, cellID= sample.nread$cellID)
    sample.library=data.frame(Library=rep(sample.ID, nrow(sample.nread)), cellID=sample.nread$cellID)
    sample.nVars= data.frame(nVars=rep(NA, nrow(sample.nread)), cellID=sample.nread$cellID, row.names = sample.nread$cellID)
    sample.nVars[names(nvars.per.cell), "nVars"] <- nvars.per.cell

    #prob_mat = s.ProbMat
    #hc <- hclust(dist(prob_mat))
    #nba.m <- as_data_frame(prob_mat[hc$order,]) %>%
    #  dplyr::mutate(cell = rownames(prob_mat[hc$order,])) %>%
    #  gather(key = "clone", value = "prob", -cell)
    
    #nba.m <- dplyr::mutate(nba.m, cell = factor(
    #  cell, levels = rownames(prob_mat[hc$order,])))
    
    #ggplot(nba.m, aes(clone, cell, fill = prob)) + 
    #geom_tile(show.legend = TRUE) +
    #scale_fill_gradient(low = "white", high = "firebrick4",
    #                      name = "posterior probability of assignment") +
    #ylab(paste("Single cells")) + 
    #cardelino:::heatmap.theme(size = 16) + #cardelino:::pub.theme() +
    #theme(axis.title.y = element_text(size = 20), legend.position = "bottom",
    #      legend.text = element_text(size = 12), legend.key.size = unit(0.05, "npc"))
    
#    sample.metadata <- as.data.frame(cbind(TotalnReads=sample.nread, nExonReads=sample.exon.reads, nIntronReads=sample.intron.reads, 
#                                           nIntergenicReads=sample.intergenic.reads, nUnmappedReads=sample.Unmapped.reads, nAmbiguityReads=sample.Ambiguity.reads,
#                                           nUMI=sample.nUMI, nGenes=sample.ng, Species=sample.species, library=rep(sample.ID, nrow(sample.nread)), donor_assigned=donor.assigned[rownames(sample.nread)], nVar=nvars.per.cell[rownames(sample.nread)]))
    
    sample.hsap.metadata <- join_all(list(sample.hsap.nread, sample.hsap.exon.reads,sample.hsap.intron.reads, 
                                           sample.hsap.intergenic.reads, sample.hsap.Unmapped.reads, hsap.sample.Ambiguity.reads, sample.hsap.Multimap.reads,
                                           sample.hsap.nUMI, sample.hsap.ng, sample.species, sample.library, sample.donor, nVar=sample.nVars), by = "cellID", type = 'full')
    rownames(sample.hsap.metadata) <- sample.hsap.metadata$cellID
    sample.hsap.metadata <- sample.hsap.metadata[,-2]
    hsap.metadata.list[[sample.ID]] <- sample.hsap.metadata
    
    sample.mmus.metadata <- join_all(list(sample.mmus.nread, sample.mmus.exon.reads,sample.mmus.intron.reads, 
                                          sample.mmus.intergenic.reads, sample.mmus.Unmapped.reads, mmus.sample.Ambiguity.reads, sample.mmus.Multimap.reads,
                                          sample.mmus.nUMI, sample.mmus.ng, sample.species, sample.library, sample.donor, nVar=sample.nVars), by = "cellID", type = 'full')
    rownames(sample.mmus.metadata) <- sample.mmus.metadata$cellID
    sample.mmus.metadata <- sample.mmus.metadata[,-2]
    mmus.metadata.list[[sample.ID]] <- sample.hsap.metadata
    
    vcfs.list[[sample.ID]] <- sample.vcf
    
    sample.hsapExp <- as.data.frame(as.matrix(sample.hsapExp.obj$intron.exon$umicounts))
    colnames(sample.hsapExp) <- create_cell_IDs(colnames(sample.hsapExp), id.type = "cell_Barcode",tech = technology, lib = sample.ID)
    sample.hsapExp$rn <- rownames(sample.hsapExp)
    hsap.ExpsMat.list[[sample.ID]] <- sample.hsapExp
    #sample.hsapMetadata <- as.data.frame(cbind(nGenes=sample.nread, nUMIs=sample.nUMI, nReads=sample.nread, species=sample.species, library=rep(sample.ID, nrow(sample.nread))))

    sample.mmusExp <- as.data.frame(as.matrix(sample.mmusExp.obj$intron.exon$umicounts))
    colnames(sample.mmusExp) <- create_cell_IDs(colnames(sample.mmusExp), id.type = "cell_Barcode",tech = technology, lib = sample.ID)
    sample.mmusExp$rn <- rownames(sample.mmusExp)
    mmus.ExpsMat.list[[sample.ID]] <- sample.mmusExp
    #sample.hsapMetadata <- as.data.frame(cbind(nGenes=sample.nread, nUMIs=sample.nUMI, nReads=sample.nread, species=sample.species, library=rep(sample.ID, nrow(sample.nread))))
    
    #sample.sce <- SingleCellExperiment(assays = list(counts = as.matrix(sample.mixedExp)), colData = sample.metadata)
    #SCEobj.list[[sample]] <- sample.sce
    #names(SCEobj.list)[sample] <- sample.ID
    
  }
  print("loop is done")
  all.hsapExp <- join_all(hsap.ExpsMat.list, by = "rn", type = 'full')
  rownames(all.hsapExp) <- all.hsapExp$rn
  all.hsapExp <- all.hsapExp[,!(names(all.hsapExp) %in% c("rn"))]
  all.hsapExp[is.na(all.hsapExp)] <- 0
  
  all.mmusExp <- join_all(mmus.ExpsMat.list, by = "rn", type = 'full')
  rownames(all.mmusExp) <- all.mmusExp$rn
  all.mmusExp <- all.mmusExp[,!(names(all.mmusExp) %in% c("rn"))]
  all.mmusExp[is.na(all.mmusExp)] <- 0
  
  all.hsap.metadata <- do.call("rbind", hsap.metadata.list)
  rownames(all.hsap.metadata) <- lapply(rownames(all.hsap.metadata), function (x) unlist(strsplit(x, "[.]"))[2])
  
  all.mmus.metadata <- do.call("rbind", mmus.metadata.list)
  rownames(all.mmus.metadata) <- lapply(rownames(all.mmus.metadata), function (x) unlist(strsplit(x, "[.]"))[2])
  
  print("Making the SingleCellExperiment Objects")
  full.SCE.hsap <- SingleCellExperiment(assays = list(counts = as.matrix(all.hsapExp)), colData = all.hsap.metadata[colnames(all.hsapExp),])
  full.SCE.hsap <- calculateQCMetrics(full.SCE.hsap)
  libsize.drop.hsap <- isOutlier(full.SCE.hsap$total_counts, nmads=3, type="lower", log=TRUE)
  feature.drop.hsap <- isOutlier(full.SCE.hsap$total_features, nmads=3, type="lower", log=TRUE)
  full.SCE.filterred.hsap <- full.SCE.hsap[,!(libsize.drop.hsap | feature.drop.hsap)]
  data.frame(ByLibSize=sum(libsize.drop.hsap), ByFeature=sum(feature.drop.hsap), Remaining=ncol(full.SCE.filterred.hsap))
  
  full.SCE.mmus <- SingleCellExperiment(assays = list(counts = as.matrix(all.mmusExp)), colData = all.mmus.metadata[colnames(all.hsapExp),])
  full.SCE.mmus <- calculateQCMetrics(full.SCE.mmus)
  libsize.drop.mmus <- isOutlier(full.SCE.mmus$total_counts, nmads=3, type="lower", log=TRUE)
  feature.drop.mmus <- isOutlier(full.SCE.mmus$total_features, nmads=3, type="lower", log=TRUE)
  full.SCE.filterred.mmus <- full.SCE.mmus[,!(libsize.drop.mmus | feature.drop.mmus)]
  data.frame(ByLibSize=sum(libsize.drop.mmus), ByFeature=sum(feature.drop.mmus), Remaining=ncol(full.SCE.filterred.mmus))
  
  #save(full.SCE.filterred.hsap, file = paste(output_SCEobj,"/hsap.full.SCE.Robj", sep = ""))
  #save(full.SCE.filterred.mmus, file = paste(output_SCEobj,"/mmus.full.SCE.Robj", sep = ""))
  save(full.SCE.hsap, file = paste(output_SCEobj,"/hsap.full.SCE.Robj", sep = ""))
  save(full.SCE.mmus, file = paste(output_SCEobj,"/mmus.full.SCE.Robj", sep = ""))
  
  return("Done")
  
}

#counts(all.sce)[1:5,1:5]
#colData(all.sce)
#rowData(all.sce)

create_cell_IDs <- function(cell.IDs, id.type = "standard", tech = technology, lib = sample.ID){
  if (id.type == "standard"){
    pre.ids <- as.character(cell.IDs[,1])
    IDs.parts <- lapply(pre.ids, function(x) unlist(strsplit(x, split = ".", fixed = T)))
    new.IDs <- lapply(IDs.parts, function(x) paste(x[1], x[3], x[4], sep = "_"))
    return(as.character(new.IDs))}
  if (id.type == "cell_Barcode"){
    pre.ids <- cell.IDs
    new.IDs <- lapply(pre.ids, function(x) paste(tech, lib, x, sep = "_"))
    return(as.character(new.IDs))}
  if (id.type == "vcf_id"){
    pre.ids <- cell.IDs
    IDs.parts1 <- lapply(pre.ids, function(x) unlist(strsplit(x, split = "/", fixed = T))[[4]])
    IDs.parts2 <- lapply(IDs.parts1, function(x) unlist(strsplit(x, split = ".", fixed = T)))
    new.IDs <- lapply(IDs.parts2, function(x) paste(x[1], x[3], x[4], sep = "_"))
  }
}


## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")

#message("working directory: ", getwd(), "\n")
#message("input file: ", opt$input_file, "\n")
#message("output file: ", opt$output_file, "\n")


## Run main function
main(opt$hsapExp, opt$mmusExp, opt$hnReads, opt$hnUMI, opt$hnGene, opt$hnFeatures, 
     opt$mnReads, opt$mnUMI, opt$mnGene, opt$mnFeatures, 
     opt$species, opt$vcfs, opt$output_SCEobj, opt$technology)

