
"
Wrap all the output from previous steps of the pipline and produce a single SCE object for each sc technology

Usage:
Robject_wrapper.R --hsapExp <DATA_IN> --mmusExp <DATA_IN> --nReads <DATA_IN> --species <DATA_IN> --vcfs <DATA_IN> --output_SCEobj <FILE_OUT> [--technology <sc_technology> ]

Options:
-hsapExp --hsapExp <DATA_IN>                list of expression object from Human mapping for all pools
-mmusExp --mmusExp <DATA_IN>                list of expression object from Mouse mapping for all pools
-nReads --nReads <DATA_IN>                  list of Number of read files for all pools
-species --species <DATA_IN>                list of species infor from species deconvolution step from mixed mapping for all pools
-vcfs --vcfs <DATA_IN>                      list of VCF files for all pools
-Output --output_SCEobj <FILE_OUT>          The output location for A R SCE object that has expression data mapped to human ref
--technology <tech_id>                      The scRNASeq technology
-h --help                                   show this
-v --version                                print version and stop

" -> doc

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(cardelino))

main <- function(hsapExp,mmusExp,nReads, species, 
                 vcfs, output_SCEobj, technology) {
  
  #ng_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/GeneNumber_list.txt", what="", sep="\n")
  #nUMI_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/UMINumber_list.txt", what="", sep="\n")
  #nread_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/ReadNumber_list.txt", what="", sep="\n")
  #species_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/species_list.txt", what="", sep="\n")
  #hsapExp_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/Hsap_expression_list.txt", what="", sep=" ")
  #mmusExp_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/Mmus_expression_list.txt", what="", sep="\n")
  #vcf_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/VCF_list.txt", what="", sep="\n")
  
  print("Wrapping starts")
  hsapExp_list <- scan(hsapExp, what="", sep=" ")
  mmusExp_list <- scan(mmusExp, what="", sep=" ")
  nread_list <- scan(nReads, what="", sep=" ")
  species_list <- scan(species, what="", sep=" ")
  vcf_list <- scan(vcfs, what="", sep=" ")
  number_of_samples <- length(hsapExp_list)
  
  #SCEobj.list <- vector("list", number_of_samples)
  hsap.ExpsMat.list <- list()
  mmus.ExpsMat.list <- list()
  metadata.list <- list()
  vcfs.list <- list()
  for (sample in 1:number_of_samples){
    #sample.ng.file <- read.table(ng_list[sample], header = T)
    #sample.nUMI.file <- read.table(nUMI_list[sample], header = T)
    sample.nread.file <- read.table(nread_list[sample], header = T)
    sample.species.file <- read.table(species_list[sample])
    sample.hsapExp.obj <- readRDS(hsapExp_list[sample])
    sample.mmusExp.obj <- readRDS(mmusExp_list[sample])
    #sample.vcf <- read.vcfR(vcf_list[sample])
    sample.vcf <- read_vcf(vcf_list[sample], genome = "GRCh38")
    sample.ID.pre <- unlist(strsplit(unlist(strsplit(nread_list[sample], "/"))[2],"_")) #for local 6
    sample.ID <- paste(sample.ID.pre[-length(sample.ID.pre)], collapse = "_")
    
    #sample.ng <- as.data.frame(sample.ng.file[which(sample.ng.file$featureType == "exons"), "Count"], row.names = create_cell_IDs(sample.ng.file[which(sample.ng.file$featureType == "exons"), "SampleID"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.ng) <- "nGenes"
    #sample.nUMI <- as.data.frame(sample.nUMI.file[which(sample.nUMI.file$featureType == "exons"), "Count"], row.names = create_cell_IDs(sample.nUMI.file[which(sample.nUMI.file$featureType == "exons"), "SampleID"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.nUMI) <- "nUMIs"
    sample.nread <- as.data.frame(sample.nread.file$Total, row.names = create_cell_IDs(sample.nread.file$XC, id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.nread) <- "nReads"
    sample.species <- as.data.frame(sample.species.file[,2], row.names = create_cell_IDs(sample.species.file, id.type = "standard"), col.names = c("species")); colnames(sample.species) <- "species"
    sample.species <- as.data.frame(sample.species[rownames(sample.nread),],row.names = rownames(sample.nread)); colnames(sample.species) <- "Species"
    
    ids <- donor_id(sample.vcf, n_donor= 4, n_vars_threshold = 5)
    #table(ids$assigned$donor_id)
    #head(ids$assigned)
    #donor.assigned <- ids$assigned$donor_id
    #names(donor.assigned) <- create_cell_IDs(ids$assigned$cell, id.type = "vcf_id", tech = technology, lib = sample.ID)
    nvars.per.cell <- ids$assigned$n_vars
    names(nvars.per.cell) <- create_cell_IDs(ids$assigned$cell, id.type = "vcf_id", tech = technology, lib = sample.ID)
    s.ProbMat <- ids$prob
    rownames(s.ProbMat) <- create_cell_IDs(rownames(s.ProbMat), id.type = "vcf_id", tech = technology, lib = sample.ID)
    donor.assigned= rep("Not_human", nrow(sample.nread))
    names(donor.assigned) <- rownames(sample.nread)
    for (cellBC in rownames(s.ProbMat)){
      donorID = "unassigned"
      for (donor in colnames(s.ProbMat)){
        if (s.ProbMat[cellBC, donor] > 0.5){
          donorID = donor}
        donor.assigned[cellBC] <- donorID}}

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
    
    sample.metadata <- as.data.frame(cbind(nReads=sample.nread, species=sample.species, library=rep(sample.ID, nrow(sample.nread)), donor_assigned=donor.assigned[rownames(sample.nread)], nVar=nvars.per.cell[rownames(sample.nread)]))
    metadata.list[[sample.ID]] <- sample.metadata
    vcfs.list[[sample.ID]] <- sample.vcf
    
    sample.hsapExp <- as.data.frame(as.matrix(sample.hsapExp.obj$intron.exon$umicounts))
    colnames(sample.hsapExp) <- create_cell_IDs(colnames(sample.hsapExp), id.type = "cell_Barcode",tech = technology, lib = sample.ID)
    sample.hsapExp$rn <- rownames(sample.hsapExp)
    hsap.ExpsMat.list[[sample.ID]] <- sample.hsapExp
    #sample.hsapMetadata <- as.data.frame(cbind(nGenes=sample.nread, nUMIs=sample.nUMI, nReads=sample.nread, species=sample.species, library=rep(sample.ID, nrow(sample.nread))))
    if (sample.ID == "5570AA"){
      print("here some statistics about 5570AA:")
      print(dim(sample.nread))
      print(dim(sample.species))
      print(dim(sample.hsapExp.obj$intron.exon$umicounts))
    }
    
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
  print("we reached line 128")
  rownames(all.hsapExp) <- all.hsapExp$rn
  print("we reached line 130")
  all.hsapExp <- all.hsapExp[,!(names(all.hsapExp) %in% c("rn"))]
  all.hsapExp[is.na(all.hsapExp)] <- 0
  
  all.mmusExp <- join_all(mmus.ExpsMat.list, by = "rn", type = 'full')
  rownames(all.mmusExp) <- all.mmusExp$rn
  all.mmusExp <- all.mmusExp[,!(names(all.mmusExp) %in% c("rn"))]
  all.mmusExp[is.na(all.mmusExp)] <- 0
  
  print("we reached line 139")
  all.metadata <- do.call("rbind", metadata.list)
  rownames(all.metadata) <- lapply(rownames(all.metadata), function (x) unlist(strsplit(x, "[.]"))[2])
  
  print("we reached line 143")
  
  print(dim(all.metadata))
  print(dim(all.hsapExp))
  #print(rownames(all.metadata))
  #print(colnames(all.hsapExp))
  extra.col = setdiff(rownames(all.metadata), colnames(all.hsapExp))
  print(extra.col)
  full.SCE.hsap <- SingleCellExperiment(assays = list(counts = as.matrix(all.hsapExp)), colData = all.metadata[,])
  full.SCE.hsap <- calculateQCMetrics(full.SCE.hsap)
  libsize.drop.hsap <- isOutlier(full.SCE.hsap$total_counts, nmads=3, type="lower", log=TRUE)
  feature.drop.hsap <- isOutlier(full.SCE.hsap$total_features, nmads=3, type="lower", log=TRUE)
  full.SCE.filterred.hsap <- full.SCE.hsap[,!(libsize.drop.hsap | feature.drop.hsap)]
  data.frame(ByLibSize=sum(libsize.drop.hsap), ByFeature=sum(feature.drop.hsap), Remaining=ncol(full.SCE.filterred.hsap))
  
  full.SCE.mmus <- SingleCellExperiment(assays = list(counts = as.matrix(all.mmusExp)), colData = all.metadata)
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
main(opt$hsapExp, opt$mmusExp, opt$nReads, opt$species, 
     opt$vcfs, opt$output_SCEobj, opt$technology)
