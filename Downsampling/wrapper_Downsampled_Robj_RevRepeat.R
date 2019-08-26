"
Wrap the Downsampleded expression matrices

Usage:
wrapper_Downsampled_Robj.R --hsapExp <DATA_IN> --RepNumber <Repeat_number> --output_SCEobj <FILE_OUT> [--technology <sc_technology> ]

Options:
-hsapExp --hsapExp <DATA_IN>                list of Downsampled expression object from Human mapping for all pools
-Output --output_SCEobj <FILE_OUT>          The output location for R SCE object that has expression data mapped to human ref and mouse Ref reparately
--technology <tech_id>                      The scRNASeq technology
--RepNumber <Repeat_number>                 The number of the repetition
-h --help                                   show this
-v --version                                print version and stop

" -> doc

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(data.table))

main <- function(hsapExp, output_SCEobj, RepNumber, technology) {
  
  print("Wrapping starts")
  hsapExp_list <- scan(hsapExp, what="", sep=" ")
  number_of_samples <- length(hsapExp_list)
  init.hsapExp.obj <- readRDS(hsapExp_list[1])
  #DS.ranges <- names(init.hsapExp.obj$intron.exon$downsampled) #change ds.data.mixed to init.hsapExp.obj
  DS.ranges <- names(init.hsapExp.obj$umicount$inex$downsampling) #change ds.data.mixed to init.hsapExp.obj
  
  hsap.DS.ExpsMat.UMI.all <- vector("list", length(DS.ranges))
  names(hsap.DS.ExpsMat.UMI.all) <- DS.ranges
  
  hsap.DS.ExpsMat.Reads.all <- vector("list", length(DS.ranges))
  names(hsap.DS.ExpsMat.Reads.all) <- DS.ranges
  
  output.readcount.umicount.joint.mats <- list()
  
  for (sample in 1:number_of_samples){
    print ("in the loop")
    sample.hsapExp.obj <- readRDS(hsapExp_list[sample])
    sample.ID.pre <- unlist(strsplit(unlist(strsplit(hsapExp_list[sample], "/"))[2],"_")) #for local 6
    print(sample.ID.pre)
    sample.ID <- paste(sample.ID.pre[-length(sample.ID.pre)], collapse = ".")
    print(sample.ID)
    
    # for UMI counts
    for (range in 1:length(DS.ranges)){
      UMI.mat <- as.data.frame(as.matrix(sample.hsapExp.obj$umicount$inex$downsampling[[range]]))
      colnames(UMI.mat) <- create_cell_IDs(colnames(UMI.mat), id.type = "cell_Barcode",tech = technology, lib = sample.ID)
      UMI.mat$rn <- rownames(UMI.mat)
      hsap.DS.ExpsMat.UMI.all[[range]][[sample.ID]] <- UMI.mat #change ds.data.mixed to init.hsapExp.obj
      #hsap.DS.ExpsMat.UMI.all[[range]][[sample.ID]] <- c(hsap.DS.ExpsMat.UMI.all[[range]][[sample.ID]], list(mat))#change ds.data.mixed to init.hsapExp.obj
      #hsap.DS.ExpsMat.UMI.all[[range]] <- list(hsap.DS.ExpsMat.UMI.all[[range]], mat)#change ds.data.mixed to init.hsapExp.obj
      
    }
    
    #for Read counts
    for (range in 1:length(DS.ranges)){
      Reads.mat <- as.data.frame(as.matrix(sample.hsapExp.obj$readcount$inex$downsampling[[range]]))
      colnames(Reads.mat) <- create_cell_IDs(colnames(Reads.mat), id.type = "cell_Barcode",tech = technology, lib = sample.ID)
      Reads.mat$rn <- rownames(Reads.mat)
      hsap.DS.ExpsMat.Reads.all[[range]][[sample.ID]] <- Reads.mat #change ds.data.mixed to init.hsapExp.obj
      #hsap.DS.ExpsMat.UMI.all[[range]][[sample.ID]] <- c(hsap.DS.ExpsMat.UMI.all[[range]][[sample.ID]], list(mat))#change ds.data.mixed to init.hsapExp.obj
      #hsap.DS.ExpsMat.UMI.all[[range]] <- list(hsap.DS.ExpsMat.UMI.all[[range]], mat)#change ds.data.mixed to init.hsapExp.obj
      
    }
    
  }
  
  #UMI final
  final.sample.merged.mat.UMI <- list()
  final.sample.merged.mat.Reads <- list()
  for (i in 1:length(DS.ranges)){
    ds.name = DS.ranges[i]
    
    #join for UMI mats
    joint.mat.UMI <- join_all(hsap.DS.ExpsMat.UMI.all[[i]], by = "rn", type = 'full')
    rownames(joint.mat.UMI) <- joint.mat.UMI$rn
    joint.mat.UMI <- joint.mat.UMI[,!(names(joint.mat.UMI) %in% c("rn"))]
    joint.mat.UMI[is.na(joint.mat.UMI)] <- 0
    final.sample.merged.mat.UMI[[ds.name]] <- joint.mat.UMI
    
    #join for UMI mats
    joint.mat.Reads <- join_all(hsap.DS.ExpsMat.Reads.all[[i]], by = "rn", type = 'full')
    rownames(joint.mat.Reads) <- joint.mat.Reads$rn
    joint.mat.Reads <- joint.mat.Reads[,!(names(joint.mat.Reads) %in% c("rn"))]
    joint.mat.Reads[is.na(joint.mat.Reads)] <- 0
    final.sample.merged.mat.Reads[[ds.name]] <- joint.mat.Reads
  }
  
  output.readcount.umicount.joint.mats[["UMI"]] <- final.sample.merged.mat.UMI
  output.readcount.umicount.joint.mats[["Reads"]] <- final.sample.merged.mat.Reads
  #save(hsap.DS.ExpsMat.UMI.all, file = paste(output_SCEobj,"/", technology,".hsap.full.SCE.Robj", sep = ""))
  save(output.readcount.umicount.joint.mats, file = paste(output_SCEobj,"/", technology,".hsap.full.SCE.jointDSmat_",RepNumber,".Robj", sep = ""))
  #save(full.SCE.mmus, file = paste(output_SCEobj,"/", technology,".mmus.full.SCE.Robj", sep = ""))
  
  return("Done")
  
}


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
main(opt$hsapExp, opt$output_SCEobj,opt$RepNumber, opt$technology)