"
Wrap the Downsampleded expression matrices

Usage:
wrapper_Downsampled_Robj.R --hsapExp <DATA_IN> --output_SCEobj <FILE_OUT> [--technology <sc_technology> ]

Options:
-hsapExp --hsapExp <DATA_IN>                list of Downsampled expression object from Human mapping for all pools
-Output --output_SCEobj <FILE_OUT>          The output location for R SCE object that has expression data mapped to human ref and mouse Ref reparately
--technology <tech_id>                      The scRNASeq technology
-h --help                                   show this
-v --version                                print version and stop

" -> doc

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(data.table))

main <- function(hsapExp, output_SCEobj, technology) {
  
  print("Wrapping starts")
  hsapExp_list <- scan(hsapExp, what="", sep=" ")
  number_of_samples <- length(hsapExp_list)
  init.hsapExp.obj <- readRDS(hsapExp_list[1])
  DS.ranges <- names(init.hsapExp.obj$intron.exon$downsampled) #change ds.data.mixed to init.hsapExp.obj
  
  hsap.DS.ExpsMat.all <- vector("list", length(DS.ranges))
  names(hsap.DS.ExpsMat.all) <- DS.ranges
  for (sample in 1:number_of_samples){
    print ("in the loop")
    sample.hsapExp.obj <- readRDS(hsapExp_list[sample])
    sample.ID.pre <- unlist(strsplit(unlist(strsplit(hsapExp_list[sample], "/"))[2],"_")) #for local 6
    print(sample.ID.pre)
    sample.ID <- paste(sample.ID.pre[-length(sample.ID.pre)], collapse = ".")
    print(sample.ID)
    
    for (range in 1:length(DS.ranges)){
      mat <- as.data.frame(as.matrix(sample.hsapExp.obj$intron.exon$downsampled[[range]][[1]]))
      colnames(mat) <- create_cell_IDs(colnames(mat), id.type = "cell_Barcode",tech = technology, lib = sample.ID)
      mat$rn <- rownames(mat)
      hsap.DS.ExpsMat.all[[range]][[sample.ID]] <- c(hsap.DS.ExpsMat.all[[range]][[sample.ID]], mat)#change ds.data.mixed to init.hsapExp.obj
    }
  }
  save(hsap.DS.ExpsMat.all, file = paste(output_SCEobj,"/", technology,".hsap.full.SCE.Robj", sep = ""))
  #save(full.SCE.hsap, file = paste(output_SCEobj,"/", technology,".hsap.full.SCE.Robj", sep = ""))
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
main(opt$hsapExp, opt$output_SCEobj, opt$technology)


#hsap.DS.ExpsMat.all
#number_of_samples = names(hsap.DS.ExpsMat.all[[1]])
#final.list = list()
#for (DS.range in 1:length(hsap.DS.ExpsMat.all)){
#  range.number = names(hsap.DS.ExpsMat.all[[DS.range]])
#  range.all.samples = join_all(hsap.DS.ExpsMat.all[[DS.range]])
  
#}
