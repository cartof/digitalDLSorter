options(stringsAsFactors = FALSE)

library(optparse,quietly=T, warn.conflicts=F, verbose=F)

##############################
## DEFINE OPTION PARAMETERS ##
##############################

opls <- list()

opls[[1]] <- make_option(c('--probMatrixNames','-m'), action="store", type="character", dest="probMatrixNamesFile", help="File containing the names of the cells to mix for each bulk profile in rows [default: %default]")
opls[[2]] <- make_option(c('--simCounts','-s'), action="store", type="character", dest="simCountsFile", help="File containing counts per gene in rows and cell in columns [default: %default]")
opls[[3]] <- make_option(c('--outputPath','-o'), action="store", type="character", dest="outputPath", default=getwd(), help="Path to write output files to [default: %default]")
opls[[4]] <- make_option(c('--prefix','-p'), action="store", type="character", dest="prefix", help="Prefix for the output file [default: %default]")
opls[[5]] <- make_option(c('--nCores','-n'), action="store", type="integer", dest="nCores", default=1, help="Number of cores to parallelize [default: %default]")

opts <- OptionParser(usage = "usage: %prog [options]",
                     option_list = opls, add_help_option = TRUE,
                     prog = "generateBulkSamples.R", description = "", epilogue = "")

args <- parse_args(opts,
                   args = commandArgs(trailingOnly = TRUE),
                   print_help_and_exit = TRUE,
                   positional_arguments = FALSE)

cat("\n")
cat("COMMAND LINE ARGUMENTS\n")
cat(paste0(names(args),": ",args,"\n"))
cat("\n")
cat("\n")

if (length(args) < 6) {
  message("Incomplete arguments list!!!!!!!")
  print_help(opts)
  quit(status = 1)
}

suppressMessages(library(gtools))
suppressMessages(library(Matrix))
suppressMessages(library(pbapply))
suppressMessages(library(edgeR))

simCountsFile <- args$simCountsFile
probMatrixNamesFile <- args$probMatrixNamesFile
outputPath <- args$outputPath
prefix <- args$prefix
nCores <- args$nCores

setwd(outputPath)

cat("Load simCounts File\n")
if (grepl(".tsv",simCountsFile)) {
  if (grepl(".tsv.gz",simCountsFile)) {
    cat("\tTab Gzipped format\n")
    simCounts <- read.delim(file = gzfile(simCountsFile),sep = "\t",header = T,stringsAsFactors = F)
  } else {
    cat("\tTab plain format\n")
    simCounts <- read.delim(file = simCountsFile,sep = "\t",header = T,stringsAsFactors = F)
  }
} else if (grepl(".rds",simCountsFile)) {
  cat("\tRDS R object\n")
  simCounts <- readRDS(file = simCountsFile)
} else if (grepl(".mtx",simCountsFile)) {
  cat("\tSparse Matrix Market format\n")
  baseDir <- dirname(simCountsFile)
  if (!file.exists(file.path(baseDir,"genes.tsv"))) {
    message("No genes.tsv file with the matrix.mtx")
    quit(status=1)
  }
  if (!file.exists(file.path(baseDir,"barcodes.tsv"))) {
    message("No barcodes.tsv file with the matrix.mtx")
    quit(status = 1)
  } 
  simCounts <- readMM(simCountsFile)
  geneNames <- read.delim(file.path(baseDir,"genes.tsv"),header = F,sep = "\t",stringsAsFactors = F)
  rownames(simCounts) <- geneNames$V1
  cellNames <- read.delim(file.path(baseDir,"barcodes.tsv"),header = F,sep = "\t",stringsAsFactors = F)
  colnames(simCounts) <- cellNames$V1
}
cat(paste0(paste(c("Genes","Cells"),dim(simCounts))),"\n")
head(simCounts)[,1:6]
cat("\n")
simCounts <- cpm(simCounts)

cat("Load probMatrixNames File\n")
if (grepl(".tsv",probMatrixNamesFile)) {
  if (grepl(".tsv.gz",probMatrixNamesFile)) {
    cat("\tTab Gzipped format\n")
    probMatrixNames <- read.delim(gzfile(probMatrixNamesFile), header = T, sep = "\t")
  } else {
    cat("\tTab plain format\n")
    probMatrixNames <- read.delim(probMatrixNamesFile, header = T, sep = "\t")
  }
} else if (grepl(".rds",probMatrixNamesFile)) {
  cat("\tRDS R object\n")
  probMatrixNames <- readRDS(probMatrixNamesFile)  
}
cat(paste0(paste(c("Bulks","MixedCells"),dim(probMatrixNames))),"\n")
head(probMatrixNames)
cat("\n")



### GENERATE BULK COUNTS
setBulks <- function (x,c,i) {
  return(rowSums(c[,x]))
}

cat("GENERATE BULK COUNTS\n")
cat(paste("Cores",nCores,"\n"))
bulkCounts <- pbapply(probMatrixNames,1,FUN=setBulks,c=simCounts, cl = nCores)
cat(paste0(paste(c("Genes","Samples"),dim(bulkCounts))),"\n")
cat("DONE\n")
cat("\n")

colnames(bulkCounts) <- paste("Bulk",seq(dim(bulkCounts)[2]),sep = "_")

cat("WRITE BULK COUNTS\n")
saveRDS(bulkCounts,file = file.path(outputPath,paste(prefix,"simBulkCounts","rds",sep=".")))
gz <- gzfile(file.path(outputPath,paste(prefix,"simBulkCounts","tsv","gz",sep=".")),"w")
write.table(bulkCounts
            , file = gz
            , sep = "\t"
            , quote = F
            , row.names = T
            , col.names = T
            )
close(gz)
cat("DONE\n")
cat("\n")

bulkCounts <- NULL
rm(bulkCounts)

