options(stringsAsFactors = FALSE)

library(optparse,quietly=T, warn.conflicts=F, verbose=F)

##############################
## DEFINE OPTION PARAMETERS ##
##############################

opls <- list()

opls[[1]] <- make_option(c('--probMatrix','-m'), action="store", type="character", dest="probMatrixFile", help="File containing the prob vectors for each simulated bulk mix in rows. Header with names of cell types [default: %default]")
opls[[2]] <- make_option(c('--cellSetList','-c'), action="store", type="character", dest="cellSetListFile", help="File containing the names of the cells selected to generate the bulk mixes [default: %default]")
opls[[3]] <- make_option(c('--bulkSimCounts','-b'), action="store", type="character", dest="bulkSimCountsFile", help="File containing the counts of the genes in rows for each bulk mix in columns. rownames with gene names and header with bulk sample names [default: %default]")
opls[[4]] <- make_option(c('--scSimCounts','-s'), action="store", type="character", dest="scSimCountsFile", help="File containing the counts of the genes in rows for each single cell. Rows with gene names and header with cell names [default: %default]")
opls[[5]] <- make_option(c('--outputPath','-o'), action="store", type="character", dest="outputPath", default=getwd(), help="Path to write output files to [default: %default]")
opls[[6]] <- make_option(c('--prefix','-p'), action="store", type="character", dest="prefix", help="Prefix for the output file [default: %default]")

opts <- OptionParser(usage = "usage: %prog [options]",
                     option_list = opls, add_help_option = TRUE,
                     prog = "generateCombinedScaledSuffledSet.R", description = "", epilogue = "")

args <- parse_args(opts,
                   args = commandArgs(trailingOnly = TRUE),
                   print_help_and_exit = TRUE,
                   positional_arguments = FALSE)

cat("\n")
cat("COMMAND LINE ARGUMENTS\n")
cat(paste0(names(args),": ",args,"\n"))
cat("\n")
cat("\n")

if (length(args) < 7) {
  message("Incomplete arguments list!!!!!!!")
  print_help(opts)
  quit(status = 1)
}

suppressMessages(library(Matrix))
suppressMessages(library(edgeR))

bulkSimCountsFile <- args$bulkSimCountsFile
probMatrixFile <- args$probMatrixFile
scSimCountsFile <- args$scSimCountsFile
cellSetListFile <- args$cellSetListFile
outputPath <- args$outputPath
prefix <- args$prefix

cat("Load Bulk Simulated Set\n")
bulkSimCounts <- readRDS(bulkSimCountsFile)
bulkProbsMatrix <- readRDS(probMatrixFile)
cellSetList <- readRDS(cellSetListFile)
cat(paste(c("Genes:","BulkSamples:"),dim(bulkSimCounts)),"\n")
cat(paste(c("BulkSamples:","CellTypeProbs:"),dim(bulkProbsMatrix)),"\n")
cat(paste("Single Cells Used:",length(unlist(cellSetList))),"\n")

cat("Load Single Cell simCounts File\n")
if (grepl(".tsv",scSimCountsFile)) {
  if (grepl(".tsv.gz",scSimCountsFile)) {
    cat("\tTab Gzipped format\n")
    scSimCounts <- read.delim(file = gzfile(scSimCountsFile),sep = "\t",header = T,stringsAsFactors = F)
  } else {
    cat("\tTab plain format\n")
    scSimCounts <- read.delim(file = scSimCountsFile,sep = "\t",header = T,stringsAsFactors = F)
  }
} else if (grepl(".rds",scSimCountsFile)) {
  cat("\tRDS R object\n")
  scSimCounts <- readRDS(file = scSimCountsFile)
} else if (grepl(".mtx",scSimCountsFile)) {
  cat("\tSparse Matrix Market format\n")
  baseDir <- dirname(scSimCountsFile)
  if (!file.exists(file.path(baseDir,"genes.tsv"))) {
    message("No genes.tsv file with the matrix.mtx")
    quit(status=1)
  }
  if (!file.exists(file.path(baseDir,"barcodes.tsv"))) {
    message("No barcodes.tsv file with the matrix.mtx")
    quit(status = 1)
  } 
  scSimCounts <- readMM(scSimCountsFile)
  geneNames <- read.delim(file.path(baseDir,"genes.tsv"),header = F,sep = "\t",stringsAsFactors = F)
  rownames(scSimCounts) <- geneNames$V1
  cellNames <- read.delim(file.path(baseDir,"barcodes.tsv"),header = F,sep = "\t",stringsAsFactors = F)
  colnames(scSimCounts) <- cellNames$V1
}
cat(paste(c("Genes:","Sinlge Cells:"),dim(scSimCounts)),"\n")
head(scSimCounts)[,1:6]
cat("\n")

### COMBINE bulkCounts and scCounts
geneList <- intersect(rownames(bulkSimCounts),rownames(scSimCounts))
cat(paste("Common Genes",length(geneList)),"\n")
Counts <- cbind(scSimCounts[geneList,unlist(cellSetList)],bulkSimCounts[geneList,])

saveRDS(Counts,file = file.path(outputPath,paste(prefix,"combinedCounts","rds",sep=".")))
gz <- gzfile(file.path(outputPath,paste(prefix,"combinedCounts","tsv","gz",sep=".")),"w")
write.table(as.matrix(Counts)
            ,file = gz
            ,sep = "\t"
            ,quote = F
            ,row.names = T
            ,col.names = T
)
close(gz)

### Generate SC probs matrix
tpsm <- matrix(unlist(sapply(names(cellSetList),function (x,l) {
  v <- rep(0,length(l))
  names(v) <- names(l)
  v[x] <- 100
  return(rep(v,length(l[[x]])))
}
,l=cellSetList
)),ncol=length(cellSetList),byrow=T)
colnames(tpsm) <- names(cellSetList)

tpsm <- tpsm[,colnames(bulkProbsMatrix)]

rownames(tpsm) <- unlist(cellSetList)

probsMatrix <- rbind(tpsm,bulkProbsMatrix)/100

rownames(probsMatrix) <- c(rownames(tpsm),colnames(bulkSimCounts))

saveRDS(probsMatrix,file = file.path(outputPath,paste(prefix,"combinedProbsMatrix","rds",sep=".")))

gz <- gzfile(file.path(outputPath,paste(prefix,"combinedProbsMatrix","tsv","gz",sep=".")),"w")
write.table(probsMatrix
            ,file = gz
            ,sep = "\t"
            ,quote = F
            ,row.names = T
            ,col.names = T
)
close(gz)

### Scale Counts Matrix
Counts <- cpm(Counts,log = T,prior.count = 1)
Counts <- scale(Counts)

s <- sample(seq(dim(Counts)[2]))

Counts <- Counts[,s]

sampleNames <- colnames(Counts)

probsMatrix <- probsMatrix[s,]

Counts <- t(Counts)
paste(c("Genes:","Samples:"),dim(Counts))

saveRDS(Counts, file = file.path(outputPath,paste(prefix,"combinedCounts","log2CPMScaledShuffledTransposed","rds",sep = ".")))
gz <- gzfile(file.path(outputPath,paste(prefix,"combinedCounts","log2CPMScaledShuffledTransposed","tsv","gz",sep = ".")),"w")
write.table(Counts
            ,file = gz
            ,sep = "\t"
            ,row.names = T
            ,col.names = T
            ,quote = F
            )
close(gz)

saveRDS(probsMatrix, file = file.path(outputPath,paste(prefix,"combinedProbsMatrix","log2CPMScaledShuffledTransposed","rds",sep = ".")))
gz <- gzfile(file.path(outputPath,paste(prefix,"combinedProbsMatrix","log2CPMScaledShuffledTransposed","tsv","gz",sep = ".")),"w")
write.table(probsMatrix
            ,file = gz
            ,sep = "\t"
            ,row.names = T
            ,col.names = T
            ,quote = F
            )
close(gz)

saveRDS(sampleNames,file = file.path(outputPath,paste(prefix,"combinedCounts","ShuffledSampleNames","rds",sep = ".")))
gz <- gzfile(file.path(outputPath,paste(prefix,"combinedCounts","ShuffledSampleNames","tsv","gz",sep = ".")),"w")
write.table(sampleNames
            ,file = gz
            ,sep = "\t"
            ,row.names = F
            ,col.names = F
            ,quote = F
            )
close(gz)

saveRDS(geneList,file = file.path(outputPath,paste(prefix,"combinedCounts","geneList","rds",sep = ".")))
gz <- gzfile(file.path(outputPath,paste(prefix,"combinedCounts","geneList","tsv","gz",sep = ".")),"w")
write.table(geneList
            ,file = gz
            ,sep = "\t"
            ,row.names = F
            ,col.names = F
            ,quote = F
)
close(gz)

