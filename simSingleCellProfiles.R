options(stringsAsFactors = FALSE)

library(optparse,quietly=T, warn.conflicts=F, verbose=F)

##############################
## DEFINE OPTION PARAMETERS ##
##############################

opls <- list()

opls[[1]] <- make_option(c('--setType','-t'), action="store", type="character", dest="st", default="All", help="Cell Type to evaluate [default: %default]")
opls[[2]] <- make_option(c('--nCells','-n'), action="store", type="numeric", dest="nCells", default=100, help="Number of cells per cell type [default: %default]")
opls[[3]] <- make_option(c('--zinbParams','-z'), action="store", type="character", dest="zinbParamsFile", default=NULL, help="File containing parameters of Zinb model (.rds) [default: %default]")
opls[[4]] <- make_option(c('--cellsMetadata','-c'), action="store", type="character", dest="cellsMetadataFile", default=NULL, help="File containing metadata about cells in rows separated by tabs (.tsv, .tsv.gz, .rds) [default: %default]")
opls[[5]] <- make_option(c('--cellTypeColumn','-q'), action="store", type="character", dest="cellTypeColumn", default=NULL, help="Name of the column containing the cell classification in the cellsMetadataFile [default: %default]")
opls[[6]] <- make_option(c('--counts','-m'), action="store", type="character", dest="countsFile", default=NULL, help="File containing counts per gene in rows and cell in columns [default: %default]")
opls[[7]] <- make_option(c('--outputPath','-o'), action="store", type="character", dest="outputPath", default=getwd(), help="Path to write output files to [default: %default]")
opls[[8]] <- make_option(c('--prefix','-p'), action="store", type="character", dest="prefix", default=NULL, help="Prefix for the output file [default: %default]")

opts <- OptionParser(usage = "usage: %prog [options]",
                     option_list = opls, add_help_option = TRUE,
                     prog = "simSingleCellProfiles.R", description = "", epilogue = "")

args <- parse_args(opts,
                   args = commandArgs(trailingOnly = TRUE),
                   print_help_and_exit = TRUE,
                   positional_arguments = FALSE)

cat("\n")
cat("COMMAND LINE ARGUMENTS\n")
cat(paste0(names(args),": ",args,"\n"))
cat("\n")
cat("\n")

if (length(args) < 9) {
  message("Incomplete arguments list!!!!!!!")
  print_help(opts)
  quit()
}

suppressMessages(library(splatter))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(zinbwave))


setType <- args$st
nCells <- args$nCells
cellsMetadataFile <- args$cellsMetadataFile
origCountsFile <- args$countsFile
zinbParamsFile <- args$zinbParamsFile
cellTypeColumnName <- args$cellTypeColumn
outputPath <- args$outputPath
prefix <- args$prefix

cat("Load Cells Metadata\n")
if (grepl(".tsv",cellsMetadataFile)) {
  if (grepl(".tsv.gz",cellsMetadataFile)) {
    cat("Tab Gzipped format\n")
    cellsMetadata <- read.delim(file = gzfile(cellsMetadataFile),sep = "\t",header = T,stringsAsFactors = F)
  } else {
    cat("Tab plain format\n")
    cellsMetadata <- read.delim(file = cellsMetadataFile, sep = "\t",header = T,stringsAsFactors = F)
  }
} else if (grepl(".rds",cellsMetadataFile)) {
  cat("RDS R object\n")
  cellsMetadata <- readRDS(file = cellsMetadataFile)
}
rownames(cellsMetadata) <- paste(cellsMetadata[,cellTypeColumnName],cellsMetadata$Cell_ID,sep="_")
cat(paste(c("Cells","Metadata"),dim(cellsMetadata)),"\n")
head(cellsMetadata)
cat("\n")

cat("Load Original Counts Matrix\n")
if (grepl(".tsv",origCountsFile)) {
  if (grepl(".tsv.gz",origCountsFile)) {
    cat("Tab Gzipped format\n")
    origCounts <- read.delim(file = gzfile(origCountsFile),sep = "\t",header = T,stringsAsFactors = F)
  } else {
    cat("Tab plain format\n")
    origCounts <- read.delim(file = origCountsFile,sep = "\t",header = T,stringsAsFactors = F)
  }
} else if (grepl(".rds",origCountsFile)) {
  cat("RDS R object\n")
  origCounts <- readRDS(file = origCountsFile)
} else if (grepl(".mtx",origCountsFile)) {
  cat("Sparse Matrix Market format\n")
  baseDir <- dirname(origCountsFile)
  if (!file.exists(file.path(baseDir,"genes.tsv"))) exit("No genes.tsv file with the matrix.mtx")
  if (!file.exists(file.path(baseDir,"barcodes.tsv"))) exit("No barcodes.tsv file with the matrix.mtx")
  origCounts <- readMM(origCountsFile)
  geneNames <- read.delim(file.path(baseDir,"genes.tsv"),header = F,sep = "\t",stringsAsFactors = F)
  rownames(origCounts) <- geneNames$V1
  cellNames <- read.delim(file.path(baseDir,"barcodes.tsv"),header = F,sep = "\t",stringsAsFactors = F)
  colnames(origCounts) <- cellNames$V1
}
cat("origCounts\n")
cat(paste(c("Genes","Cells"),dim(origCounts)),"\n")
head(origCounts)[,1:6]
cat("\n")

cat("Filter Original Counts File for cells used using cellsMetadata\n")
origCounts <- ceiling(origCounts[,cellsMetadata$Cell_ID])
colnames(origCounts) <- paste(cellsMetadata[,cellTypeColumnName],cellsMetadata$Cell_ID,sep = "_")
cat(paste(c("Genes","Cells"),dim(origCounts)),"\n")
head(origCounts)[,1:6]
cellsMetadata$simCellName <- paste(cellsMetadata[,cellTypeColumnName],cellsMetadata$Cell_ID,sep = "_")
cellsMetadata$Simulated <- FALSE

cat("Load ZinbParams RDS Object\n")
zinbParams <- readRDS(zinbParamsFile)

cellSetNames <- NULL
modelCellTypes <- grep(cellTypeColumnName,colnames(zinbParams@model@X),value = T)
cellTypeNames <- sub(cellTypeColumnName,"",modelCellTypes)
names(cellTypeNames) <- modelCellTypes
cat("Cell Types in Model\n")
cat(modelCellTypes,"\n")
cat("\n")

for (s in modelCellTypes) {
  cellTypeName <- cellTypeNames[s]
  cat(paste(s,cellTypeName),"\n")
  cellIndex <- rownames(zinbParams@model@X)[which(zinbParams@model@X[,s] == 1)]
  nams <- sample(cellIndex,size = nCells,replace = T)
  if (is.null(cellSetNames)) {
    cellSetNames <- nams
    names(cellSetNames) <- paste(cellTypeName,"_S",seq(from = 1, to = nCells),sep="")
  } else {
    ns <- names(cellSetNames)
    cellSetNames <- c(cellSetNames,nams)
    names(cellSetNames) <- c(ns,paste(cellTypeName,"_S",seq(from = length(ns) + 1, to = length(ns) + nCells),sep=""))
  }
}

interceptCellTypeName <- setdiff(levels(factor(cellsMetadata[,cellTypeColumnName])),cellTypeNames)
cat(interceptCellTypeName,"\n") # To get the intercept cell type the rowSum of all FinalCellType columns should be 0
cellIndex <- rownames(zinbParams@model@X)[rowSums(zinbParams@model@X[,grep(cellTypeColumnName,colnames(zinbParams@model@X),value = T)]) == 0]
nams <- sample(cellIndex,size = nCells,replace = T)
ns <- names(cellSetNames)
cellSetNames <- c(cellSetNames,nams)
names(cellSetNames) <- c(ns,paste(interceptCellTypeName,"_S",seq(from = length(ns) + 1, to = length(ns) + nCells),sep=""))

nCells <- length(cellSetNames)

cat("Get params from model\n")
mu <- getMu(zinbParams@model) # rows are cells
cat(paste("mu",dim(mu)),"\n")
pi <- getPi(zinbParams@model) # rows are cells
cat(paste("pi",dim(pi)),"\n")

mu <- mu[cellSetNames,]
pi <- pi[cellSetNames,]

theta <- getTheta(zinbParams@model) # for genes
cat(paste("Theta",length(theta)),"\n")

n <- length(cellSetNames) # as.numeric(nCells)
J <- nFeatures(zinbParams@model)
cat("Simulated Matrix dimensions\n")
cat(paste("n",n),"\n")
cat(paste("J",J),"\n")
i <- seq(n*J)
cat(paste("i",length(i)),"\n")

datanb <- rnbinom(length(i), mu = mu[i], size = theta[ceiling(i/n)])
data.nb <- matrix(datanb, nrow = n)

datado <- rbinom(length(i), size=1, prob = pi[i])
data.dropout <- matrix(datado, nrow = n)

simCounts <- data.nb * (1 - data.dropout)
colnames(simCounts) <- colnames(mu)
simCounts <- t(simCounts)

simCounts <- Matrix(simCounts,dimnames = list(rownames=rownames(zinbParams@model@V),colnames=names(cellSetNames)))

cat(paste(c("Genes","Cells"),dim(simCounts)),"\n")
head(simCounts)[,1:5]

simCellsMetadata <- cellsMetadata[cellSetNames,]
simCellsMetadata$simCellName <- names(cellSetNames)
simCellsMetadata$Simulated <- TRUE

# Save cells Metadata together with the original real cells
simCellsMetadata <- rbind(cellsMetadata,simCellsMetadata)
cat(paste("Save cellsMetadata in",file.path(outputPath,paste(prefix,setType,"simCellsMetadata",nCells,"tsv","gz",sep="."))),"\n")
saveRDS(simCellsMetadata,file = file.path(outputPath,paste(prefix,setType,"simCellsMetadata",nCells,"rds",sep=".")))
gz1 <- gzfile(file.path(outputPath,paste(prefix,setType,"simCellsMetadata",nCells,"tsv","gz",sep=".")), "w")
write.table(simCellsMetadata,file = gz1, sep = "\t",row.names = F,col.names = T,quote = F)
close(gz1)

### Write SimData in SparseMatrix together with the original real cells
cat(paste("Save Matrix Market Counts in",file.path(outputPath,paste(prefix,setType,"simCellsCountsMtx",nCells,sep = "."))),"\n")
simCounts <- Matrix(as.matrix(simCounts),sparse = T)
origCounts <- Matrix(as.matrix(origCounts),sparse = T)
simCounts <- cbind(origCounts[rownames(simCounts),],simCounts)

dir.create(file.path(outputPath,paste(prefix,setType,"simCellsCountsMtx",nCells,sep = ".")))
writeMM(simCounts
        ,file = file.path(outputPath,paste(prefix,setType,"simCellsCountsMtx",nCells,sep = "."),"matrix.mtx"))

write.table(x = rownames(simCounts),file = file.path(outputPath,paste(prefix,setType,"simCellsCountsMtx",nCells,sep = "."),"genes.tsv"),sep="\t",quote = F,col.names = F, row.names = F)

write.table(x = colnames(simCounts),file = file.path(outputPath,paste(prefix,setType,"simCellsCountsMtx",nCells,sep = "."),"barcodes.tsv"),sep="\t",quote = F,col.names = F, row.names = F)

cat("DONE\n")
