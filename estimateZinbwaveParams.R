options(stringsAsFactors = FALSE)

library(optparse,quietly=T, warn.conflicts=F, verbose=F)

##############################
## DEFINE OPTION PARAMETERS ##
##############################

opls <- list()

opls[[1]] <- make_option(c('--setType','-t'), action="store", type="character", dest="st", default="All", help="Cell Type to evaluate [default: %default]")
opls[[2]] <- make_option(c('--cellsMetadataFile','-m'), action="store", type="character", dest="cellsMetadataFile", default=NULL, help="File containing metadata about cells in rows separated by tabs (.tsv, .tsv.gz, .rds) [default: %default]")
opls[[3]] <- make_option(c('--cellTypeColumn'), action="store", type="character", dest="cellTypeColumn", default=NULL, help="Name of the column containing the cell classification in the cellsMetadataFile [default: %default]")
opls[[4]] <- make_option(c('--cellIDColumn'), action="store", type="character", dest="cellIDColumn", default=NULL, help="Name of the column containing the cell ID in the cellsMetadataFile, cell IDs should be the same as column names in counts file [default: %default]")
opls[[5]] <- make_option(c('--cellCovColumns'), action="store", type="character", dest="cellCovColumns", default=NULL, help="List separated by commas of column names from cellsMetadataFile to consider as covariates for the model [default: %default]")
opls[[6]] <- make_option(c('--genesAnnotFile','-g'), action="store", type="character", dest="genesAnnotFile", default=NULL, help="File containing metadata about genes in rows separated by tabs (.tsv, .tsv.gz, .rds) [default: %default]")
opls[[7]] <- make_option(c('--geneIDColumn'), action="store", type="character", dest="geneIDColumn", default=NULL, help="Name of the column containing the gene ID in the genesAnnotFile, gene IDs should be the same as row names in counts file [default: %default]")
opls[[8]] <- make_option(c('--geneCovColumns'), action="store", type="character", dest="geneCovColumns", default=NULL, help="List separated by commas of column names from genesAnnotFile to consider as covariates for the model [default: %default]")
opls[[9]] <- make_option(c('--countsFile','-c'), action="store", type="character", dest="countsFile", default=NULL, help="File containing counts per gene in rows and cells in columns (.tsv, .tsv.gz, .mtx, .rds) [default: %default]")
opls[[10]] <- make_option(c('--minCounts'), action="store", type="numeric", dest="minCounts", default=0, help="Minimum gene counts to filter [default: %default]")
opls[[11]] <- make_option(c('--minCells'), action="store", type="numeric", dest="minCells", default=1, help="Minimum of cells with more than minCounts [default: %default]")
opls[[12]] <- make_option(c('--outputPath','-o'), action="store", type="character", dest="outputPath", default=getwd(), help="Path to write output files to [default: %default]")
opls[[13]] <- make_option(c('--prefix','-p'), action="store", type="character", dest="prefix", default=NULL, help="Prefix for the output file [default: %default]")
opls[[14]] <- make_option(c('--threads','-d'), action="store", type="numeric", dest="np", default=1, help="Number of threads [default: %default]")

opts <- OptionParser(usage = "usage: %prog [options]",
                     option_list = opls, add_help_option = TRUE,
                     prog = "estimateZinbwaveParams.R", description = "", epilogue = "")

args <- parse_args(opts,
                   args = commandArgs(trailingOnly = TRUE),
                   print_help_and_exit = TRUE,
                   positional_arguments = FALSE)

cat("\n")
cat("COMMAND LINE ARGUMENTS\n")
cat(paste0(names(args),": ",args,"\n"))
cat("\n")
cat("\n")

if (length(args) < 10) {
  message("Incomplete arguments list!!!!!!!")
  print_help(opts)
  quit()
}

suppressMessages(library(splatter))

st <- args$st
samplesMetadataFile <- args$samplesMetadataFile
cellsMetadataFile <- args$cellsMetadataFile
cellTypeColumn <- args$cellTypeColumn
cellIDColumn <- args$cellIDColumn
cellCovColumns <- args$cellCovColumns

genesAnnotFile <- args$genesAnnotFile
geneIDColumn <- args$geneIDColumn
geneCovColumns <- args$geneCovColumns

countsFile <- args$countsFile
minCounts <- args$minCounts
minCells <- args$minCells

outputPath <- args$outputPath
prefix <- args$prefix
np <- args$np

cellCovColumns <- strsplit(x = cellCovColumns,split = ",")
geneCovColumns <- strsplit(x = geneCovColumns,split = ",")

cat(paste("Parallel",np,"setType",st),"\n")

cat("Load Cells Metadata\n")
if (grepl(".tsv",cellsMetadataFile)) {
  if (grepl(".tsv.gz",cellsMetadataFile)) {
    cat("\tTab Gzipped format\n")
    cellsMetadata <- read.delim(file = gzfile(cellsMetadataFile),sep = "\t",header = T,stringsAsFactors = F)
  } else {
    cat("\tTab plain format\n")
    cellsMetadata <- read.delim(file = cellsMetadataFile,sep = "\t",header = T,stringsAsFactors = F)
  }
} else if (grepl(".rds",cellsMetadataFile)) {
  cat("\tRDS R object\n")
  cellsMetadata <- readRDS(file = cellsMetadataFile)
}

if (!cellIDColumn %in% colnames(cellsMetadata)) {stop(paste(cellIDColumn, "column is not present in",cellsMetadataFile))}

if (!cellTypeColumn %in% colnames(cellsMetadata)) {stop(paste(cellTypeColumn, "column is not present in",cellsMetadataFile))}

if (any(!cellCovColumns %in% colnames(cellsMetadata))) {stop(paste("Some columns of",cellCovColumns,"are not present in",cellsMetadataFile))}

rownames(cellsMetadata) <- paste(cellsMetadata[,cellTypeColumn],cellsMetadata[,cellIDColumn],sep = "_")

cat(paste0(paste(c("Cells","Metadata"),dim(cellsMetadata))),"\n")
table(cellsMetadata[,cellTypeColumn])
cat("\n")



cat("Load Genes Metadata\n")
if (grepl(".tsv",genesAnnotFile)) {
  if (grepl(".tsv.gz",genesAnnotFile)) {
    cat("\tTab Gzipped format\n")
    genesAnnot <- read.delim(file = gzfile(genesAnnotFile),sep = "\t",header = T,stringsAsFactors = F)
  } else {
    cat("\tTab plain format\n")
    genesAnnot <- read.delim(file = genesAnnotFile,sep = "\t",header = T,stringsAsFactors = F)
  }
} else if (grepl(".rds",genesAnnotFile)) {
  cat("\tRDS R object\n")
  genesAnnot <- readRDS(file = genesAnnotFile)
}

cat(paste0(paste(c("Genes","Metadata"),dim(genesAnnot))),"\n")
head(genesAnnot)
cat("\n")

if (!geneIDColumn %in% colnames(genesAnnot)) {stop(paste(geneIDColumn, "column is not present in",genesAnnotFile))}

if (any(!geneCovColumns %in% colnames(genesAnnot))) {stop(paste("Some columns of",geneCovColumns,"are not present in",genesAnnotFile))}

cat("Load Counts\n")
if (grepl(".tsv",countsFile)) {
  if (grepl(".tsv.gz",countsFile)) {
    cat("\tTab Gzipped format\n")
    counts <- read.delim(file = gzfile(countsFile),sep = "\t",header = T,stringsAsFactors = F)
  } else {
    cat("\tTab plain format\n")
    counts <- read.delim(file = countsFile,sep = "\t",header = T,stringsAsFactors = F)
  }
} else if (grepl(".rds",countsFile)) {
  cat("\tRDS R object\n")
  counts <- readRDS(file = countsFile)
} else if (grepl(".mtx",countsFile)) {
  cat("\tSparse Matrix Market format\n")
  baseDir <- dirname(countsFile)
  if (!file.exists(file.path(baseDir,"genes.tsv"))) {
    message("No genes.tsv file with the matrix.mtx")
    quit(status = 1)
  } 
  if (!file.exists(file.path(baseDir,"barcodes.tsv"))) {
    message("No barcodes.tsv file with the matrix.mtx")
    quit(status = 1)
  } 
  counts <- readMM(countsFile)
  geneNames <- read.delim(file.path(baseDir,"genes.tsv"),header = F,sep = "\t",stringsAsFactors = F)
  rownames(counts) <- geneNames$V1
  cellNames <- read.delim(file.path(baseDir,"barcodes.tsv"),header = F,sep = "\t",stringsAsFactors = F)
  colnames(counts) <- cellNames$V1
}
cat(paste("Loaded Counts Matrix from",countsFilem,"\n"))
cat(paste0(paste(c("Genes","Cells"),dim(counts))),"\n")
head(counts)[,1:6]
cat("\n")

commonGenes <- intersect(rownames(counts),geneAnnot[,geneIDColumn])

if (length(commonGenes) < min(dim(counts)[1],dim(genesAnnot)[1])) {stop(paste("Some genes from",genesAnnotFile,"and",countsFile,"do not match"))}

counts <- counts[commonGenes,]
counts <- as.matrix(counts)
cat(paste("Counts Matrix after selecting for genes in",genesAnnotFile,"\n"))
cat(paste(c("Genes","Cells"),dim(counts)),"\n")
cat("\n")

commonCells <- intersect(colnames(counts),cellsMetadata[,cellIDColumn])
if (length(commonCells) < min(dim(counts)[2],dim(cellsMetadataFile)[1])) {stop(paste("Some cells from",cellsMetadataFile,"and",countsFile,"do not match"))}

counts <- counts[,commonCells]
counts <- as.matrix(counts)
cat(paste("Counts Matrix after selecting for cells in",cellsMetadataFile,"\n"))
cat(paste(c("Genes","Cells"),dim(counts)),"\n")
cat("\n")

minCounts <- 0
minCells <- 1

dir.create(file.path(outputPath))

if (file.exists(file = file.path(outputPath,paste(prefix,"genesSelected","txt","gz",sep=".")))) {
  
  cat(paste("Load list of selected genes:from",paste(prefix,"genesSelected","txt","gz",sep="."),"\n"))
  genesSelected <- read.delim(file = file.path(outputPath,paste(prefix,"genesSelected","txt","gz",sep=".")),sep = "\t",header = F)
  cat(paste("Genes:",length(genesSelected),"\n"))
  head(genesSelected)
  cat("\n")
  
} else {
  
  cat(paste("Filter Genes by minCounts:",minCounts,"in minCells:",minCells,"\n"))
  counts <- counts[rowSums(counts > minCounts) >= minCells,]
  cat("Filtered genes and selected cells\n")
  cat(paste(c("Genes","Cells"),dim(counts)),"\n")
  
  genesSelected <- rownames(counts)
  genesSelected <- genesSelected[!genesSelected == "NA"]
  cat("Genes Selected: ",length(genesSelected),"\n")
  
  write.table(genesSelected,file = file.path(outputPath,paste(prefix,"genesSelected","txt","gz",sep=".")),sep = "\t",quote = "F",row.names = F,col.names = F)  
  cat("Saved to: ",file.path(outputPath,paste(prefix,"genesSelected","txt","gz",sep=".")),"\n")
  cat("\n")
  
}

cat("Filter counts on genes selected\n")
counts <- counts[genesSelected,]
cat(paste(c("Genes","Cells"),dim(counts)),"\n")
cat("\n")

cat("Estimate zinb Parameters for",st,"\n")

suppressMessages(library(BiocParallel))

cat("Set parallel env to",np,"\n")
snowParam <- SnowParam(workers = np, type = "SOCK")
cat("\n")


if (st == "All") {
  
  cat("Estimate parameters for the whole experiment\n")
  counts <- as.matrix(counts)
  counts <- counts[rowSums(counts) > 0,]
  cat(paste(dim(counts)),"\n")
  cat("\n")
  
  cat("Estimate parameters for experiment with model matrix\n")
  cat(paste("Create model matrix based on:",cellCovColumns,"and",cellTypeColumn,"\n"))
  sdm <- model.matrix(as.formula(paste("~",paste(cellCovColumns,cellTypeColumn,sep="+"))), data = cellsMetadata[match(colnames(counts),cellsMetadata[,cellIDColumn]),])
  cat(paste(dim(sdm)),"\n")
  print(head(sdm))
  cat("\n")
  
  if (length(geneCovColumns)) {
    cat(paste("Create gene model matrix with",geneCovColumns,"Covariates\n"))
    gdm <- model.matrix(as.formula(paste("~",paste(geneCovColumns,sep="+"))), data = genesAnnot[match(rownames(counts),genesAnnot[,geneIDColumn]),])
  } else {
    cat(paste("Create gene model matrix without Covariates\n"))
    gdm <- model.matrix(~1, data = genesAnnot[match(rownames(counts),genesAnnot[,geneIDColumn]),])  
  }
  
  rownames(gdm) <- rownames(counts)
  cat(paste(dim(gdm)),"\n")
  print(head(gdm))
  cat("\n")
  
  cat("Run estimation process\n")
  start_time <- Sys.time()
  cat(paste(start_time,"\n"))
  zinbParamsAll <- zinbEstimate(ceiling(counts)
                                , BPPARAM = snowParam
                                , design.samples = sdm
                                , design.genes = gdm
                                , O_mu = matrix(0,nrow = ncol(counts), ncol = nrow(counts)
                                                , dimnames = list(rownames = seq(ncol(counts))
                                                                  ,colnames = rownames(counts)
                                                )
                                )
                                , O_pi = matrix(0,nrow = ncol(counts), ncol = nrow(counts)
                                                , dimnames = list(rownames = seq(ncol(counts))
                                                                  ,colnames = rownames(counts)
                                                )
                                )
                                , beta_mu = matrix(0,nrow = ncol(sdm), ncol = nrow(counts)
                                                   , dimnames = list(rownames = colnames(sdm)
                                                                     ,colnames = rownames(counts)
                                                   )
                                )
                                , beta_pi = matrix(0,nrow = ncol(sdm), ncol = nrow(counts)
                                                   , dimnames = list(rownames = colnames(sdm)
                                                                     ,colnames = rownames(counts)
                                                   )
                                )
                                , alpha_mu = matrix(0,nrow = 0, ncol = nrow(counts)
                                                    , dimnames = list(rownames = NULL
                                                                      ,colnames = rownames(counts)
                                                    )
                                )
                                , alpha_pi = matrix(0,nrow = 0, ncol = nrow(counts)
                                                    , dimnames = list(rownames = NULL
                                                                      ,colnames = rownames(counts)
                                                    )
                                )
  )
  
  cat("Save model parameters for",st,"cells to:",file.path(outputPath,paste(prefix,st,"zinbParams","rds",sep=".")),"\n")
  saveRDS(zinbParamsAll,file=file.path(outputPath,paste(prefix,st,"zinbParams","rds",sep=".")))
  
  end_time <- Sys.time()
  
  cat("DONE\n")
  cat((end_time - start_time),"\n")
  cat("\n")
  
} else {
  
  cat(paste("Estimate parameters for",st,"from the experiment\n"))
  cat("Collect counts for",st,"cells\n")
  cellIDs <- cellsMetadata[which(cellsMetadata[,cellTypeColumn] == st),cellIDColumn]
  counts <- counts[,cellIDs]
  counts <- as.matrix(counts)
  counts <- counts[rowSums(counts) > 0,]
  cat(paste(c("Genes","Cells"),dim(counts)),"\n")
  cat("\n")
  
  cat("Run estimation process\n")
  start_time <- Sys.time()
  cat(start_time,"\n")
  
  cat("Create model matrix based on Gene Length","\n")
  if (length(geneCovColumns)) {
    cat(paste("Create gene model matrix with",geneCovColumns,"Covariates\n"))
    gdm <- model.matrix(as.formula(paste("~",paste(geneCovColumns,sep="+"))), data = genesAnnot[match(rownames(counts),genesAnnot[,geneIDColumn]),])
  } else {
    cat(paste("Create gene model matrix without Covariates\n"))
    gdm <- model.matrix(~1, data = genesAnnot[match(rownames(counts),genesAnnot[,geneIDColumn]),])  
  }
  
  rownames(gdm) <- rownames(counts)
  print(head(gdm))
  cat("\n")
  
  zinbParams <- zinbEstimate(ceiling(counts)
                             , BPPARAM = snowParam
                             , design.genes = gdm
                             , O_mu = matrix(0,nrow = ncol(counts), ncol = nrow(counts)
                                             , dimnames = list(rownames = seq(ncol(counts))
                                                               ,colnames = rownames(counts)
                                             )
                             )
                             , O_pi = matrix(0,nrow = ncol(counts), ncol = nrow(counts)
                                             , dimnames = list(rownames = seq(ncol(counts))
                                                               ,colnames = rownames(counts)
                                             )
                             )
                             , beta_mu = matrix(0,nrow = 1, ncol = nrow(counts)
                                                , dimnames = list(rownames = seq(1)
                                                                  ,colnames = rownames(counts)
                                                )
                             )
                             , beta_pi = matrix(0,nrow = 1, ncol = nrow(counts)
                                                , dimnames = list(rownames = seq(1)
                                                                  ,colnames = rownames(counts)
                                                )
                             )
                             , alpha_mu = matrix(0,nrow = 0, ncol = nrow(counts)
                                                 , dimnames = list(rownames = NULL
                                                                   ,colnames = rownames(counts)
                                                 )
                             )
                             , alpha_pi = matrix(0,nrow = 0, ncol = nrow(counts)
                                                 , dimnames = list(rownames = NULL
                                                                   ,colnames = rownames(counts)
                                                 )
                             )
  )
  
  cat("Save model parameters for",st,"cells to:",file.path(outputPath,paste(prefix,st,"zinbParams","rds",sep=".")),"\n")
  saveRDS(zinbParams,file=file.path(outputPath,paste(prefix,st,"zinbParams","rds",sep=".")))
  
  end_time <- Sys.time()
  
  cat("DONE\n")
  cat((end_time - start_time),"\n")
  cat("\n")
  
}

