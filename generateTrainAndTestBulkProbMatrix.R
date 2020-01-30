options(stringsAsFactors = FALSE)

library(optparse,quietly=T, warn.conflicts=F, verbose=F)

##############################
## DEFINE OPTION PARAMETERS ##
##############################

opls <- list()

opls[[1]] <- make_option(c('--freq','-f'), action="store", type="numeric", dest="trainFreq", default=0.65, help="Proportion of cells used for training set [default: %default]")
opls[[2]] <- make_option(c('--cellsMetadata','-c'), action="store", type="character", dest="cellsMetadataFile", default=NULL, help="File containing metadata about cells in rows separated by tabs (.tsv, .tsv.gz, .rds) [default: %default]")
opls[[3]] <- make_option(c('--genesAnnotFile','-g'), action="store", type="character", dest="geneAnnotFile", default=NULL, help="File containing metadata about genes in rows separated by tabs (.tsv, .tsv.gz, .rds) [default: %default]")
opls[[4]] <- make_option(c('--simCounts','-s'), action="store", type="character", dest="simCountsFile", default=NULL, help="File containing counts per gene in rows and cell in columns [default: %default]")
opls[[5]] <- make_option(c('--outputPath','-o'), action="store", type="character", dest="outputPath", default=getwd(), help="Path to write output files to [default: %default]")
opls[[6]] <- make_option(c('--design','-d'), action="store", type="character", dest="probsDesignFile", default=NULL, help="File containing probability range for each cell type (column headers: cellType, from, to) [default: %default]")
opls[[7]] <- make_option(c('--prefix','-p'), action="store", type="character", dest="prefix", help="Prefix for the output file [default: %default]")
opls[[8]] <- make_option(c('--cellTypeColumn','-q'), action="store", type="character", dest="cellTypeColumn", default=NULL, help="Name of the column containing the cell classification in the cellsMetadataFile [default: %default]")
#opls[[9]] <- make_option(c('--simCells','-n'), action="store", type="numeric", dest="sCells", default=NULL, help="Number of simulated cells [default: %default]")

opts <- OptionParser(usage = "usage: %prog [options]",
                     option_list = opls, add_help_option = TRUE,
                     prog = "generateTrainAndTestBulkProbMatrix.R", description = "", epilogue = "")

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

suppressMessages(library(gtools))
suppressMessages(library(Matrix))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))

adjustHundred <- function (x) {
  d <- abs(sum(x)-100)
  if (sum(x) > 100) {
    x[which(x >= d)[1]] <- x[which(x >= d)[1]] - d
  } else if (sum(x) < 100) {
    x[which(x < (100 - d))[1]] <- x[which(x < (100-d))[1]] + d
  }
  return(x)
}

trainFreq <- args$trainFreq
cellsMetadataFile <- args$cellsMetadataFile
cellTypeColumn <- args$cellTypeColumn
geneAnnotFile <- args$geneAnnotFile
simCountsFile <- args$simCountsFile
outputPath <- args$outputPath
probsDesignFile <- args$probsDesignFile
prefix <- args$prefix
#sCells <- args$sCells

dir.create(outputPath)
dir.create(file.path(outputPath,"Plots"),showWarnings = F)

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
cat(paste0(paste(c("Cells","Metadata"),dim(cellsMetadata))),"\n")
cat("Simulated Data")
table(cellsMetadata$Simulated,cellsMetadata[,cellTypeColumn])
cat("\n")



cat("Load analysisGeneAnnotations\n")
if (grepl(".tsv",geneAnnotFile)) {
  if (grepl(".tsv.gz",geneAnnotFile)) {
    cat("\tTab Gzipped format\n")
    analysisGeneAnnotations <- read.delim(file = gzfile(geneAnnotFile),sep = "\t",header = T,stringsAsFactors = F)
  } else {
    cat("\tTab plain format\n")
    analysisGeneAnnotations <- read.delim(file = geneAnnotFile,sep = "\t",header = T,stringsAsFactors = F)
  }
} else if (grepl(".rds",geneAnnotFile)) {
  cat("\tRDS R object\n")
  analysisGeneAnnotations <- readRDS(file = geneAnnotFile)
}
cat(paste0(paste(c("Genes","Metadata"),dim(analysisGeneAnnotations))),"\n")
head(analysisGeneAnnotations)
cat("\n")

cat("Load simCounts\n")
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
    quit(status = 1)
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

sCells <- dim(simCounts)[2]

# Split data into training and test sets
trainSet <- sample(colnames(simCounts),size = round(dim(simCounts)[2]*trainFreq))
trainTypes <- sub("_\\w+","",trainSet,perl = T)
trainSetList <- list()
for (ts in levels(factor(trainTypes))) {
  trainSetList[[ts]] <- trainSet[trainTypes == ts]
}
cat("Train Set cells by type\n")
tb <- unlist(lapply(trainSetList,length))
print(tb)
cat("\n")

testSet <- colnames(simCounts)[!colnames(simCounts) %in% trainSet]
testTypes <- sub("_\\w+","",testSet,perl = T)
testSetList <- list()
for (ts in levels(factor(testTypes))) {
  testSetList[[ts]] <- testSet[testTypes == ts]
}
cat("Test Set cells by type\n")
tb <- unlist(lapply(testSetList,length))
print(tb)
cat("\n")

probsDesign <- read.delim(probsDesignFile,sep = "\t",header = T)

probsLists <- apply(probsDesign,1,function (x) {return(seq(from = x['from'], to = x['to']))})

names(probsLists) <- probsDesign$cellType

jpeg(file.path(outputPath,"Plots",paste(prefix,"PredefineFreqRangeByCellType.jpg",sep = ".")))
df <- melt(probsLists)
colnames(df) <- c("Perc","CellType")
df$CellType <- factor(df$CellType,levels = names(probsLists))
ggplot(df,aes(x=CellType,y=Perc)) +
  geom_boxplot() + 
  ggtitle("Predefine Range of cell fractions")
dev.off()


### TRAIN SET 1
n <- ceiling(1000 * sCells/1000)
trainProbsMatrix <- matrix(rep(0,10),nrow=1,byrow = T)
while (dim(trainProbsMatrix)[1] < n+1) trainProbsMatrix <- rbind(trainProbsMatrix,unlist(lapply(probsLists,sample,1)))
trainProbsMatrix <- trainProbsMatrix[-1,]
trainProbsMatrix <- round(trainProbsMatrix*100/rowSums(trainProbsMatrix))

trainProbsMatrix <- t(apply(trainProbsMatrix,1,adjustHundred))

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set1Probs.Violin.jpg",sep = ".")))
df <- melt(trainProbsMatrix)
colnames(df) <- c("Sample","CellType","Prob")
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_violin() +
  ggtitle("Bulk Probility Dist. Set 1")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set1Probs.BoxPlot.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_boxplot() +
  ggtitle("Bulk Probility Dist. Set 1")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set1Probs.Lines.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob,group=Sample)) +
  geom_line(colour="grey60") +
  ggtitle("Bulk Probility Dist. Set 1")
dev.off()

dummy <- t(apply(trainProbsMatrix,1,sort,decreasing=T))
df <- melt(dummy)
colnames(df) <- c("Sample","MaxProb","Prob")
df$MaxProb <- factor(df$MaxProb)
jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set1Probs.Lines.Sorted.jpg",sep = ".")))
ggplot(df,aes(x=MaxProb,y=Prob)) +
  geom_boxplot()  +
  ggtitle("Bulk Probility Dist. Set 1")
dev.off()

### TRAIN SET 2
trainProbs <- list()
n <- ceiling(3000 * sCells/1000)
while(length(trainProbs) < n) {
  p <- rep(0,10)
  i <- 1
  while (sum(p) < 99) {
    p[i] <- p[i] + sample(seq(100 - sum(p)),size = 1)
    i <- i + 1
    if (i > 10) {
      i <- 1
    }
  }
  p[0] <- p[0] + 1
  if (sum(p==0) <= 9) {
    p <- sample(p)
    trainProbs[[length(trainProbs)+1]] <- p
  }
}

trainProbs <- matrix(unlist(trainProbs),nrow = n, byrow = T)
colnames(trainProbs) <- colnames(trainProbsMatrix)
trainProbs <- round(trainProbs*100/rowSums(trainProbs))

trainProbs <- t(apply(trainProbs,1,adjustHundred))

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set2Probs.Violin.jpg",sep = ".")))
df <- melt(trainProbs)
colnames(df) <- c("Sample","CellType","Prob")
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_violin() +
  ggtitle("Bulk Probility Dist. Set 2")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set2Probs.BoxPlot.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_boxplot() +
  ggtitle("Bulk Probility Dist. Set 2")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set2Probs.Lines.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob,group=Sample)) +
  geom_line(colour="grey60") +
  ggtitle("Bulk Probility Dist. Set 2")
dev.off()

dummy <- t(apply(trainProbs,1,sort,decreasing=T))
df <- melt(dummy)
colnames(df) <- c("Sample","MaxProb","Prob")
df$MaxProb <- factor(df$MaxProb)
jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set2Probs.Lines.Sorted.jpg",sep = ".")))
ggplot(df,aes(x=MaxProb,y=Prob)) +
  geom_boxplot()  +
  ggtitle("Bulk Probility Dist. Set 1")
dev.off()


trainProbsMatrix <- rbind(trainProbsMatrix,trainProbs)

### TRAIN SET 3
trainProbs <- list()
n <- ceiling(500 * sCells/1000)
while(length(trainProbs) < n) {
  p <- rep(0,10)
  names(p) <- names(probsLists)
  i <- 1
  while (sum(p) < 99) {
    dp <- 101
    while (dp > max(probsLists[[i]])) dp <- sample(probsLists[[i]],size = 1)
    p[i] <- dp
    i <- i+1
    if (i > 9) i <- 1
  }
  p[1] <- p[1] + 1
  if (sum(p==0) <= 7) {
    p <- sample(p)
    trainProbs[[length(trainProbs)+1]] <- p
  }
  # n <- n+1
}
trainProbs <- lapply(trainProbs,function (x) {return(x[names(probsLists)])})
trainProbs <- matrix(unlist(trainProbs),nrow = n, byrow = T)
colnames(trainProbs) <- colnames(trainProbsMatrix)
trainProbs <- round(trainProbs*100/rowSums(trainProbs))

trainProbs <- t(apply(trainProbs,1,adjustHundred))

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set3Probs.Violin.jpg",sep = ".")))
df <- melt(trainProbs)
colnames(df) <- c("Sample","CellType","Prob")
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_violin() +
  ggtitle("Bulk Probility Dist. Set 3")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set3Probs.BoxPlot.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_boxplot() +
  ggtitle("Bulk Probility Dist. Set 3")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set3Probs.Lines.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob,group=Sample)) +
  geom_line(colour="grey60") +
  ggtitle("Bulk Probility Dist. Set 3")
dev.off()

dummy <- t(apply(trainProbs,1,sort,decreasing=T))
df <- melt(dummy)
colnames(df) <- c("Sample","MaxProb","Prob")
df$MaxProb <- factor(df$MaxProb)
jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set3Probs.Lines.Sorted.jpg",sep = ".")))
ggplot(df,aes(x=MaxProb,y=Prob)) +
  geom_boxplot()  +
  ggtitle("Bulk Probility Dist. Set 1")
dev.off()

trainProbsMatrix <- rbind(trainProbsMatrix,trainProbs)

### TRAIN SET 4
trainProbs <- list()
n <- ceiling(1000 * sCells/1000)
while(length(trainProbs) < n) {
  p <- rep(0,10)
  names(p) <- names(probsLists)
  i <- 1
  while (sum(p) < 99) {
    dp <- sample(probsLists[[i]],size = 1)
    p[i] <- dp
    i <- i+1
    if (i > 9) i <- 1
  }
  p[1] <- p[1] + 1
  if (sum(p==0) <= 7) {
    p <- sample(p)
    trainProbs[[length(trainProbs)+1]] <- p
  }
  # n <- n+1
}
trainProbs <- lapply(trainProbs,sample)
trainProbs <- matrix(unlist(trainProbs),nrow = n, byrow = T)
colnames(trainProbs) <- colnames(trainProbsMatrix)
trainProbs <- round(trainProbs*100/rowSums(trainProbs))

trainProbs <- t(apply(trainProbs,1,adjustHundred))

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set4Probs.Violin.jpg",sep = ".")))
df <- melt(trainProbs)
colnames(df) <- c("Sample","CellType","Prob")
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_violin() +
  ggtitle("Bulk Probility Dist. Set 4")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set4Probs.BoxPlot.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_boxplot() +
  ggtitle("Bulk Probility Dist. Set 4")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set4Probs.Lines.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob,group=Sample)) +
  geom_line(colour="grey60") +
  ggtitle("Bulk Probility Dist. Set 4")
dev.off()

dummy <- t(apply(trainProbs,1,sort,decreasing=T))
df <- melt(dummy)
colnames(df) <- c("Sample","MaxProb","Prob")
df$MaxProb <- factor(df$MaxProb)
jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set4Probs.Lines.Sorted.jpg",sep = ".")))
ggplot(df,aes(x=MaxProb,y=Prob)) +
  geom_boxplot()  +
  ggtitle("Bulk Probility Dist. Set 1")
dev.off()


trainProbsMatrix <- rbind(trainProbsMatrix,trainProbs)

### TRAIN SET 5
n <- ceiling(3000 * sCells/1000)
d <- rdirichlet(n,rep(1,10))
d <- round(d*100)

d <- t(apply(d,1,adjustHundred))

colnames(d) <- colnames(trainProbsMatrix)

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set5Probs.Violin.jpg",sep = ".")))
df <- melt(trainProbs)
colnames(df) <- c("Sample","CellType","Prob")
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_violin() +
  ggtitle("Bulk Probility Dist. Set 5")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set5Probs.BoxPlot.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_boxplot() +
  ggtitle("Bulk Probility Dist. Set 5")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set5Probs.Lines.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob,group=Sample)) +
  geom_line(colour="grey60") +
  ggtitle("Bulk Probility Dist. Set 5")
dev.off()

dummy <- t(apply(d,1,sort,decreasing=T))
df <- melt(dummy)
colnames(df) <- c("Sample","MaxProb","Prob")
df$MaxProb <- factor(df$MaxProb)
jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set5Probs.Lines.Sorted.jpg",sep = ".")))
ggplot(df,aes(x=MaxProb,y=Prob)) +
  geom_boxplot()  +
  ggtitle("Bulk Probility Dist. Set 1")
dev.off()


trainProbsMatrix <- rbind(trainProbsMatrix,d)

### TRAIN SET 6
trainProbs <- list()
n <- ceiling(2000 * sCells/1000)
while (length(trainProbs) < n) trainProbs[[length(trainProbs)+1]] <- unlist(lapply(probsLists,sample,1))

trainProbs <- lapply(trainProbs,function(x){return(round(x*100/sum(x)))})

trainProbs <- lapply(trainProbs,sample)

trainProbs <- lapply(trainProbs,adjustHundred)

trainProbs <- matrix(unlist(trainProbs),nrow = n, byrow = T)
colnames(trainProbs) <- colnames(trainProbsMatrix)

# colnames(trainProbsMatrix) <- names(probsLists)

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set6Probs.Violin.jpg",sep = ".")))
df <- melt(trainProbs)
colnames(df) <- c("Sample","CellType","Prob")
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_violin() +
  ggtitle("Bulk Probility Dist. Set 6")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set6Probs.BoxPlot.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_boxplot() +
  ggtitle("Bulk Probility Dist. Set 6")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set6Probs.Lines.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob,group=Sample)) +
  geom_line(colour="grey60") +
  ggtitle("Bulk Probility Dist. Set 6")
dev.off()

dummy <- t(apply(trainProbs,1,sort,decreasing=T))
df <- melt(dummy)
colnames(df) <- c("Sample","MaxProb","Prob")
df$MaxProb <- factor(df$MaxProb)
jpeg(file.path(outputPath,"Plots",paste(prefix,"Trainning.Set6Probs.Lines.Sorted.jpg",sep = ".")))
ggplot(df,aes(x=MaxProb,y=Prob)) +
  geom_boxplot()  +
  ggtitle("Bulk Probility Dist. Set 1")
dev.off()


trainProbsMatrix <- rbind(trainProbsMatrix,trainProbs)

print(paste(c("bulks","types"),dim(trainProbsMatrix)))
head(trainProbsMatrix)


### TEST SET 1
n <- ceiling(700 * sCells/1000)
testProbsMatrix <- matrix(rep(0,10),nrow=1,byrow = T)
while (dim(testProbsMatrix)[1] < n+1) testProbsMatrix <- rbind(testProbsMatrix,unlist(lapply(probsLists,sample,1)))
testProbsMatrix <- testProbsMatrix[-1,]
testProbsMatrix <- round(testProbsMatrix*100/rowSums(testProbsMatrix))

testProbsMatrix <- t(apply(testProbsMatrix,1,adjustHundred))

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set1Probs.Violin.jpg",sep = ".")))
df <- melt(testProbsMatrix)
colnames(df) <- c("Sample","CellType","Prob")
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_violin() +
  ggtitle("Bulk Probility Dist. Set 1")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set1Probs.BoxPlot.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_boxplot() +
  ggtitle("Bulk Probility Dist. Set 1")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set1Probs.Lines.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob,group=Sample)) +
  geom_line(colour="grey60") +
  ggtitle("Bulk Probility Dist. Set 1")
dev.off()

### TEST SET 2
testProbs <- list()
n <- ceiling(2000 * sCells/1000)
while(length(testProbs) < n) {
  p <- rep(0,10)
  i <- 1
  while (sum(p) < 99) {
    p[i] <- p[i] + sample(seq(100 - sum(p)),size = 1)
    i <- i + 1
    if (i > 10) {
      i <- 1
    }
  }
  p[0] <- p[0] + 1
  if (sum(p==0) <= 9) {
    p <- sample(p)
    testProbs[[length(testProbs)+1]] <- p
  }
}

testProbs <- matrix(unlist(testProbs),nrow = n, byrow = T)
colnames(testProbs) <- colnames(testProbsMatrix)
testProbs <- round(testProbs*100/rowSums(testProbs))

testProbs <- t(apply(testProbs,1,adjustHundred))

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set2Probs.Violin.jpg",sep = ".")))
df <- melt(testProbs)
colnames(df) <- c("Sample","CellType","Prob")
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_violin() +
  ggtitle("Bulk Probility Dist. Set 2")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set2Probs.BoxPlot.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_boxplot() +
  ggtitle("Bulk Probility Dist. Set 2")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set2Probs.Lines.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob,group=Sample)) +
  geom_line(colour="grey60") +
  ggtitle("Bulk Probility Dist. Set 2")
dev.off()

testProbsMatrix <- rbind(testProbsMatrix,testProbs)

### TEST SET 3
testProbs <- list()
n <- ceiling(350 * sCells/1000)
while(length(testProbs) < n) {
  p <- rep(0,10)
  names(p) <- names(probsLists)
  i <- 1
  while (sum(p) < 99) {
    dp <- 101
    while (dp > max(probsLists[[i]])) dp <- sample(probsLists[[i]],size = 1)
    p[i] <- dp
    i <- i+1
    if (i > 9) i <- 1
  }
  p[1] <- p[1] + 1
  if (sum(p==0) <= 7) {
    p <- sample(p)
    testProbs[[length(testProbs)+1]] <- p
  }
  # n <- n+1
}
testProbs <- lapply(testProbs,function (x) {return(x[names(probsLists)])})
testProbs <- matrix(unlist(testProbs),nrow = n, byrow = T)
colnames(testProbs) <- colnames(testProbsMatrix)
testProbs <- round(testProbs*100/rowSums(testProbs))

testProbs <- t(apply(testProbs,1,adjustHundred))

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set3Probs.Violin.jpg",sep = ".")))
df <- melt(testProbs)
colnames(df) <- c("Sample","CellType","Prob")
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_violin() +
  ggtitle("Bulk Probility Dist. Set 3")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set3Probs.BoxPlot.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_boxplot() +
  ggtitle("Bulk Probility Dist. Set 3")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set3Probs.Lines.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob,group=Sample)) +
  geom_line(colour="grey60") +
  ggtitle("Bulk Probility Dist. Set 3")
dev.off()

testProbsMatrix <- rbind(testProbsMatrix,testProbs)

### TEST SET 4
testProbs <- list()
n <- ceiling(700 * sCells/1000)
while(length(testProbs) < n) {
  p <- rep(0,10)
  names(p) <- names(probsLists)
  i <- 1
  while (sum(p) < 99) {
    dp <- sample(probsLists[[i]],size = 1)
    p[i] <- dp
    i <- i+1
    if (i > 9) i <- 1
  }
  p[1] <- p[1] + 1
  if (sum(p==0) <= 7) {
    p <- sample(p)
    testProbs[[length(testProbs)+1]] <- p
  }
}
testProbs <- lapply(testProbs,sample)
testProbs <- matrix(unlist(testProbs),nrow = n, byrow = T)
colnames(testProbs) <- colnames(testProbsMatrix)
testProbs <- round(testProbs*100/rowSums(testProbs))

testProbs <- t(apply(testProbs,1,adjustHundred))

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set4Probs.Violin.jpg",sep = ".")))
df <- melt(testProbs)
colnames(df) <- c("Sample","CellType","Prob")
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_violin() +
  ggtitle("Bulk Probility Dist. Set 4")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set4Probs.BoxPlot.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_boxplot() +
  ggtitle("Bulk Probility Dist. Set 4")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set4Probs.Lines.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob,group=Sample)) +
  geom_line(colour="grey60") +
  ggtitle("Bulk Probility Dist. Set 4")
dev.off()

testProbsMatrix <- rbind(testProbsMatrix,testProbs)

### TEST SET 5
n <- ceiling(2000 * sCells/1000)
d <- rdirichlet(n,rep(1,10))
d <- round(d*100)
d <- t(apply(d,1,adjustHundred))

colnames(d) <- colnames(testProbsMatrix)

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set5Probs.Violin.jpg",sep = ".")))
df <- melt(testProbs)
colnames(df) <- c("Sample","CellType","Prob")
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_violin() +
  ggtitle("Bulk Probility Dist. Set 5")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set5Probs.BoxPlot.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_boxplot() +
  ggtitle("Bulk Probility Dist. Set 5")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set5Probs.Lines.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob,group=Sample)) +
  geom_line(colour="grey60") +
  ggtitle("Bulk Probility Dist. Set 5")
dev.off()

testProbsMatrix <- rbind(testProbsMatrix,d)

### TEST SET 6
testProbs <- list()
n <- ceiling(1400 * sCells/1000)
while (length(testProbs) < n) testProbs[[length(testProbs)+1]] <- unlist(lapply(probsLists,sample,1))

testProbs <- lapply(testProbs,function(x){return(round(x*100/sum(x)))})

testProbs <- lapply(testProbs,sample)

testProbs <- lapply(testProbs, adjustHundred)

testProbs <- matrix(unlist(testProbs),nrow = n, byrow = T)
colnames(testProbs) <- colnames(testProbsMatrix)

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set6Probs.Violin.jpg",sep = ".")))
df <- melt(testProbs)
colnames(df) <- c("Sample","CellType","Prob")
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_violin() +
  ggtitle("Bulk Probility Dist. Set 6")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set6Probs.BoxPlot.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob)) +
  geom_boxplot() +
  ggtitle("Bulk Probility Dist. Set 6")
dev.off()

jpeg(file.path(outputPath,"Plots",paste(prefix,"Test.Set6Probs.Lines.jpg",sep = ".")))
ggplot(df,aes(x=CellType,y=Prob,group=Sample)) +
  geom_line(colour="grey60") +
  ggtitle("Bulk Probility Dist. Set 6")
dev.off()

testProbsMatrix <- rbind(testProbsMatrix,testProbs)

print(paste(c("bulks","types"),dim(testProbsMatrix)))
head(testProbsMatrix)


### GENERATE PROBS MATRIX NAMES
cat(paste(names(trainSetList)),"\n")

setCount <- function (x, setList, sn) {
  names(x) <- sn
  sc <- c() 
  for (cType in names(x)) {
    n <- ceiling(x[cType])
    if (n > 0) {
      # if (cType == "Tumour") {cType <- sample(c("HER2","luminalA","luminalB","TNBC"),1)}
      repl <- ifelse(n > length(setList[[cType]]),TRUE,FALSE)
      sc <- c(sc,sample(setList[[cType]],size = n,replace = repl))  
    }
  }
  return(sc[seq(100)])
}

trainProbMatrixNames <- apply(trainProbsMatrix,1,setCount,setList=trainSetList,sn = colnames(trainProbsMatrix))

trainProbMatrixNames <- t(trainProbMatrixNames)


testProbMatrixNames <- apply(testProbsMatrix,1,setCount,setList=testSetList,sn = colnames(testProbsMatrix))

testProbMatrixNames <- t(testProbMatrixNames)

cat("WRITE OUTPUTS\n")

saveRDS(trainProbsMatrix,file = file.path(outputPath,paste(prefix,"trainProbsMatrix","rds",sep=".")))
saveRDS(trainProbMatrixNames,file = file.path(outputPath,paste(prefix,"trainProbMatrixNames","rds",sep=".")))
saveRDS(trainSetList,file = file.path(outputPath,paste(prefix,"trainSetList","rds",sep=".")))
saveRDS(trainSet,file = file.path(outputPath,paste(prefix,"trainSet","rds",sep=".")))

saveRDS(testProbsMatrix,file = file.path(outputPath,paste(prefix,"testProbsMatrix","rds",sep=".")))
saveRDS(testProbMatrixNames,file = file.path(outputPath,paste(prefix,"testProbMatrixNames","rds",sep=".")))
saveRDS(testSetList,file = file.path(outputPath,paste(prefix,"testSetList","rds",sep=".")))
saveRDS(testSet,file = file.path(outputPath,paste(prefix,"testSet","rds",sep=".")))

cat("DONE\n")
