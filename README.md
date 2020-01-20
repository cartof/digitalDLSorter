# digitalDLSorter
A pipeline to generate a Deep Neural Network cell type deconvolution model for bulk RNASeq samples from single cell rna-seq data. As described in [Torroja and Sanchez-Cabo, 2019](https://www.frontiersin.org/articles/10.3389/fgene.2019.00978/full). 

# Installation
The pipeline runs in R 3.5.1 and python 3.6.

Install R with the following packages: 
- splatter
- zimbwave
- SingleCellExperiment
- Matrix
- dplyr
- gtools
- ggplot2
- resahpe2
- edgeR
- pbapply
- optparse

Install Python 3 with the following packages:
- numpy
- pandas
- argparse
- matplotlib
- keras

Using conda enviroment:
```
conda create -n digitalDLSorter
conda activate digitalDLSorter
conda install tensorflow=1.5.0 python=3.6 numpy pandas argparse matplotlib keras
conda install r=3.5.1 bioconductor-splatter bioconductor-zinbwave bioconductor-SingleCellExperiment r-gtools r-dplyr r-ggplot2 r-reshape2  bioconductor-edgeR r-pbapply r-optparse
```

# Pipeline summary

![Pipeline Graph](PipelineScript.png)

# Running the pipeline

Each script has its own help instructions.
```
Rscript estimateZinbwaveParams.R --help
python ~/Documents/GitLabProjects/uploadedDLSorter/digitalDLSorter.py --help
```
## Input data
Under the data folder there are counts, cellsMetada, genesMetada and probPriors files to test the pipeline. These data has been obtained from [GSE81861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81861) ([Li et al., 2017](https://doi.org/10.1038/ng.3818)), a collection of single cell RNA-Seq profiles from 11 CRC patients.

There is also a a set of bulk RNA-Seq profiles from TCGA colorectal samples obtained from [TCGA project data portal](https://portal.gdc.cancer.gov/).

countsMatrix.tsv.gz: File with single cell counts. Genes in rows and Cells in columns.

```bash
zcat data/countsMatrix.tsv.gz | cut -f1-5 | head -n6
```
RHC3934 | RHC3944 | RHC3962 | RHC4005 | RHC4220
------- | ------- | ------- | ------- | -------
7SK | 0 | 0 | 0 | 0
A1BG | 0 | 0 | 0 | 0
A1BG-AS1 | 0 | 0 | 0 | 0
A1CF | 0 | 0 | 0.075 | 0
A2M | 0 | 0 | 0 | 0

cellsMetadata.tsv.gz: file with annotations for each cell (rows) in countsMatrix
```bash
zcat cellsMetadata.tsv.gz | cut -f1-4,12 | head -n6
```
Cell_ID | Patient | issue_type | gender | Cell_Type
------- | ------- | ---------- | ------ | ---------
RHC3934 | CRC03 | NM | Female | pB
RHC3944 | CRC03 | NM | Female | gB
RHC3962 | CRC03 | NM | Female | CD8Gn
RHC4005 | CRC03 | NM | Female | gB
RHC4220 | CRC04 | NM | Female | pB

genesMetadata.tsv.gz: file with annotations for each gene (rows) in countsMatrix
```
zcat data/genesMetadata.tsv.gz | head -n 6
```
ensembl_gene_id | external_gene_name | gene_length
--------------- | ------------------ | -----------
ENSG00000000003 | TSPAN6 | 2959
ENSG00000000005 | TNMD | 1603
ENSG00000000419 | DPM1 | 1197
ENSG00000000457 | SCYL3 | 5561
ENSG00000000460 | C1ORF112 | 4936

probPriors.txt: file with the frequency ranges expected for each cell type annotated in cellsMatadata. This information can be estimated from literature or from the single cell experiment itself. It is a priory guess to build the simulated bulk profiles.
```
head -n6 data/probPriors.txt
```
cellType | from | to
-------- | ---- | --
CRC | 30 | 70
CD4 | 0 | 15
CD8Gn | 0 | 15
CD8Gp | 0 | 15
Ep | 1 | 50

## Simulating new Single Cell profiles

In cases where there is low numbers of cells in the experiment or there is a particular cell type which is under represented we can increase the sample size by simulating new data based on the real one. This pipeline uses the ZinbWave framework implemented in R to estimate the parameters to simulate single cell profiles ([Risso et al., 2018](https://www.frontiersin.org/articles/10.3389/fgene.2019.00978/full#B36)).

### Estimate simulation parameters

We can adjust the cell type model `--cellTypeColumn=Cell_Type` based on our experimental design adding covariates such as patient and gender for example `--cellCovColumns=Patient,gender`. Also we can add covariates at gene level, for example the gene length `--geneCovColumns=gene_length`.

We do also filter genes based on minimum expression at a minimum number of cells `--minCounts=1 --minCells=1`.

```bash
conda activate digitalDLSorter
Rscript estimateZinbwaveParams.R \
-t All \
-m data/cellsMetadata.tsv.gz --cellTypeColumn=Cell_Type --cellIDColumn=Cell_ID --cellCovColumns=Patient,gender \
-g data/genesMetadata.tsv.gz --geneIDColumn=external_gene_name --geneCovColumns=gene_length \
-c data/countsMatrix.tsv.gz --minCounts=1 --minCells=1 \
-o digitalDLSorterDemo -p Demo -d 4
```
Be patient this may take a few hours to run. It is not very optimized. `-d 4` will use 4 threads for some steps of the estimation. Adjust it depending on your resources.
Results will be stored at `digitalDLSorterDemo/` folder.

### Simulating Single Cells

Using the `Demo.All.zinbParams.rds` r object file with the estimated parameters we will generate new single cell profiles.

We will simulate 100 cell profiles per cell type `-t All -n 100`. In our dataset we have 10 cells types so we will end up with 1000 simulated single cell profiles.

```bash
conda activate digitalDLSorter
Rscript simSingleCellProfiles.R \
-t All -n 100 \
-z digitalDLSorterDemo/Demo.All.zinbParams.rds \
-c data/cellsMetadata.tsv.gz -q Cell_Type \
-m data/countsMatrix.tsv.gz \
-o digitalDLSorterDemo -p Demo
```
 The script generates a counts matrix stored as a sparse matrix in this folder `digitalDLSorterDemo/Demo.All.simCellsCountsMtx.1000` and a metadata file for the simulated cells in two formats, a gzip tsv file `Demo.All.simCellsMetadata.1000.tsv.gz` and an R object `Demo.All.simCellsMetadata.1000.rds`. This files contains both , the original single cell profiles and the simulated ones. 

## Simulate Bulk profiles

To train the digitalDLSorter deconvolution model we will first generate bulk samples with known proportions of cells.

### Generate vectors of cell type proportions for training and test set

First we will separate the cells into training and test set and create a list of vectors with the cell type fractions for each bulk sample.

We will split the Single Cell profiles in 65% for training and 35% for testing with `-f 0.65`. The script will need the simulated metadata file in tsv or r object format `-c digitalDLSorterDemo/Demo.All.simCellsMetadata.1000.rds`, the column that contains the target classification we want to train on `-q Cell_Type` and its corresponding counts matrix `-s digitalDLSorterDemo/Demo.All.simCellsCountsMtx.1000/matrix.mtx` from the previous step.

```bash
conda activate digitalDLSorter
Rscript generateTrainAndTestBulkProbMatrix.R \
-f 0.65 -g  data/genesMetadata.tsv.gz -d data/probPriors.txt \
-c digitalDLSorterDemo/Demo.All.simCellsMetadata.1000.rds -q Cell_Type \
-s digitalDLSorterDemo/Demo.All.simCellsCountsMtx.1000/matrix.mtx \
-o digitalDLSorterDemo -p Demo
```

This script will generate several files separated in Train and Test sets.

`Demo.trainProbsMatrix` and `Demo.testProbMatrix` with the list of vectors with cell type fractions to simulate. Cell types in columns and simulated sample in rows.

`Demo.trainProbMatrixNames` and `Demo.testProbMatrixNames` with the names of the cells that have been sampled from the training and test pools repectively to build the bulk samples. Each sample in rows made of a 100 cells.

`Demo.trainSetList` and `Demo.testSetList` files containing the list of cell types with the pool of cells belonging to each cell type

### Generate Train and Test Bulk Samples

Using the vectors of cell type fractions created for Train and Test and the Single Cell matrix counts, we will generate the counts profiles for the simulated Bulk samples.

We will run this step twice, one for the training set and one for the test set.

Training Set
```bash
conda activate digitalDLSorter
Rscript generateBulkSamples.R \
-m digitalDLSorterDemo/Demo.trainProbMatrixNames.rds \
-s digitalDLSorterDemo/Demo.All.simCellsCountsMtx.1000/matrix.mtx \
-o digitalDLSorterDemo -p Demo.Train -n 6
```

Test set
```bash
conda activate digitalDLSorter
Rscript generateBulkSamples.R \
-m digitalDLSorterDemo/Demo.testProbMatrixNames.rds \
-s digitalDLSorterDemo/Demo.All.simCellsCountsMtx.1000/matrix.mtx \
-o digitalDLSorterDemo -p Demo.Test -n 6
```

This step will generate the `Demo.Train.simBulkCounts and `Demo.Test.simBulkCounts` matrix counts in tsv.gz and r .rds object. With genes in rows and samples in columns.

## Prepare Single Cell And Bulk profiles for training

Now we will prepare the data to feed the digitalDLSorter model by combining the simulated bulk and single cell profiles, normalize, shuffle and scale them.

We will need the `Demo.trainProbsMatrix` or `Demo.testProbsMatrix` file with the vectors cell type fractions, the `Demo.trainSetList` or `Demo.testSetList` file with the list of cells selected for each set and the counts matrix for the single cell profiles `Demo.All.simCellsCountsMtx.1000/matrix.mtx` and the bulk profiles `Demo.Train.simBulkCounts` or `Demo.Test.simBulkCounts`.

Again, we will run this step twice, one for the training set and one for the test set.

Training Set
```bash
conda activate digitalDLSorter
Rscript generateCombinedScaledShuffledSet.R \
-m digitalDLSorterDemo/Demo.trainProbsMatrix.rds -c digitalDLSorterDemo/Demo.trainSetList.rds \
-s digitalDLSorterDemo/Demo.All.simCellsCountsMtx.1000/matrix.mtx -b digitalDLSorterDemo/Demo.Train.simBulkCounts.rds \
-o digitalDLSorterDemo -p Demo.Train
```

Test Set
```bash
conda activate digitalDLSorter
Rscript generateCombinedScaledShuffledSet.R \
-m digitalDLSorterDemo/Demo.testProbsMatrix.rds -c digitalDLSorterDemo/Demo.testSetList.rds \
-s digitalDLSorterDemo/Demo.All.simCellsCountsMtx.1000/matrix.mtx -b digitalDLSorterDemo/Demo.Test.simBulkCounts.rds \
-o digitalDLSorterDemo -p Demo.Test
```

This step will generate two versions of the data, one with all combined single cell and bulk samples and counts without further transformation (`combinedCounts` and  `combinedProbsMatrix`) and one with the samples normalized, scaled and shuffled (`combinedCounts.log2CPMScaledShuffledTransposed` and `combinedProbsMatrix.log2CPMScaledShuffledTransposed`). It also provides the list of genes used `combinedCounts.geneList`.


## Train digitalDLSorter model

With the `log2CPMScaledShuffledTransposed` files of data (samples in rows and genes in columns) we will train the model.

We need to provide the number of training `--num_train_samples 16056` and test samples `--num_test_samples 10800` and the final number of genes used after filtering steps `--num_genes 36477`.

You could use this bash lines to find out the numbers
```bash
zcat digitalDLSorterDemo/Demo.Train.combinedCounts.log2CPMScaledShuffledTransposed.tsv.gz | wc -l | awk '{print $0-1}'

zcat digitalDLSorterDemo/Demo.Test.combinedCounts.log2CPMScaledShuffledTransposed.tsv.gz | wc -l | awk '{print $0-1}'

zcat digitalDLSorterDemo/Demo.Train.combinedCounts.geneList.tsv.gz | wc -l
```

The number of cell types `--num_classes 10`,  batch size to process the data in chunks `--batch_size 100 ` and number of rounds for training `--num_epochs 20`.

Finaly we provide the loss function to reduce `--loss kullback_leibler_divergence` for this demo.

```bash
conda activate digitalDLSorter
python digitalDLSorter.py \
--trainCountsFile digitalDLSorterDemo/Demo.Train.combinedCounts.log2CPMScaledShuffledTransposed.tsv.gz \
--trainProbsFile digitalDLSorterDemo/Demo.Train.combinedProbsMatrix.log2CPMScaledShuffledTransposed.tsv.gz \
--testCountsFile digitalDLSorterDemo/Demo.Test.combinedCounts.log2CPMScaledShuffledTransposed.tsv.gz \
--testProbsFile digitalDLSorterDemo/Demo.Test.combinedProbsMatrix.log2CPMScaledShuffledTransposed.tsv.gz \
--num_classes 10 --num_genes 36477 --num_train_samples 16056 --num_test_samples 10800 \
--batch_size 100 --num_epochs 10 --loss kullback_leibler_divergence \
--outputPath digitalDLSorterDemo --prefix Demo
```

The training step will generate an h5 file with the model `Demo.digitalDLSorterTrainedModel.kullback_leibler_divergence.h5` and its weights `Demo.digitalDLSorterTrainedWeightsModel.kullback_leibler_divergence`, the list of genes for the input `Demo.digitalDLSorterTrainedModel.inputGeneList`, the target class names or cell types `Demo.digitalDLSorterTrainedModel.targetClassNames`, the model performance evaluation `Demo.digitalDLSorterTrainedModel.evalStats` and the predictions on the test set `Demo.digitalDLSorterTrainedModel.DeconvTestPredictions`.

## Predict samples with the model

Now with the model trained we can provide new samples to be deconvoluted. In example, these colorectal samples obtained from TCGA site. The counts matrix of bulk samples should have samples in rows and genes in columns. If the data is not normalized and scaled we can do so by `--normData True`. If we have trained the model starting from TPMs we should feed it with TPMs.

```bash
conda activate digitalDLSorter
python digitalDLSorterModelDeconv.py \
--modelFile digitalDLSorterDemo/Demo.digitalDLSorterTrainedModel.kullback_leibler_divergence.h5 \
--modelGenesList digitalDLSorterDemo/Demo.digitalDLSorterTrainedModel.inputGeneList.kullback_leibler_divergence.txt.gz \
--modelClassNames digitalDLSorterDemo/Demo.digitalDLSorterTrainedModel.targetClassNames.kullback_leibler_divergence.txt.gz \
--countsFile data/TCGA.Colon.Counts.transpose.tsv.gz \
--batch_size 100 --num_samples 521 --normData \
--outputPath digitalDLSorterDemo --prefix Demo.PredictTCGA
```
