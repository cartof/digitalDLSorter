# digitalDLSorter
A pipeline to generate a Deep Neural Network cell type deconvolution model for bulk RNASeq samples from single cell rna-seq data.

# Installation
The pipeline runs in R 3.51 and python 3.6.

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

Under the data folder there are counts, cellsMetada and genesMetada files to test the pipeline



