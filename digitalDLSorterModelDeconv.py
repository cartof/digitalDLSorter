import keras as keras
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import threading
from keras import metrics
from keras.models import load_model
import sys
import argparse
import itertools
from sklearn.preprocessing import scale

class predictIter(keras.utils.Sequence):
	def __init__(self, data_file, num_samples, batch_size, inputGeneNames, normData=False):
		sys.stderr.write('generator initiated\n')
		self.batch_size = batch_size
		self.data_file = data_file
		self.num_samples = num_samples
		self.normData = normData
		self.inputGeneNames = inputGeneNames
		self.data = pd.read_csv(data_file,compression='gzip',sep='\t',chunksize=batch_size)
		self.nb = 0
		self.sample_names = []
	def __len__(self):
		return int(np.ceil(float(self.num_samples)/float(self.batch_size)))
	def __getitem__(self, index):
		for X in self.data:
			self.nb += 1
			self.sample_names.append(X.index.get_values())
			sys.stderr.write('\ngenerator yielded a batch %d \n' % self.nb)
			sys.stderr.write('size: %d \n' % X.shape[0])
			X = X.transpose().reindex(self.inputGeneNames).fillna(0).transpose()
			if self.normData:
				X = self.cpm(X,doLog=True)
				sc = scale(X, axis=1, with_mean=True, with_std=True)
				X = pd.DataFrame(sc,index=X.index.get_values(),columns=X.columns.get_values())
			return X
	def get_samples_names(self):
		return list(itertools.chain.from_iterable(self.sample_names))
	def cpm(self,x,doLog=True):
		x = x.div(x.sum(axis=1), axis="rows").multiply(1000000)
		if doLog:
			x = x.add(1).transform(np.log2)
		return x

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--modelFile', metavar='modelFile',
                   help='File containing the trained digitalDLSorter model (h5 format)')
parser.add_argument('--modelGenesList', metavar='modelGenesList',
                   help='File containing the list of genes used to train the model')
parser.add_argument('--modelClassNames', metavar='modelClassNames',
                   help='File containing the names of the cell types used to train the model')   
parser.add_argument('--countsFile', metavar='countsFile',
                  help='File containing counts matrix separated by tabs WITH GENES IN COLUMNS and SAMPLES IN ROWS')
parser.add_argument('--batch_size', metavar='batch_size', type=int, default=100,
                  help='Number of samples per batch')
parser.add_argument('--num_samples', metavar='num_samples', type=int,
                  help='Number of samples in counts matrix')
parser.add_argument('--normData', metavar='normData', type=bool, default=False,
                  help='Whether to log normalize and scale the data if not done before')
parser.add_argument('--outputPath', metavar='outputPath',
                  help='output path to write the results to')
parser.add_argument('--prefix', metavar='prefix',
                  help='file prefix to write the results to')

args = parser.parse_args()

outputPath = args.outputPath
prefix = args.prefix

modelFile = args.modelFile
genesListFile = args.modelGenesList
classNamesFile = args.modelClassNames
countsFile = args.countsFile

num_samples = args.num_samples
batch_size = args.batch_size
normData = args.normData

model = keras.models.load_model(modelFile)
geneList = pd.read_csv(genesListFile,compression='gzip',sep='\t',names=["inputGeneNames"])
classNames = pd.read_csv(classNamesFile,compression='gzip',sep='\t',names=["targetClassNames"])

print("Predict Cell Fractions\n")
samplesIter = predictIter(countsFile,num_samples,batch_size,geneList.inputGeneNames,normData=normData)

print("Sample Set:"+str(samplesIter.num_samples)+" Steps:"+str(samplesIter.__len__())+"\n")
cellFractionPredictions = model.predict_generator(samplesIter,steps=samplesIter.__len__(),verbose=1, use_multiprocessing=False)

print("Save Results\n")
predictSampleNames = samplesIter.get_samples_names()

cellFractionPredictions = pd.DataFrame(cellFractionPredictions,index=predictSampleNames,columns=classNames.targetClassNames)

cellFractionPredictions.to_csv(outputPath+"/"+prefix+".CellFractionPredictions.txt.gz",compression="gzip",sep='\t')
print("DONE\n")

