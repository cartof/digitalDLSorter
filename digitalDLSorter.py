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


class dataGenerator(keras.utils.Sequence):
	def __init__(self, data_file, probs_file, num_samples, batch_size, shuffle=True):
		sys.stderr.write('generator initiated\n')
		self.batch_size = batch_size
		self.data_file = data_file
		self.probs_file = probs_file
		self.num_samples = num_samples
		self.shuffle = shuffle
		self.data = pd.read_csv(data_file,compression='gzip',sep='\t',chunksize=batch_size)
		self.probs = pd.read_csv(probs_file,compression='gzip',sep='\t',chunksize=batch_size)
		self.nb = 0
	def __len__(self):
		return int(np.ceil(float(self.num_samples)/float(self.batch_size)))
	def __getitem__(self, index):
		while 1:
			if self.nb == self.__len__():
				self.data.close()
				self.probs.close()
				self.data = pd.read_csv(self.data_file,compression='gzip',sep='\t',chunksize=self.batch_size)
				self.probs = pd.read_csv(self.probs_file,compression='gzip',sep='\t',chunksize=self.batch_size)
				self.nb = 0
			for X,y in zip(self.data,self.probs):
				if self.shuffle == True:
					dataIndex = np.arange(X.shape[0])
					np.random.shuffle(dataIndex)
					X = X.iloc[dataIndex,:]
					y = y.iloc[dataIndex,:]
				self.nb += 1
				sys.stderr.write('\ngenerator yielded a batch %d \n' % self.nb)
				sys.stderr.write('size: %d \n' % X.shape[0])
				sys.stderr.write('\n')
				return X, y

class predictDataGenerator(keras.utils.Sequence):
	def __init__(self, data_file, num_samples, batch_size):
		sys.stderr.write('generator initiated\n')
		self.batch_size = batch_size
		self.data_file = data_file
		self.num_samples = num_samples
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
			return X
	def get_sample_names(self):
		return list(itertools.chain.from_iterable(self.sample_names))


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--trainCountsFile', metavar='trainCountsFile',
                   help='File containing counts matrix separated by tabs with header (cell/bulk names) and rownames (genes)')
parser.add_argument('--trainProbsFile', metavar='trainProbsFile',
                  help='File containing probs matrix separated by tabs with header (cell types) and rownames (cell/bulk names)')
parser.add_argument('--testCountsFile', metavar='testCountsFile',
                   help='File containing counts matrix separated by tabs with header (cell/bulk names) and rownames (genes)')
parser.add_argument('--testProbsFile', metavar='testProbsFile',
                  help='File containing probs matrix separated by tabs with header (cell types) and rownames (cell/bulk names)')
parser.add_argument('--num_classes', metavar='num_classes', type=int,
                  help='Number of cell types defined')
parser.add_argument('--num_genes', metavar='num_genes', type=int,
                  help='Number of genes in the counts matrix')
parser.add_argument('--batch_size', metavar='batch_size', type=int, default=100,
                  help='Number of samples per batch')
parser.add_argument('--num_train_samples', metavar='num_train_samples', type=int,
                  help='Number of samples in train counts matrix')
parser.add_argument('--num_test_samples', metavar='num_test_samples', type=int,
                  help='Number of samples in test counts matrix')
parser.add_argument('--num_epochs', metavar='num_epochs', type=int, default=50,
                  help='Number of epochs during trainning')
parser.add_argument('--loss', metavar='loss', default='kullback_leibler_divergence',
                  help='loss function for trainning (mean_absolute_error,mean_absolute_percentage_error,categorical_crossentropy,kullback_leibler_divergence)')
parser.add_argument('--outputPath', metavar='outputPath',
                  help='output path to write the model to')
parser.add_argument('--prefix', metavar='prefix',
                  help='file prefix to write the model to')

args = parser.parse_args()

trainCountsFile = args.trainCountsFile
trainProbsFile = args.trainProbsFile
testCountsFile = args.testCountsFile
testProbsFile = args.testProbsFile
outputPath = args.outputPath
prefix = args.prefix

num_classes = args.num_classes
num_genes = args.num_genes
batch_size = args.batch_size
num_train_samples = args.num_train_samples
num_test_samples = args.num_test_samples
num_epochs = args.num_epochs
loss = args.loss

d = pd.read_csv(trainCountsFile,nrows=2,compression='gzip',sep='\t')
inputGeneList = d.columns.get_values()
np.savetxt(outputPath+"/"+prefix+".digitalDLSorterTrainedModel.inputGeneList."+loss+".txt.gz",inputGeneList,fmt='%s')

d = pd.read_csv(trainProbsFile,nrows=2,compression='gzip',sep='\t')
targetClassNames = d.columns.get_values()
np.savetxt(outputPath+"/"+prefix+".digitalDLSorterTrainedModel.targetClassNames."+loss+".txt.gz",targetClassNames,fmt='%s')

predictSampleNames = []

print("Build model")
model = keras.models.Sequential()
model.add(keras.layers.Dense(200, input_dim=num_genes, name="Dense1"))
model.add(keras.layers.BatchNormalization())
model.add(keras.layers.Activation('relu'))
model.add(keras.layers.Dropout(0.25))
model.add(keras.layers.Dense(200, name="Dense2"))
model.add(keras.layers.BatchNormalization())
model.add(keras.layers.Activation('relu'))
model.add(keras.layers.Dropout(0.25))
model.add(keras.layers.Dense(num_classes, name="Dense3"))
model.add(keras.layers.BatchNormalization())
model.add(keras.layers.Activation('softmax'))

model.summary()

json_string = model.to_json()

fh = open(outputPath+"/"+prefix+".digitalDLSorterModel."+loss+".json","w")
fh.writelines(json_string)
fh.close()

print("Compile model")
model.compile(loss=loss, optimizer='adam', metrics=[metrics.mae,metrics.mean_absolute_percentage_error,metrics.kullback_leibler_divergence,metrics.categorical_accuracy])

print("")

print("Train model")
train = dataGenerator(trainCountsFile, trainProbsFile, num_train_samples, batch_size, shuffle=True)

print("Train with "+str(train.num_samples)+" samples in "+str(train.__len__())+" Steps per epoch")
hist = model.fit_generator(train
, epochs=num_epochs
, verbose = 1
, steps_per_epoch=train.__len__()
, use_multiprocessing=False
)

model.save(outputPath+"/"+prefix+".digitalDLSorterTrainedModel."+loss+".h5")

model.save_weights(outputPath+"/"+prefix+'.digitalDLSorterTrainedWeightsModel.'+loss+'.h5')

#fh = open(outputPath+"/"+prefix+".digitalDLSorterTrainedModel.History."+loss+".json","w")
#fh.writelines(hist.history)
#fh.close()

print(hist.history)
print("")

print("Evaluate model")
test = dataGenerator(testCountsFile,testProbsFile,num_test_samples,batch_size, shuffle=False)

print("Evaluate "+str(test.num_samples)+" Test samples in "+str(test.__len__())+" steps\n")
testEval = model.evaluate_generator(test
, steps=test.__len__()
, use_multiprocessing=False
)

np.savetxt(outputPath+"/"+prefix+".digitalDLSorterTrainedModel.evalStats."+loss+".txt.gz",testEval,delimiter='\t',newline='\n')

print("Predict Tests")
test = predictDataGenerator(testCountsFile,num_test_samples,batch_size)

print("Samples to predict "+str(test.num_samples)+" in "+str(test.__len__())+" steps\n")
testPredictions = model.predict_generator(test,steps=np.ceil(num_test_samples/batch_size),verbose=1)

predictSampleNames = test.get_sample_names()

testPredictions = pd.DataFrame(testPredictions,index=predictSampleNames,columns=targetClassNames)

testPredictions.to_csv(outputPath+"/"+prefix+".digitalDLSorterTrainedModel.DeconvTestPredictions."+loss+".txt.gz",compression="gzip",sep='\t')

print("DONE")
