import keras as keras
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import threading
from keras import metrics
from keras.models import load_model
import sys
import argparse

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
                  help='file prefix to write the model to')
args = parser.parse_args()

trainCountsFile = args.trainCountsFile
trainProbsFile = args.trainProbsFile
testCountsFile = args.testCountsFile
testProbsFile = args.testProbsFile
outputPath = args.outputPath

num_classes = args.num_classes
num_genes = args.num_genes
batch_size = args.batch_size
num_train_samples = args.num_train_samples
num_test_samples = args.num_test_samples
num_epochs = args.num_epochs

def trainIter (countsFile,probsFile,size):
	sys.stderr.write('generator initiated\n')
	nb = 1
	while 1:
		counts = pd.read_csv(countsFile,compression='gzip',sep='\t',chunksize=size)
		probs = pd.read_csv(probsFile,compression='gzip',sep='\t',chunksize=size)
		for x,y in zip(counts,probs):
			countsIndex = np.arange(x.shape[0])
			np.random.shuffle(countsIndex)
			x = x.iloc[countsIndex,:]
			y = y.iloc[countsIndex,:]
			sys.stderr.write('\ngenerator yielded a batch %d \n' % nb)
			yield(x,y)
			nb += 1

def predictIter (countsFile,size):
	sys.stderr.write('generator initiated\n')
	nb = 1
	while 1:
		counts = pd.read_csv(countsFile,compression='gzip',sep='\t',chunksize=size)
		for x in counts:
			sys.stderr.write('\ngenerator yielded a batch %d \n' % nb)
			yield(x)
			nb += 1

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

fh = open(outputPath+".digitalDLSorterModel."+loss+".json","w")
fh.writelines(json_string)
fh.close()

print("Compile model")
model.compile(loss=loss, optimizer='adam', metrics=[metrics.mae,metrics.mean_absolute_percentage_error,metrics.kullback_leibler_divergence,metrics.categorical_accuracy])

print("")

train = trainIter(trainCountsFile,trainProbsFile,batch_size)

print("Train model")
print("Steps per epoch"+str(np.ceil(float(num_train_samples)/float(batch_size))))
hist = model.fit_generator(train
, epochs=num_epochs
, verbose = 1
, steps_per_epoch=np.ceil(float(num_train_samples)/float(batch_size))
, use_multiprocessing=True
)

model.save(outputPath+".digitalDLSorterTrainedModel."+loss+".h5")

model.save_weights(outputPath+'.digitalDLSorterTrainedWeightsModel.'+loss+'.h5')

print(hist.history)
print("")

print("Evaluate model")
print("Steps per epoch"+str(np.ceil(float(num_test_samples)/float(batch_size))))

test = trainIter(testCountsFile,testProbsFile,batch_size)

testEval = model.evaluate_generator(test,steps=np.ceil(float(num_test_samples)/float(batch_size)))

np.savetxt(outputPath+".digitalDLSorterTrainedModel.evalStats."+loss+".txt",testEval,delimiter='\t',newline='\n')

print("Predict Tests")

test = predictIter(testCountsFile,batch_size)

testPredictions = model.predict_generator(test,steps=np.ceil(num_test_samples/batch_size),verbose=1)

np.savetxt(outputPath+".digitalDLSorterTrainedModel.DeconvTestPredictions."+loss+".txt",testPredictions,delimiter='\t',newline='\n')
print("DONE")
