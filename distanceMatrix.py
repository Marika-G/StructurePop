from datetime import datetime
import functools as fct
import multiprocessing as mp
import numpy as np
import pandas as pd
import scipy.spatial.distance as sciDist
import scipy.stats as sciStt
import sys


'''
	:param 1 : Path to data
	:oaram 2 : Boolean, normalize or not? Normalize=[True/False]
	:param 3 : Path where data should be saved
'''


def distance(N, distType, nbCol, colInd):
	if distType == "pearson":
		dist = fct.partial(pearsonDistFromZScores, N, normalizedSharedMatrix[:, colInd])
	else:
		dist = fct.partial(sciDist.euclidean, normalizedSharedMatrix[:, colInd])

	for i in range(colInd, nbCol):
		# Add value to array of distance, given that dist(a, b) == dist(b, a)
		distSharedMatrix[colInd][i] = distSharedMatrix[i][colInd] = dist(normalizedSharedMatrix[:, i])


def distMatrix(distType):
	'''
	:param mat: implicit, global M x N matrix from which to compute distance between columns
	:param distType: euclidean or pearson
	:return: modifies a global shared matrix of distance
	'''
	shape = normalizedSharedMatrix.shape
	with mp.Pool(23) as P:
		P.map(fct.partial(distance, shape[0], distType, shape[1]), range(shape[1]))


def normalize(index):
	normalizedSharedMatrix[:, index] = sciStt.zscore(mat.iloc[:, index].values, ddof=1)


def normalizeByCol():
	with mp.Pool(23) as P:
		P.map(normalize, range(mat.shape[1]))


def pearsonDistFromZScores(N, u, v):
	# 1 - pearson's r, where r is expressed as the sum of the product of each corresponding z score,
	# divided by the size-1 (correction for sample)
	return 1.0 - np.divide(np.sum(u*v), N-1)


def readData(path):
	return pd.read_csv(path, header=0, index_col=0, float_precision='high')


def transfer(index):
	normalizedSharedMatrix[:, index] = mat.iloc[:, index].values


def transferByCol():
	with mp.Pool(23) as P:
		P.map(transfer, range(mat.shape[1]))


def main():

	print(mat.head())

	if(toNormalize):
		#### Normalize ####
		print("Normalizing data... " + str(datetime.now()))
		# Acts on global variable
		normalizeByCol()
	else:
		#### Transfer data to shared object ####
		print("Preparing data... " + str(datetime.now()))
		# Acts on global variable
		transferByCol()

	#### Compute distance ####
	# Pearson inherently normalizes data, so if we do not want a normalized distance, Pearson should not be used
	if(toNormalize):
		print("Computing pearson distance..." + str(datetime.now()))
		distMatrix("pearson")
		rDist = pd.DataFrame(distSharedMatrix, columns=mat.columns, index=mat.columns)
		print(rDist.head())
		rDist.to_pickle("%srDist_%s" % (savePath, dataPath.split('/')[-1].split('.')[0]))


	print("Computing euclidean distance..." + str(datetime.now()))
	distMatrix("euclidean")
	euclDist = pd.DataFrame(distSharedMatrix, columns=mat.columns, index=mat.columns)
	print(euclDist.head())
	normalizationFlag = '' if toNormalize else 'notNormed_'
	euclDist.to_pickle("%s%seuclideanDist_%s"%(savePath, normalizationFlag, dataPath.split('/')[-1].split('.')[0]))

if __name__ == '__main__':

	#### Read args ####
	dataPath = sys.argv[1]
	toNormalize = sys.argv[2].split('=')[1] == 'True'
	savePath = sys.argv[3]

	#### Importing data to be read by child processes and
	# creating global shared array to be writen into by child processes ####
	print("Importing data... " + str(datetime.now()))
	mat = pd.read_pickle(dataPath)
	normalizedSharedMatrix = np.frombuffer(mp.Array('d', mat.shape[0] * mat.shape[1], lock=False)).reshape((mat.shape[0], mat.shape[1]))
	distSharedMatrix = np.frombuffer(mp.Array('d', mat.shape[1] * mat.shape[1], lock=False)).reshape((mat.shape[1], mat.shape[1]))

	main()
