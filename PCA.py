from datetime import datetime
import multiprocessing as mp
import pandas as pd
import numpy as np
import scipy.stats as sciStt
import sklearn.decomposition as decomp
import sys


#TODO : merge with distanceMatrix


'''
	:param 1: Path to pickled methylation data
	:param 2: Boolean, normalize or not? Normalize=[True/False]
	:param 3: Save path
'''


def normalize(index):
	normalizedSharedMatrix[:, index] = sciStt.zscore(mat.iloc[:, index].values, ddof=1)


def normalizeByDimension():
	with mp.Pool(23) as P:
		P.map(normalize, range(mat.shape[1]))


def PCA():
	pca = decomp.PCA(n_components=2)  ##NOTE: the randomized svd_solver will be used because nb components << dimensions
	transformed = pca.fit_transform(np.transpose(normalizedSharedMatrix)) #PCA takes arrays shape(n_samples, n_components)

	with open("%sexplained_variance_%s.txt"%(savePath, methPath.split('/')[-1].split('.')[0]), 'w') as file:
		file.write("Explained variance ratio:\n" + str(pca.explained_variance_ratio_))

	return pd.DataFrame(transformed, index=mat.columns, columns=[1,2])


def transfer(index):
	normalizedSharedMatrix[:, index] = mat.iloc[:, index].values


def transferByCol():
	with mp.Pool(23) as P:
		P.map(transfer, range(mat.shape[1]))


def main():
	if (toNormalize):
		#### Normalize ####
		print("Normalizing data...\t" + str(datetime.now()))
		# Acts on global variable
		normalizeByDimension()
	else:
		#### Transfer data to shared object ####
		print("Preparing data... " + str(datetime.now()))
		# Acts on global variable
		transferByCol()

	#### PCA ####
	print("Performing PCA...\t" + str(datetime.now()))
	transformed = PCA()
	print(transformed.head())
	transformed.to_pickle("%sPCA_%s"%(savePath, methPath.split('/')[-1].split('.')[0]))


if __name__ == '__main__':
	methPath = sys.argv[1]
	toNormalize = sys.argv[2].split('=')[1] == 'True'
	savePath= sys.argv[3]

	#### Importing data to be read by child processes and
	# creating global shared array to be writen into by child processes ####
	print("Importing data... " + str(datetime.now()))
	mat = pd.read_pickle(methPath)
	print(mat.head())

	normalizedSharedMatrix = np.frombuffer(mp.Array('d', mat.shape[0] * mat.shape[1], lock=False)).reshape(
		(mat.shape[0], mat.shape[1]))

	main()