from datetime import datetime
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sklearn.manifold as skman
import sys


'''
	:param 1: Path to pickled distance matrix
	:param 2: Path to where scatter should be saved
'''


def scatter3D(x,y,z):
	fig = plt.figure()
	ax = Axes3D(fig)
	ax.scatter(x,y,z)
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	#fig.savefig("%sscatter3D_%s.svg" % (savePath, dataName), format='svg')
	plt.close()


def main():
	print("\nLoading distance matrix...\t%s" % (str(datetime.now())))
	distMat = pd.read_pickle(distMatPath)

	print("\nComputing points position in new space...\t%s" % (str(datetime.now())))
	MDS = skman.MDS(n_components=3, n_jobs=4, dissimilarity='precomputed')
	fit = MDS.fit_transform(distMat.values)
	print(fit)

	print("\nGenerating scatter plot...\t%s" % (str(datetime.now())))
	xCoord, yCoord, zCoord = zip(*fit)
	scatter3D(xCoord, yCoord, zCoord)

	print("\nEnd\t%s" % (str(datetime.now())))

if __name__ == '__main__':
	distMatPath = sys.argv[1]
	savePath = sys.argv[2]
	dataName = distMatPath.split('/')[-1]
	main()
