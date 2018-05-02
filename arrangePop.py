from datetime import datetime
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.cluster.hierarchy as sciHi
import scipy.spatial.distance as sciDist
import sklearn.manifold as skman
import sys


'''
	:param 1 : Path to pickled pandas data frame, a N x N distance matrix
	:param 2 : Type of distance used to generate initial data ('pearson' ot 'euclidean')
	:param 3 : Path where output fig should be saved
'''

#Adapted from http://www.nxn.se/valent/extract-cluster-elements-by-color-in-python
def get_cluster_classes(den, label='ivl'):
	cluster_idxs = dict()
	for c, pi in zip(den['color_list'], den['icoord']):
		for leg in pi[1:3]:
			i = (leg - 5.0) / 10.0
			if abs(i - int(i)) < 1e-5:
				if c in cluster_idxs.keys():
					cluster_idxs[c].append(int(i))
				else:
					cluster_idxs[c] = [int(i)]

	colorGroup = dict()
	for c, l in cluster_idxs.items():
		i_l = [int(den[label][i]) for i in l]
		colorGroup[c] = i_l

	return colorGroup


def makeDendro(flatDist):
	clusters = sciHi.linkage(flatDist, metric=distMetric, method='average')
	print("Linkage:")
	print(clusters)

	dendro = sciHi.dendrogram(clusters)
	for i in dendro:
		print(i, dendro[i])
	plt.savefig('%sdendro_%s.svg' % (savePath, distDataPath.split('/')[-1]), format='svg')
	# dendro = sciHi.dendrogram(clusters, truncate_mode='level', p=20)
	# plt.savefig('%sdendro_p20_%s.svg' % (savePath, distDataPath.split('/')[-1]), format='svg')
	plt.close()

	return dendro


def makeHeatMap(methMatrix):
	plt.pcolormesh(methMatrix.values, cmap='hot')
	plt.axis([0, 452, 0, 452])
	plt.gca().invert_yaxis()
	plt.gca().xaxis.tick_top()
	plt.colorbar()
	plt.savefig('%sheatmap_%s.svg'%(savePath,distDataPath.split('/')[-1]), format='svg')
	plt.close()


def reorderData(dendroLeaves, matrix):
	colNames = matrix.columns
	indNames = matrix.index
	# Add numerical level to column index
	matrix.columns = pd.MultiIndex.from_arrays([range(len(colNames)), colNames], names=['Num', 'Id'])
	matrix.index = pd.MultiIndex.from_arrays([range(len(indNames)), indNames], names=['Num', 'Id'])
	# Sort columns by the numerical level, based on dendrogram order
	return matrix.reindex(labels=dendroLeaves, axis='columns', level='Num').reindex(labels=dendroLeaves, axis='index', level='Num')


def scatter2D(coordList, groupDict):
	for color, colNum in groupDict.items():
		x, y = zip(*[coordList[i] for i in colNum])
		plt.scatter(x,y,color=color)
	plt.savefig("%sscatter2D_colored_%s.svg" % (savePath, distDataPath.split('/')[-1]), format='svg')
	plt.close()


def scatter3D(coordList, groupDict):
	fig = plt.figure()
	ax = Axes3D(fig)
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	for color, colNum in groupDict.items():
		x, y, z = zip(*[coordList[i] for i in colNum])
		ax.scatter(x,y,z,color=color)
	fig.savefig("%sscatter3D_colored_%s.svg" % (savePath, distDataPath.split('/')[-1]), format='svg')
	plt.close()


def toFlatDistance(distMat):
	return sciDist.squareform(distMat, force='tovector', checks=False)


def main():
	#### Import data ####
	print("\nLoading distance matrix...\t%s"%(str(datetime.now())))
	dist = pd.read_pickle(distDataPath)
	print(dist.shape)
	print(dist.head())

	#### From square dist matrix to condensed ####
	print("\nGenerating condensed distance matrix...\t%s"%(str(datetime.now())))
	flatDist = toFlatDistance(dist)
	print(flatDist)

	#### Make dendrogram ####
	print("\nGenerating dendrogram...\t%s"%(str(datetime.now())))
	dendro = makeDendro(flatDist)
	#print("Leaves:")
	#print(dendro['leaves'])
	clust = get_cluster_classes(dendro)


	for color in clust.keys():
		print("%s : %s\t"%(color, [dist.columns[int(leave)] for leave in clust[color]]))
		#print("%s : %s\t" % (color, [leave for leave in clust[color]]))

	#### Reorder methylation matrix based on the resulting clustering and generate heatmap ####
	# print("\nReordering distance matrix...\t%s" % (str(datetime.now())))
	# reorderedDist = reorderData(dendro['leaves'], dist)
	# print(reorderedDist.shape)
	# print(reorderedDist.head(10))
	#
	# print("\nGenerating heatmap...\t%s"%(str(datetime.now())))
	# makeHeatMap(reorderedDist)

	print("\nComputing points position in new space...\t%s" % (str(datetime.now())))
	MDS3D = skman.MDS(n_components=3, n_jobs=4, dissimilarity='precomputed')
	fit3D = MDS3D.fit_transform(dist.values)

	MDS2D = skman.MDS(n_components=2, n_jobs=4, dissimilarity='precomputed')
	fit2D = MDS2D.fit_transform(dist.values)

	print("\nGenerating scatter plots...\t%s" % (str(datetime.now())))
	scatter3D(fit3D, clust)
	scatter2D(fit2D, clust)

	print("End\t%s"%(str(datetime.now())))


if __name__ == '__main__':
	distDataPath = sys.argv[1]
	distMetric = 'correlation' if sys.argv[2] != 'euclidean' else sys.argv[2]
	savePath = sys.argv[3]

	main()