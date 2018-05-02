## Imports
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import sys
import seaborn as sns


'''
	:param 1: Path to placental DNAm data (pickled pandas df)
	:param 2: Path to placental meta data (csv)
	:param 3: Path to cord blood DNAm data (pickled pandas df)
	:param 4: Path to cord blood meta data (csv)
	:param 5: m for M-values and b for beta values
	:param 6: Path where to save figures
'''


def bivariatedist(df1Means, df2Means):
	jp = sns.jointplot(x=df1Means.values, y=df2Means.values)
	jp.set_axis_labels("Placenta", "Cord blood")
	plt.savefig("%sdistributionMoyennes_%s.svg"%(savePath,valueType), format="svg")
	plt.close()


def boxplotOfMeans(df1Means, df2Means):
	toPlot = pd.DataFrame([df1Means.values, df2Means.values])
	toPlot = toPlot.T
	toPlot.columns = ['Placenta', 'Cord blood']
	ax = sns.boxplot(data=toPlot)
	ax.set(title='DNA methylation levels are lower placenta than in cord blood', xlabel='Tissue',
		   ylabel='Per individual DNAm mean (M-values)')
	plt.savefig("%sboxplotMoyennes_%s.svg"%(savePath,valueType), format="svg")
	plt.close()


def intersectMethData(df1, meta1, df2, meta2):
	### Create multiindex
	df1.columns = pd.MultiIndex.from_tuples([(id, meta1['Sample_Name'][id]) for id in df1.columns],
											  names=['ID', 'Sample_Name'])
	df2.columns = pd.MultiIndex.from_tuples([(id, meta2['Sample_Name'][id]) for id in df2.columns],
											   names=['ID', 'Sample_Name'])

	### Participants
	intersectPart = [name for name in df1.columns.get_level_values('Sample_Name') if
					 name in df2.columns.get_level_values('Sample_Name')]

	### Sites
	intersectSites = [site for site in df1.index if site in df2.index]

	### Filter
	df1Filtered = df1.loc[intersectSites, (slice(None), intersectPart)]
	df2Filtered = df2.loc[intersectSites, (slice(None), intersectPart)]

	return df1Filtered, df2Filtered


def plotMeansDistribution(df1, df2):
	## Distribution of means by column
	df1Means = df1.mean(axis=0)
	df1Means.sort_index(level='Sample_Name', inplace=True)

	df2Means = df2.mean(axis=0)
	df2Means.sort_index(level='Sample_Name', inplace=True)

	bivariatedist(df1Means, df2Means)
	boxplotOfMeans(df1Means, df2Means)


def main():
	## Data loading
	print("Importing data...\t" + str(datetime.now()))
	pMeth = pd.read_pickle(pMethPath)
	pMeta = pd.read_csv(pMetaPath,header=0, index_col=0, float_precision='high')

	cbMeth = pd.read_pickle(cbMethPath)
	cbMeta = pd.read_csv(cbMetaPath,header=0, index_col=0, float_precision='high')

	## Keep only participants and CpG sites present in both tissues
	print("Creating intersection of data...\t" + str(datetime.now()))
	pMeth, cbMeth = intersectMethData(pMeth, pMeta, cbMeth, cbMeta)
	print("Placenta :\n%s"%(pMeth.head()))
	print("Cord blood :\n%s"%(cbMeth.head()))

	## Generate bi-variate distribution and boxplot of the means
	print("Generating distributions of the means...\t" + str(datetime.now()))
	plotMeansDistribution(pMeth, cbMeth)

	print("End\t" + str(datetime.now()))


if __name__ == '__main__':
	pMethPath = sys.argv[1]
	pMetaPath = sys.argv[2]
	cbMethPath = sys.argv[3]
	cbMetaPath = sys.argv[4]
	valueType = "Mvals" if sys.argv[5] in ["M", "m"] else "Bvals" if sys.argv[5] in ["B", "b"] else sys.argv[5]
	savePath = sys.argv[6]

	main()

