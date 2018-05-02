## Imports
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.append("/home/marikagr/.local/lib/python3.6/site-packages/")
import seaborn as sns


def bivariatedist(df1Means, df2Means):
	jp = sns.jointplot(x=df1Means.values, y=df2Means.values)
	jp.set_axis_labels("Placenta", "Cord blood")
	plt.savefig("/home/marikagr/Gen3G/DistributionDonnees/distributionMoyennes_Mvals_filtered_18-04-13.svg", format="svg")
	plt.close()


def boxplotOfMeans(df1Means, df2Means):
	toPlot = pd.DataFrame([df1Means.values, df2Means.values])
	toPlot = toPlot.T
	toPlot.columns = ['Placenta', 'Cord blood']
	ax = sns.boxplot(data=toPlot)
	ax.set(title='DNA methylation levels are lower placenta than in cord blood', xlabel='Tissue',
		   ylabel='Per individual DNAm mean (M-values)')
	plt.savefig("/home/marikagr/Gen3G/DistributionDonnees/boxplotMoyennes_Mvals_filtered_18-04-13.svg", format="svg")
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
	pMeta = pd.read_csv(
		"/mnt/parallel_scratch_mp2_wipe_on_december_2018/jacques/marikagr/Gen3G/Filtered_18-04-13/placenta_meta.csv",
		header=0, index_col=0, float_precision='high')
	pMeth = pd.read_pickle(
		"/mnt/parallel_scratch_mp2_wipe_on_december_2018/jacques/marikagr/Gen3G/Filtered_18-04-13/ComBat_Mvals_placenta_filtered_18-04-13")
	cbMeta = pd.read_csv(
		"/mnt/parallel_scratch_mp2_wipe_on_december_2018/jacques/marikagr/Gen3G/Filtered_18-04-13/cordblood_meta.csv",
		header=0, index_col=0, float_precision='high')
	cbMeth = pd.read_pickle(
		"/mnt/parallel_scratch_mp2_wipe_on_december_2018/jacques/marikagr/Gen3G/Filtered_18-04-13/ComBat_Mvals_cordblood_filtered_18-04-13")

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

	main()

