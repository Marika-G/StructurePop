import csv
from datetime import datetime
import pandas as pd
import sys


'''
	:param 1: Path to CpG genomic coordinate for the methylation data
	:param 2: Path to pickled methylation data
	:param 3: Sample name to extract
	:param 4: Tissue
	:param 5: Path to save bedGraph
'''


def main():
	## Import data
	print("Importing data...\t" + str(datetime.now()))
	coord = pd.read_csv(coordPath, header=0, index_col=0)
	meth = pd.read_pickle(dataPath)
	print("Genomic coordinates:\n%s"%(coord.head()))
	print("Methylation data:\n%s"%(meth.head()))

	##Duplicate position column (range is 1)
	coord['pos2'] = coord['pos']

	##Join methylation data
	print("Joining data...\t" + str(datetime.now()))
	bedInfo = coord.join(meth[sample], how='inner')

	##Saving bedGraph
	print("Saving bedGraph...\t" + str(datetime.now()))
	# Flags for the browser
	header = ["track type=bedGraph name=\"%s:%s\" visibility=full"%(tissue, sample), "", "", ""]
	bedInfo.to_csv("%s%s_%s.bedgraph"%(savePath, sample, dataPath.split('/')[-1].split('.')[0]), sep='\t', header=header, index=False, quoting=csv.QUOTE_NONE)

	print("End\t" + str(datetime.now()))


if __name__ == '__main__':
	coordPath = sys.argv[1]
	dataPath = sys.argv[2]
	sample = sys.argv[3]
	tissue = sys.argv[4]
	savePath = sys.argv[5]

	main()