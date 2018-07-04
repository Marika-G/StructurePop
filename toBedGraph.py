import csv
from datetime import datetime
import pandas as pd
import sys

'''
	:param 1: Path to CpG genomic coordinate for the methylation data
	:param 2: Path to pickled combiend methylation data
	:param 3: Sample name to extract
	:param 4: Header? [Y/N]
	:param 5: Path to save bedGraphs
'''


def main():
    ## Import data
    print("Importing data...\t" + str(datetime.now()))
    coord = pd.read_csv(coordPath, header=0, index_col=0)
    meth = pd.read_pickle(dataPath)
    print("Genomic coordinates:\n%s" % (coord.head()))
    print("Methylation data:\n%s" % (meth.head()))

    ##Duplicate position column (range is 1)
    coord['pos2'] = coord['pos'] + 1

    ##Join methylation data
    print("Joining data...\t" + str(datetime.now()))
    bedInfo = coord.join(meth[sample], how='inner')

    ##Split by tissue
    placenta = bedInfo.drop("Cord_blood", axis=1)
    cord = bedInfo.drop("Placenta", axis=1)

    ##Saving bedGraphs
    print("Saving bedGraphs...\t" + str(datetime.now()))

    if withHeader == "Y":
        # Flags for the browser
        headerP = ["track type=bedGraph name=\"Placenta:%s\" visibility=full" % (sample), "", "", ""]
        placenta.to_csv("%s%s_Placenta_%s.bedgraph" % (savePath, sample, dataPath.split('/')[-1].split('.')[0]), sep='\t',
                       header=headerP, index=False, quoting=csv.QUOTE_NONE)
        headerCB = ["track type=bedGraph name=\"Cord_Blood:%s\" visibility=full" % (sample), "", "", ""]
        cord.to_csv("%s%s_Cord_Blood_%s.bedgraph" % (savePath, sample, dataPath.split('/')[-1].split('.')[0]),
                        sep='\t', header=headerCB, index=False, quoting=csv.QUOTE_NONE)
    else:
        placenta.to_csv("%s%s_Placenta_%s.bedgraph" % (savePath, sample, dataPath.split('/')[-1].split('.')[0]),
                        sep='\t', index=False, quoting=csv.QUOTE_NONE)
        cord.to_csv("%s%s_Cord_Blood_%s.bedgraph" % (savePath, sample, dataPath.split('/')[-1].split('.')[0]),
                    sep='\t', index=False, quoting=csv.QUOTE_NONE)

    print("End\t" + str(datetime.now()))


if __name__ == '__main__':
    coordPath = sys.argv[1]
    dataPath = sys.argv[2]
    sample = int(sys.argv[3])
    withHeader = sys.argv[4]
    savePath = sys.argv[5]

    main()
