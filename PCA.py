from datetime import datetime
import pandas as pd
import sklearn.decomposition as decomp
import sys


'''
	:param 1: Path to pickled methylation data
	:param 2: Save path
'''


def main():
	# Importing data
	print("Importing data...\t" + str(datetime.now()))
	methData = pd.read_pickle(methPath)


if __name__ == '__main__':
	methPath = sys.argv[1]

	main()