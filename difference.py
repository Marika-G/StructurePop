import pandas as pd
import sys

'''
    :param 1: Pandas pickled combined and named meth dataframe's path
    :param 2: Save path
'''


def main():
    methData = pd.read_pickle(dataPath)
    difference = methData.groupby(level="Sample_Name", axis=1).diff(axis=1)
    difference = difference.xs("Cord_blood", level="Tissue", axis=1)
    difference.to_pickle("%sComBat_Mvals_difference-git_script-named_filtered_18-04-13"%(savePath))


if __name__ == '__main__':
    dataPath = sys.argv[1]
    savePath = sys.argv[2]
    main()