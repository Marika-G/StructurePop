from datetime import datetime
import functools as fct
import multiprocessing as mp
import numpy as np
import pandas as pd
import scipy.stats as sciStt
import sys


'''
    :Path to combined tissues
    :Correlation type [Pearson/Spearman]
    :Save path
'''


def pearson(matName, correlationsName, index):
    data = globals()[matName]
    correlations = globals()[correlationsName]

    site = data.iloc[index].unstack()
    correlations[index] = sciStt.pearsonr(site.iloc[:, 0].values, site.iloc[:, 1].values)[0]


def spearman(matName, correlationsName, index):
    data = globals()[matName]
    correlations = globals()[correlationsName]

    site = data.iloc[index].unstack()
    correlations[index] = sciStt.spearmanr(site.iloc[:, 0].values, site.iloc[:, 1].values)[0]


def main():
    print("Computing inter-tissue %s correlation for each CpG... "%(corrType) + str(datetime.now()))
    if corrType == 'Pearson':
        with mp.Pool(processes=23) as p:
            p.map(fct.partial(pearson, 'methData', 'correlations'), range(methData.shape[0]))
    elif corrType == 'Spearman':
        with mp.Pool(processes=23) as p:
            p.map(fct.partial(spearman, 'methData', 'correlations'), range(methData.shape[0]))
    else:
        raise ValueError("Correlation type must be 'Spearman' or 'Pearson'")
    global correlations
    correlations = pd.DataFrame(data=correlations, index=methData.index, columns=['Correlation'])

    print("Saving results... " + str(datetime.now()))
    print(correlations.head())
    correlations.to_pickle("%sinter-tissue_site_%s_%s"%(savePath,corrType,dataPath.split('/')[-1].split('.')[0]))


if __name__ == '__main__':
    dataPath = sys.argv[1]
    corrType = sys.argv[2]
    savePath = sys.argv[3]

    print("Importing data... " + str(datetime.now()))
    methData = pd.read_pickle(dataPath)
    print(methData.head())

    correlations = np.frombuffer(mp.Array('d', methData.shape[0], lock=False))

    main()