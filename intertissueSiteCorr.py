from datetime import datetime
import functools as fct
import multiprocessing as mp
import numpy as np
import pandas as pd
import scipy.stats as sciStt
import sys


'''
    :Path to combined tissues
    :Save path
'''


def pearson(matName, correlationsName, index):
    data = globals()[matName]
    correlations = globals()[correlationsName]

    site = data.iloc[index].unstack()
    correlations[index] = sciStt.pearsonr(site.iloc[:, 0].values, site.iloc[:, 1].values)[0]


def main():
    print("Computing inter-tissue Pearson correlation for each CpG... " + str(datetime.now()))
    with mp.Pool(processes=2) as p:
        p.map(fct.partial(pearson, 'methData', 'correlations'), range(methData.shape[0]))
    global correlations
    correlations = pd.DataFrame(data=correlations, index=methData.index, columns=['Correlation'])

    print("Saving results... " + str(datetime.now()))
    print(correlations.head())
    correlations.to_pickle("%sinter-tissue_site_correlation_%s"%(savePath,dataPath.split('/')[-1].split('.')[0]))


if __name__ == '__main__':
    dataPath = sys.argv[1]
    savePath = sys.argv[2]

    print("Importing data... " + str(datetime.now()))
    methData = pd.read_pickle(dataPath)
    print(methData.head())

    correlations = np.frombuffer(mp.Array('d', methData.shape[0], lock=False))

    main()