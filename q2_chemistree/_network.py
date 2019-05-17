# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage

from ._collate_fingerprint import collate_fingerprint
from ._semantics import CSIDirFmt

def make_network(csi_results: CSIDirFmt,
                prob_threshold: float = 0.5,
                network_distance_threshold: float = 0.2,
                distance_metric: str = 'jaccard') -> pd.DataFrame:
    '''
    This function makes a tree of relatedness between mass-spectrometry
    features using molecular substructure information.

    Parameters
    ----------
    csi_results : CSIDirFmt
        one or more CSI:FingerID output folder
    prob_threshold : float
        probability value below which a molecular substructure is
        considered absent from a feature. 'None' for no threshold.
    network_distance_threshold: float
        maximum distance for similarity network output
    distance_metric : str
        Distance metric to calculate distances between chemical fingerprints
        for making hierarchy

    Raises
    ------
    ValueError
        If ``collated_fingerprints`` is empty
        If ``prob_threshold`` is not in [0,1]

    Returns
    -------
    FingerprintNetworkEdgesFile
        network edges
    '''
    
    all_fingerprints = collate_fingerprint(csi_results, False)

    #Rounding to 0 
    if prob_threshold != None:
        for key in all_fingerprints.keys():
            fingerprint_column = all_fingerprints[key]
            all_fingerprints[key] = fingerprint_column.where(fingerprint_column < prob_threshold, 0.0)

    distmat = pairwise_distances(X=all_fingerprints, Y=None, metric=distance_metric)

    output_list = []

    for i in range(distmat.shape[0]):
        for j in range(distmat.shape[1]):
            if i == j:
                continue
            if distmat[i][j] < network_distance_threshold:
                output_list.append([all_fingerprints.index[i],
                all_fingerprints.index[j],
                distmat[i][j]])

    my_pd = pd.DataFrame(output_list, columns=["FeatureID1",
                                                    "FeatureID2",
                                                    "Distance"])

    return my_pd