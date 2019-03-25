# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage
from skbio import TreeNode
import pandas as pd


def make_hierarchy(collated_fingerprints: biom.Table,
                   prob_threshold: float = None,
                   distance_metric: str = 'jaccard') -> TreeNode:
    '''
    This function makes a tree of relatedness between mass-spectrometry
    features using molecular substructure information.

    Parameters
    ----------
    collated_fingerprints : biom.Table
        biom table containing mass-spec feature IDs (as observations)
        and molecular substructure IDs (as samples).
    prob_threshold : float
        probability value below which a molecular substructure is
        considered absent from a feature. 'None' for no threshold.
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
    skbio.TreeNode
        a tree of relatedness of molecules
    '''
    table = collated_fingerprints.to_dataframe()
    if table.shape == (0, 0):
        raise ValueError("Cannot have empty fingerprint table")
    if prob_threshold is not None:
        if not 0 <= prob_threshold <= 1:
            raise ValueError("Probability threshold is not in [0,1]")

        table = (table > prob_threshold).astype(int)
    distmat = pairwise_distances(X=table, Y=None, metric=distance_metric)
    distsq = squareform(distmat, checks=False)
    linkage_matrix = linkage(distsq, method='average')
    tree = TreeNode.from_linkage_matrix(linkage_matrix, list(table.index))

    return tree


def make_network(collated_fingerprints: biom.Table,
                prob_threshold: float = None,
                network_distance_threshold: float = 0.2,
                distance_metric: str = 'jaccard') -> pd.DataFrame:
    '''
    This function makes a tree of relatedness between mass-spectrometry
    features using molecular substructure information.

    Parameters
    ----------
    collated_fingerprints : biom.Table
        biom table containing mass-spec feature IDs (as observations)
        and molecular substructure IDs (as samples).
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
    table = collated_fingerprints.to_dataframe()
    if table.shape == (0, 0):
        raise ValueError("Cannot have empty fingerprint table")
    if prob_threshold is not None:
        if not 0 <= prob_threshold <= 1:
            raise ValueError("Probability threshold is not in [0,1]")
        table = (table > prob_threshold).astype(int)
    distmat = pairwise_distances(X=table, Y=None, metric=distance_metric)

    output_list = []

    for i in range(distmat.shape[0]):
        for j in range(distmat.shape[1]):
            if i == j:
                continue
            if distmat[i][j] < network_distance_threshold:
                output_list.append([table.index[i],
                table.index[j],
                distmat[i][j]])

    my_pd = pd.DataFrame(output_list, columns=["FeatureID1",
                                                    "FeatureID2",
                                                    "Distance"])

    return my_pd
