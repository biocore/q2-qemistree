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
from skbio import TreeNode
from q2_feature_table import merge

from ._collate_fingerprint import collate_fingerprint
from ._match import match_label
from ._semantics import CSIDirFmt


def build_tree(relabeled_fingerprints: pd.DataFrame) -> TreeNode:
    '''
    This function makes a tree of relatedness between mass-spectrometry
    features using molecular substructure fingerprints.
    '''
    distmat = pairwise_distances(X=relabeled_fingerprints,
                                 Y=None, metric='jaccard')
    distsq = squareform(distmat, checks=False)
    linkage_matrix = linkage(distsq, method='average')
    tree = TreeNode.from_linkage_matrix(linkage_matrix,
                                        relabeled_fingerprints.index.tolist())
    return tree

# def build_network(relabeled_fingerprints: pd.DataFrame) -> pd.DataFrame:
#      '''
#     This function makes a tree of relatedness between mass-spectrometry
#     features using molecular substructure fingerprints.
#     '''
#     distmat = pairwise_distances(X=relabeled_fingerprints,
#                                  Y=None, metric='jaccard')
#     distsq = squareform(distmat, checks=False)
#     linkage_matrix = linkage(distsq, method='average')
#     tree = TreeNode.from_linkage_matrix(linkage_matrix,
#                                         relabeled_fingerprints.index.tolist())
#     return tree


def merge_feature_data(fdata: pd.DataFrame):
    '''
    This function merges feature data from multiple feature tables. The
    resulting table is indexed by MD5 hash mapped to unique feature
    identifiers in the original feature tables.
    '''
    for idx, data in enumerate(fdata):
        data['table_number'] = str(idx+1)
    merged_fdata = pd.concat(fdata)
    dupl_bool = merged_fdata.index.duplicated(keep='first')
    duplicates = merged_fdata[dupl_bool].index.unique()
    if len(duplicates) == 0:
        return merged_fdata
    else:
        for idx in duplicates:
            merged_fdata.loc[idx, '#featureID'] = (',').join(
                list(merged_fdata.loc[idx, '#featureID']))
            merged_fdata.loc[idx, 'table_number'] = (',').join(
                list(merged_fdata.loc[idx, 'table_number']))
        merged_fdata = merged_fdata[~dupl_bool]
        return merged_fdata


def make_hierarchy(csi_results: CSIDirFmt,
                   feature_tables: biom.Table,
                   qc_properties: bool = True) -> (TreeNode, biom.Table,
                                                   pd.DataFrame):
    '''
    This function generates a hierarchy of mass-spec features based on
    predicted chemical fingerprints. It filters the feature table to
    retain only the features with fingerprints and relables each feature with
    a hash (MD5) of its binary fingerprint vector.

    Parameters
    ----------
    csi_results : CSIDirFmt
        one or more CSI:FingerID output folder
    feature_table : biom.Table
        one or more feature tables with mass-spec feature intensity per sample
    qc_properties : bool, default True
        flag to filter molecular properties to keep only PUBCHEM fingerprints

    Raises
    ------
    ValueError
        If ``feature_table`` in empty
        If collated fingerprint table is empty
    UserWarning
        If features in collated fingerprint table are not a subset of
        features in ``feature_table``

    Returns
    -------
    skbio.TreeNode
        a tree of relatedness of molecules
    biom.Table
        merged feature table that is filtered to contain only the
        features present in the tree
    pd.DataFrame
        merged feature data
    '''

    fps, fts, fdata = [], [], []
    if len(feature_tables) != len(csi_results):
        raise ValueError("The feature tables and CSI results should have a "
                         "one-to-one correspondance.")
    for feature_table, csi_result in zip(feature_tables, csi_results):
        if feature_table.is_empty():
            raise ValueError("Cannot have empty feature table")
        fingerprints = collate_fingerprint(csi_result, qc_properties)
        relabeled_fp, matched_ft, feature_data = match_label(fingerprints,
                                                             feature_table)
        fps.append(relabeled_fp)
        fts.append(matched_ft)
        fdata.append(feature_data)
    merged_fdata = merge_feature_data(fdata)
    merged_fps = pd.concat(fps)
    merged_fps = merged_fps[~merged_fps.index.duplicated(keep='first')]
    merged_fts = merge(fts, overlap_method='error_on_overlapping_sample')
    tree = build_tree(merged_fps)

    return tree, merged_fts, merged_fdata

def output_fingerprints(csi_results: CSIDirFmt) -> pd.DataFrame:
    all_fingerprints = collate_fingerprint(csi_results, False)
    all_fingerprints["#FeatureID"] = all_fingerprints.index
    print(all_fingerprints.head())
    return all_fingerprints
    


# def make_network(csi_results: CSIDirFmt,
#                 prob_threshold: float = None,
#                 network_distance_threshold: float = 0.2,
#                 distance_metric: str = 'jaccard') -> pd.DataFrame:
#     '''
#     This function makes a tree of relatedness between mass-spectrometry
#     features using molecular substructure information.

#     Parameters
#     ----------
#     collated_fingerprints : biom.Table
#         biom table containing mass-spec feature IDs (as observations)
#         and molecular substructure IDs (as samples).
#     prob_threshold : float
#         probability value below which a molecular substructure is
#         considered absent from a feature. 'None' for no threshold.
#     network_distance_threshold: float
#         maximum distance for similarity network output
#     distance_metric : str
#         Distance metric to calculate distances between chemical fingerprints
#         for making hierarchy

#     Raises
#     ------
#     ValueError
#         If ``collated_fingerprints`` is empty
#         If ``prob_threshold`` is not in [0,1]

#     Returns
#     -------
#     FingerprintNetworkEdgesFile
#         network edges
#     '''
#     fps = []
#     for feature_table, csi_result in zip(feature_tables, csi_results):
#         if feature_table.is_empty():
#             raise ValueError("Cannot have empty feature table")
#         fingerprints = collate_fingerprint(csi_result, qc_properties)
#         relabeled_fp, matched_ft, feature_data = match_label(fingerprints,
#                                                              feature_table)
#         fps.append(relabeled_fp)
#     merged_fdata = merge_feature_data(fdata)
#     merged_fps = pd.concat(fps)
#     merged_fps = merged_fps[~merged_fps.index.duplicated(keep='first')]

#     if table.shape == (0, 0):
#         raise ValueError("Cannot have empty fingerprint table")
#     if prob_threshold is not None:
#         if not 0 <= prob_threshold <= 1:
#             raise ValueError("Probability threshold is not in [0,1]")
#         table = (table > prob_threshold).astype(int)
#     distmat = pairwise_distances(X=table, Y=None, metric=distance_metric)

#     output_list = []

#     for i in range(distmat.shape[0]):
#         for j in range(distmat.shape[1]):
#             if i == j:
#                 continue
#             if distmat[i][j] < network_distance_threshold:
#                 output_list.append([table.index[i],
#                 table.index[j],
#                 distmat[i][j]])

#     my_pd = pd.DataFrame(output_list, columns=["FeatureID1",
#                                                     "FeatureID2",
#                                                     "Distance"])

#     return my_pd