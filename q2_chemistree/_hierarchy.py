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

from q2_chemistree._collate_fingerprint import collate_fingerprint
from q2_chemistree._match import match_label
from q2_chemistree._semantics import CSIDirFmt


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
                                        list(relabeled_fingerprints.index))
    return tree


def make_hierarchy(csi_result: CSIDirFmt,
                   feature_table: biom.Table,
                   qc_properties: bool = True) -> (TreeNode, biom.Table):
    '''
    This function generates a hierarchy of mass-spec features based on
    predicted chemical fingerprints. It filters the feature to table to
    retain only the features with fingerprints and relables each feature with
    a unique hash of its binary fingerprint vector.

    Parameters
    ----------
    csi_result : CSIDirFmt
        CSI:FingerID output folder
    qc_properties : bool
        flag to filter molecular properties to keep only PUBCHEM fingerprints
    feature_table : biom.Table
        feature table with mass-spec feature intensity per sample

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
        filtered feature table that contains only the features present in
        the tree
    '''

    if feature_table.shape == (0, 0):
        raise ValueError("Cannot have empty feature table")
    fingerprints = collate_fingerprint(csi_result, qc_properties)
    relabeled_fingerprints, matched_feature_table = match_label(fingerprints,
                                                                feature_table)
    tree = build_tree(relabeled_fingerprints)

    return tree, matched_feature_table
