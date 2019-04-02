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


def make_hierarchy(csi_results: CSIDirFmt,
                   feature_tables: biom.Table,
                   qc_properties: bool = True) -> (TreeNode, biom.Table):
    '''
    This function generates a hierarchy of mass-spec features based on
    predicted chemical fingerprints. It filters the feature table to
    retain only the features with fingerprints and relables each feature with
    a hash (MD5) of its binary fingerprint vector.

    Parameters
    ----------
    csi_results : CSIDirFmt
        CSI:FingerID output folder
    feature_table : biom.Table
        feature table with mass-spec feature intensity per sample
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
        filtered feature table that contains only the features present in
        the tree
    '''

    # CHECK the same number of feature tables and fingerprints, error otherwise

    # processed_fingerprints = []
    # processed_feature_tables = []

    # for each pair of feature Table and fingerprints
    #   check if FEATURE TABLE is empty
    #   collate fingerprints
    #   match those fingerprints and feature tables (store them somewhere)

    # MERGE THEM ALL
    #   merge fingerprints as FPS
    #       make sure that this function or code handles one single table
    #   merge feature table (using q2-feature-tabLE)

    # build_tree using the merged fingerprints

    # return tree and feature table
    fps, fts = [], []
    if len(feature_tables) != len(csi_results):
        raise ValueError("The feature tables and CSI results should have a "
                         "one-to-one correspondance.")

    for feature_table, csi_result in zip(feature_tables, csi_results):
        if feature_table.is_empty():
            raise ValueError("Cannot have empty feature table")
        fingerprints = collate_fingerprint(csi_result, qc_properties)
        relabeled_fp, matched_ft = match_label(fingerprints, feature_table)
        fps.append(relabeled_fp)
        fts.append(matched_ft)
    # merge relabeled fps pandas & remove duplicate indices
    merged_fps = pd.concat(relabeled_fps)
    merged_fps = merged_fps[~merged_fps.index.duplicated(keep='first')]
    # merge matched fts qiime merge table
    merged_fts = merge(fts, overlap_method = 'error_on_overlapping_sample')
    tree = build_tree(merged_fps)

    return tree, merged_fts #, pd.DataFrame(['contents to be determined'])
