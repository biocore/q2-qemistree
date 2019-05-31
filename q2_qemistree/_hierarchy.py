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
import numpy as np
from itertools import combinations_with_replacement
from q2_feature_table import merge

from ._collate_fingerprint import collate_fingerprint
from ._match import match_label
from ._semantics import CSIDirFmt

def jaccard_modified(u,v):
    identical = np.abs(u[-1] == v[-1])
    u = u[:-1]
    v = v[:-1]
    nonzero = np.bitwise_or(u != 0, v != 0)
    unequal_nonzero = np.bitwise_and((u != v), nonzero)
    a = np.double(unequal_nonzero.sum()) + identical
    b = np.double(nonzero.sum()) + 1

    return (a / b) if b != 0 else 0

def build_tree(merged_fps: pd.DataFrame,
               merged_fdata: pd.DataFrame,
               mz_tolerance: Float) -> TreeNode:
    '''
    This function makes a tree of relatedness between mass-spectrometry
    features using molecular substructure fingerprints.
    '''

    distmat = pairwise_distances(X=merged_fps,
                                 Y=None, metric=jaccard_modified)
    distsq = squareform(distmat, checks=False)
    linkage_matrix = linkage(distsq, method='average')
    tree = TreeNode.from_linkage_matrix(linkage_matrix,
                                        merged_fps.index.tolist())
    return tree


def merge_relabel(fps :pd.DataFrame, fts: pd.DataFrame, fdata: pd.DataFrame):
    '''
    This function merges fingerprints, feature table and feature data from
    multiple feature tables.
    '''
    for i, data in enumerate(fdata):
        n = str(i+1)
        fdata[i].index = ['table' + n + '_' + fid for fid in fdata[i].index]
        fps[i].index = ['table' + n + '_' + fid for fid in fps[i].index]
        fts[i].index = ['table' + n + '_' + fid for fid in fts[i].index]
        npfeatures = fts[idx].values
        ft_biom = biom.table.Table(
            data=npfeatures, observation_ids=fts[i].index.astype(str),
            sample_ids=fts[i].columns.astype(str))
        fts[i] = ft_biom
    merged_fps = pd.concat(fps)
    merged_fps.index.name = '#featureID'
    merged_fdata = pd.concat(fdata)
    merged_fdata.index.name = '#featureID'
    merged_ftable = merge(fts, overlap_method='error_on_overlapping_sample')
    merged_ftable.table_id = '#featureID'

    return merged_fps, merged_ftable, merged_fdata


def make_hierarchy(csi_results: CSIDirFmt,
                   feature_tables: biom.Table,
                   feature_data: pd.DataFrame,
                   mz_tolerance: float,
                   qc_properties: bool = False) -> (TreeNode, biom.Table,
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
    feature_data : pd.DataFrame
        metadata (row ID, row m/z) about features in the feature table
    mz_tolerance : float
        maximum allowable tolerance in m/z of parent ions
    qc_properties : bool, default False
        flag to filter molecular properties to keep only PUBCHEM fingerprints

    Raises
    ------
    ValueError
        If number of feature_tables don't match the number of csi_results
        If number of feature_tables don't match the number of feature_data sets
        If any feature table in empty
        If any collated fingerprint table is empty

    Returns
    -------
    skbio.TreeNode
        a tree of relatedness of molecules
    biom.Table
        matched feature table that is filtered to contain only the
        features present in the tree
    pd.DataFrame
        matched feature data
    '''

    fps, fts, fdata = [], [], []
    if len(feature_tables) != len(csi_results):
        raise ValueError("The feature tables and CSI results should have a "
                         "one-to-one correspondance.")
    if len(feature_tables) != len(feature_data):
        raise ValueError("The feature tables and feature data should have a "
                         "one-to-one correspondance.")
    for feature_table, csi_result in zip(feature_tables, csi_results):
        if feature_table.is_empty():
            raise ValueError("Cannot have empty feature table")
        fingerprints = collate_fingerprint(csi_result, qc_properties)
        bin_fps, matched_ftable, matched_fdata = match_tables(fingerprints,
                                                              feature_table,
                                                              feature_data)
        fps.append(bin_fps)
        fts.append(matched_ftable)
        fdata.append(matched_fdata)
    merged_fps, merged_ftable, merged_fdata = merge_relabel(fps, fts, fdata)
    tree = build_tree(merged_fps, merged_fdata, mz_tolerance)

    return tree, merged_fts, merged_fdata
