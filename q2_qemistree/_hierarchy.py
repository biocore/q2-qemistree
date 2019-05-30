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

def pairwise(merged_fps, merged_fdata=None, operation):
    X = merged_fps.copy()
    fnames = list(merged_fps.index)
    n_featrs = len(fnames)
    precomputed = pd.DataFrame(index=fnames, columns=fnames, dtype='float')
    iterator = list(combinations_with_replacement(fnames, 2))
    if operation == 'intersect':
        for i, j in iterator:
            precomputed.loc[i, j] = np.bitwise_and(X[i], X[j]).sum()
    if operation == 'union':
        for i, j in iterator:
            precomputed.loc[i, j] = np.bitwise_or(X[i], X[j]).sum()
    if operation == 'mzdiff':
        for i, j in iterator:
            mzi = float(merged_fdata.loc[i, 'row m/z'])
            mzj = float(merged_fdata.loc[j, 'row m/z'])
            precomputed.loc[i, j] = abs(mzi-mzj)

    #TODO: fill lower triangle!!
    return precomputed

def build_tree(merged_fps: pd.DataFrame,
               merged_fdata: pd.DataFrame,
               mz_tolerance: Float) -> TreeNode:
    '''
    This function makes a tree of relatedness between mass-spectrometry
    features using molecular substructure fingerprints.
    '''

    intersect = pairwise(merged_fps, operation='intersect')
    union = pairwise(merged_fps, operation='union')
    mzdiff = pairwise(merged_fps, merged_fdata, operation='mzdiff')

    # iterate over all pairs:
        # if mzdiff > mz_tolerance:
            # distance = 1-intersection/(union+1)
        #else:
            # distance = 1-intersection+1/(union+1)

    ## OLD code
    # distmat = pairwise_distances(X=relabeled_fingerprints,
    #                              Y=None, metric='jaccard')
    # distsq = squareform(distmat, checks=False)
    # linkage_matrix = linkage(distsq, method='average')
    # tree = TreeNode.from_linkage_matrix(linkage_matrix,
    #                                     relabeled_fingerprints.index.tolist())
    return tree


def merge_feature_data(fps :pd.DataFrame, fts: pd.DataFrame,
                       fdata: pd.DataFrame):
    '''
    This function merges fingerprints, feature table and feature data from
    multiple feature tables.
    '''
    for idx, data in enumerate(fdata):
        n_data = str(idx+1)
        fdata[idx].index = [n_data + '_' + fid for fid in fdata[idx].index]
        fps[idx].index = [n_data + '_' + fid for fid in fps[idx].index]
        fts[idx].index = [n_data + '_' + fid for fid in fts[idx].index]
        npfeatures = fts[idx].values
        ft_biom = biom.table.Table(
            data=npfeatures, observation_ids=fts[idx].index.astype(str),
            sample_ids=fts[idx].columns.astype(str))
        fts[idx] = ft_biom
    merged_fps = pd.concat(fps)
    merged_fdata = pd.concat(fdata)
    merged_ftable = merge(ftable, overlap_method='error_on_overlapping_sample')

    return merged_fps, merged_ftable, merged_fdata


def make_hierarchy(csi_results: CSIDirFmt,
                   feature_tables: biom.Table,
                   feature_data: pd.DataFrame,
                   mz_tolerance: Float,
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
        matched_ftable, matched_fdata = match_tables(fingerprints,
                                                     feature_table,
                                                     feature_data)
        fps.append(fingerprints)
        fts.append(matched_ftable)
        fdata.append(matched_fdata)
    merged_fps, merged_ftable, merged_fdata = merge_relabel(fps, fts, fdata)
    tree = build_tree(merged_fps, merged_fdata, mz_tolerance)

    return tree, merged_fts, merged_fdata
