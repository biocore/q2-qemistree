# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd
from skbio import TreeNode


def match_table(tree: TreeNode,
                feature_table: pd.DataFrame) -> biom.Table:
    '''
    This function filters the feature table to retain the features present in
        the tree

    Parameters
    ----------
    tree : TreeNode
        skbio TreeNode object representing tree of relatedness
        between molecules
    feature_table : pandas dataframe
            MS1 feature table from MZmine2 with features in rows and
            samples in columns

    Raises
    ------
    ValueError
        If ``feature_table`` has no features
        If ``feature_table`` does not have a column named '#row ID'
        If ``tree`` tips are not a subset of feature names in ``feature_table``

    Returns
    -------
    filtered_feature_table : pandas dataframe
        filtered MS1 feature table that contains only the features present in
        the tree
        Raises

    '''
    if feature_table.empty:
        raise ValueError("Cannot have empty feature table")
    if '#row ID' not in feature_table.columns:
        raise ValueError('Feature table needs to be formatted correctly')
    feature_table = feature_table.set_index('#row ID')
    allfeatrs = set(feature_table.index)
    tip_names = {node.name for node in tree.tips()}
    if not tip_names.issubset(allfeatrs):
        raise ValueError('One or more tree tips not in the feature table')
    common_features = list(allfeatrs.intersection(tip_names))
    filtered_feature_table = feature_table.loc[common_features]

    return filtered_feature_table
