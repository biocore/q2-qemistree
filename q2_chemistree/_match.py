# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import warnings
from skbio import TreeNode


def match_table(tree: TreeNode,
                feature_table: biom.Table) -> biom.Table:
    '''
    Filters the feature table to retain the features present in the tree.

    Parameters
    ----------
    tree : TreeNode
        skbio TreeNode object representing tree of relatedness
        between molecules
    feature_table : pd.DataFrame
        feature table with features in columns and samples in rows

    Raises
    ------
    ValueError
        If ``feature_table`` has no features
        If ``tree`` tips are not a subset of feature names in ``feature_table``
        If ``filtered_feature_table`` is empty

    Returns
    -------
    biom.Table
        filtered feature table that contains only the features present in
        the tree
    '''
    if feature_table.shape[0] == 0:
        raise ValueError("There are no features in the feature table!")
    allfeatrs = set(feature_table.ids(axis='observation'))
    tip_names = {node.name for node in tree.tips()}
    if not tip_names.issubset(allfeatrs):
        extra_tips = tip_names-tip_names.intersection(allfeatrs)
        warnings.warn(UserWarning('The following tips were not '
                                  'found in the feature table:\n' +
                                  ', '.join([str(i) for i in extra_tips])))
    common_features = list(allfeatrs.intersection(tip_names))
    filtered_feature_table = feature_table.filter(common_features,
                                                  axis='observation',
                                                  inplace=False)
    return filtered_feature_table
