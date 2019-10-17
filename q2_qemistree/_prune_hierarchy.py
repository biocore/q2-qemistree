# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from skbio import TreeNode


def prune_hierarchy(feature_data: pd.DataFrame,
                    tree: TreeNode,
                    prune_type: str = 'classyfire',
                    classyfire_level: str = 'class') -> TreeNode:
    '''This function prunes the tree for visualization. It retains only
    the features that have been annotated for a user-specified `prune_type`

    Paramters
    ---------
    feature_data : pd.DataFrame
        Feature data table with Classyfire annotations and/or SMILES. Output
        of make_hierarchy() or get_classyfire_taxonomy()
    tree : skbio.TreeNode
        Tree of relatedness of molecules. Output of make_hierarchy()
    prune_type : str, default 'classyfire'
        Metadata category for tree pruning ('smiles' or 'classyfire')
    classyfire_level : str, default 'class'
        One of the Classyfire levels in ['kingdom', 'superclass', 'class',
        'subclass', 'direct_parent']

    Raises
    ------
    ValueError
        If user-specified `classyfire_level` is not in `feature_data.columns`
        If `smiles` not in `feature_data.columns` when `prune_type` = 'smiles'
        If there are less than 2 annotated features for a given `prune_column`

    Returns
    -------
    skbio.TreeNode
        Pruned tree of molecules with annotated tips
    '''
    if prune_type == 'classyfire':
        if classyfire_level not in feature_data.columns:
            raise ValueError('Classyfire level = %s not present in feature '
                             'data. Molecular hierarchy could not be pruned.'
                             % classyfire_level)
        prune_column = classyfire_level
    if prune_type == 'smiles':
        if 'smiles' not in feature_data.columns:
            raise ValueError("Feature data does not contain the column "
                             "'smiles'. Molecular hierarchy could not "
                             "be pruned.")
        prune_column = 'smiles'
    values = feature_data[[prune_column]]
    values = values.fillna(-1)
    annotated_features = [v
                          for v in values.index
                          if values.loc[v, prune_column]
                          not in [-1, 'unclassified',
                                  'unexpected server response',
                                  'SMILE parse error']
                          ]
    tips = {tip.name for tip in tree.tips()}
    overlap = set(annotated_features).intersection(tips)
    if len(overlap) < 2:
        raise ValueError('Tree pruning aborted! There are less than two '
                         'tree tips with annotations. Please check if '
                         'correct feature data table was provided.')
    pruned_tree = tree.shear(overlap)
    return pruned_tree
