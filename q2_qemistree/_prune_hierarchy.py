# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from skbio import TreeNode


def prune_hierarchy(feature_data: pd.DataFrame, tree: TreeNode,
                    column: str) -> TreeNode:
    '''Prunes the tips of the tree to remove unannotated entries

    Parameters
    ----------
    feature_data : pd.DataFrame
        Feature data table with Classyfire annotations and/or SMILES. Output
        of make_hierarchy() or get_classyfire_taxonomy()
    tree : skbio.TreeNode
        Tree of relatedness of molecules. Output of make_hierarchy()
    column : str
        The column used to prune the phylogeny. All elements with no data in
        this column will be removed from the phylogeny. The choices are one of:
        `'kingdom'`, `'superclass'`, `'class'`, `'subclass'`,
        `'direct_parent'`, `'csi_smiles'` and `'ms2-smiles'`.

    Raises
    ------
    ValueError
        If the column value is not valid, or not contained in the feature data.
        If there are less than 2 annotated features for a given `column`.

    Returns
    -------
    skbio.TreeNode
        Pruned tree of molecules with annotated tips
    '''

    known = ['kingdom', 'superclass', 'class', 'subclass', 'direct_parent',
             'csi_smiles', 'ms2_smiles']

    if column not in feature_data.columns:
        raise ValueError("The feature data does not contain the column '%s'" %
                         column)

    if column not in known:
        raise ValueError('Pruning cannot be applied on column "%s". The only '
                         'options are: %s. If your feature data does not '
                         'include this information, consider using '
                         'get_classyfire_taxonomy.' % (column,
                                                       ', '.join(known)))

    failed_values = {'unclassified', 'unexpected server response',
                     'SMILE parse error'}

    # remove all NA values or missing values
    feature_data = feature_data[~(feature_data[column].isin(failed_values) |
                                  feature_data[column].isna())]

    tips = {tip.name for tip in tree.tips()}
    overlap = feature_data.index.intersection(tips)

    if len(overlap) < 2:
        raise ValueError('Tree pruning aborted! There are less than two tree '
                         'tips with annotations. Please check if the correct '
                         'feature data table was provided.')
    pruned_tree = tree.shear(overlap)
    return pruned_tree
