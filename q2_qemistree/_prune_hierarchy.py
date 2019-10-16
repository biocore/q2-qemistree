import pandas as pd
import numpy as np
from skbio import TreeNode

def prune_hierarchy(feature_data: pd.DataFrame,
                    feature_fingerprints: pd.DataFrame,
                    prune_type: str = 'classyfire',
                    classyfire_level: str = 'class') -> TreeNode:
    '''This function prunes the tree for visualization. It retains only
    the features that have been annotated for a user-specified `prune_column`
    name.

    Paramters
    ---------
    feature_data : pd.DataFrame
        Feature data table with Classyfire annotations and/or SMILES
    feature_fingerprints : pd.DataFrame
        Table of fingerprints returned by make_hierarchy()
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
        If there are less than 2 annotated features for a given level

    Returns
    -------
    skbio.TreeNode
        a tree of relatedness of molecules with Classyfire annotations
    '''
    if prune_type == 'classyfire':
        if classyfire_level not in feature_data.columns:
            raise ValueError('Classyfire level = %s not present in feature '
                             'data. Molecular hierarchy could not be pruned.'
                             %classyfire_level)
        prune_column = classyfire_level
    if prune_type == 'smiles':
        if 'smiles' not in feature_data.columns:
            raise ValueError("Feature data does not contain the column "
                             "'smiles'. Molecular hierarchy could not be "
                             "pruned.")
        prune_column = 'smiles'
    values = feature_data[[prune_column]]
    values = values.fillna(-1)
    annotated_features = [v
                          for v in values.index
                          if values.loc[v, prune_column]
                          not in [-1, 'unclassified', 'SMILE parse error',
                          'unexpected server response']
                          ]
    feature_overlap = set(annotated_features).intersection(
        feature_fingerprints.index)
    extra_annotations = set(annotated_features) - feature_overlap
    extra_fps = set(feature_fingerprints.index) - feature_overlap
    if len(extra_annotations) > 0:
        warnings.warn('The following annotated features did not have'
                      ' fingerprints:\n' +
                      ', '.join([str(i) for i in extra_annotations])
                      %len(extra_annotations))
    if len(extra_fps) > 0:
        warnings.warn('The following fingerprints did not have'
                      ' annotations:\n' +
                      ', '.join([str(i) for i in extra_fps]))
    if len(feature_overlap) < 2:
        raise ValueError('Less than two annotated features in '
                         '`feature_fingerprints` table. Unable '
                         'to prune hierarchy.')
    classified_fps = feature_fingerprints.loc[feature_overlap]
    return build_tree(classified_fps)
