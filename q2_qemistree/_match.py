# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import hashlib
import pandas as pd


def match_label(collated_fingerprints: pd.DataFrame,
                feature_table: biom.Table):
    '''
    This function filters the feature table to retain only features with
    fingerprints. It also relabels features with MD5 hash of its
    binary fingerprint vector.

    Parameters
    ----------
    collated_fingerprints : pd.DataFrame
        table containing mass-spec molecular substructures (columns) for each
        mass-spec feature (index)
    feature_table : biom.Table
        feature tables with mass-spec feature intensity per sample.

    Raises
    ------
    ValueError
        If features in collated fingerprint table are not a subset of
        features in ``feature_table``

    Returns
    -------
    pd.DataFrame
        fingerprint table with features relabeled with MD5 hash of
        its binary fingerprint vector
    biom.Table
        feature table that is filtered to contain only the
        features with predicted fingerprints. Features are labeled by MD5 hash
        of its binary fingerprint vector
    pd.DataFrame
        table that maps MD5 hash of a feature to the original feature ID in
        the input feature table
    '''

    fps = collated_fingerprints.copy()
    allfps = list(fps.index)
    if fps.empty:
        raise ValueError("Cannot have empty fingerprint table")
    allfeatrs = set(feature_table.ids('observation'))
    if not set(allfps).issubset(allfeatrs):
        extra_tips = set(allfps) - set(allfps).intersection(allfeatrs)
        raise ValueError('The following tips were not '
                         'found in the feature table:\n' +
                         ', '.join([str(i) for i in extra_tips]))
    filtered_data = [feature_table.data(i, axis='observation') for i in allfps]
    fps = (fps > 0.5).astype(int)
    list_md5 = []
    for fid in allfps:
        md5 = str(hashlib.md5(fps.loc[fid].values.tobytes()).hexdigest())
        list_md5.append(md5)
    fps['label'] = list_md5
    feature_data = pd.DataFrame(columns=['label', '#featureID'])
    feature_data['label'] = list_md5
    feature_data['#featureID'] = allfps
    feature_data.set_index('label', inplace=True)
    relabel_fps = fps.groupby('label').first()
    relabeled_table = biom.table.Table(
        filtered_data, allfps, feature_table.ids('sample'),
        [{'label': str(val)} for val in list_md5])

    def bin_f(id_, x): return x['label']
    collapsed_table = relabeled_table.collapse(
                      bin_f, norm=True, min_group_size=1, axis='observation')
    return relabel_fps, collapsed_table, feature_data
