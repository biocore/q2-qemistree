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
import warnings


def get_matched_tables(collated_fingerprints: pd.DataFrame,
                       smiles: pd.DataFrame,
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
    smiles: pd.DataFrame
        table containing smiles for each mass-spec feature (index)
    feature_table : biom.Table
        feature tables with mass-spec feature intensity per sample.

    Raises
    ------
    UserWarning
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
    table = feature_table.to_dataframe(dense=True)
    allfeatrs = set(table.index)
    overlap = list(set(allfps).intersection(allfeatrs))
    if not set(allfps).issubset(allfeatrs):
        extra_tips = set(allfps) - set(overlap)
        warnings.warn('The following fingerprints were not '
                      'found in the feature table; removed from qemistree:\n' +
                      ', '.join([str(i) for i in extra_tips]), UserWarning)
    filtered_table = table.reindex(overlap)
    filtered_fps = fps.reindex(overlap)
    list_md5 = []
    for fid in overlap:
        md5 = str(hashlib.md5(fps.loc[fid].values.tobytes()).hexdigest())
        list_md5.append(md5)
    filtered_fps['label'] = list_md5
    filtered_table['label'] = list_md5
    feature_data = pd.DataFrame(columns=['label', '#featureID','csi_smiles',
                                         'ms2_smiles','ms2_library_match',
                                         'parent_mass', 'retention_time'])
    feature_data['label'] = list_md5
    feature_data['#featureID'] = overlap
    feature_data['csi_smiles'] = list(smiles.loc[overlap, 'csi_smiles'])
    feature_data['ms2_smiles'] = list(smiles.loc[overlap, 'ms2_smiles'])
    feature_data['ms2_library_match'] = list(
        smiles.loc[overlap, 'ms2_library_match'])
    feature_data['parent_mass'] = list(smiles.loc[overlap, 'parent_mass'])
    feature_data['retention_time'] = list(smiles.loc[overlap,
                                                     'retention_time'])
    feature_data.set_index('label', inplace=True)
    relabel_fps = filtered_fps.groupby('label').first()
    matched_table = filtered_table.groupby('label').sum()
    # biom requires that ids be strings
    npfeatures = matched_table.values
    matched_table = biom.table.Table(
        data=npfeatures, observation_ids=matched_table.index.astype(str),
        sample_ids=matched_table.columns.astype(str))

    return relabel_fps, matched_table, feature_data
