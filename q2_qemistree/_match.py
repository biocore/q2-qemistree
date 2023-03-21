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
import os

from ._semantics import SiriusDirFmt

def get_matched_tables(collated_fingerprints: pd.DataFrame,
                       smiles: pd.DataFrame,
                       feature_table: biom.Table,
                       csi_result: SiriusDirFmt):
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
    if isinstance(csi_result, SiriusDirFmt):
        csi_result = str(csi_result.get_path())
    fps = collated_fingerprints.copy()
    allfps = list(fps.index)
    canopus = pd.read_csv(os.path.join(csi_result, 'canopus_compound_summary.tsv'), dtype=str, sep='\t')

    if fps.empty:
        raise ValueError("Cannot have empty fingerprint table")
    table = feature_table.to_dataframe(dense=True)
    # table.index = table.index.map(str) # Caso seja necessário converter essa merda 💩💩
    allfeatrs = set(table.index)
    overlap = list(set(allfps).intersection(allfeatrs))
    #overlap = [str(i) for i in overla] # Caso seja necessário converter essa merda 💩💩
    if not set(allfps).issubset(allfeatrs):
        extra_tips = set(allfps) - set(overlap)
        warnings.warn('The following fingerprints were not '
                      'found in the feature table; removed from qemistree:\n' +
                      ', '.join([str(i) for i in extra_tips]), UserWarning)
    filtered_table = table[table.index.isin(overlap)]
    filtered_fps = fps[fps.index.isin(overlap)]

    canopus['feature_id'] = [id.split('_')[-1] for id in canopus['id']]    # criar a coluna com o id apenas
    canopus.set_index('feature_id', inplace=True)  # transformar essa coluna em index do DF
    filtered_canopus = canopus[canopus.index.isin(overlap)]  # reindexar com base em overlap

    list_md5 = []
    list_fid = []
    df_md5 = pd.DataFrame(index=overlap)
    for fid in overlap:
        md5 = str(hashlib.md5(fps.loc[fid].values.tobytes()).hexdigest())
        list_fid.append(fid)
        list_md5.append(md5)
    # df_md5['feature_id'] = list_fid
    df_md5['label'] = list_md5

    filtered_fps = pd.merge(filtered_fps, df_md5, left_index=True, right_index=True)
    filtered_table = pd.merge(filtered_table, df_md5, left_index=True, right_index=True)
    filtered_canopus = pd.merge(filtered_canopus, df_md5, left_index=True, right_index=True)

    # list_md5 = []
    # for fid in overlap:
    #     md5 = str(hashlib.md5(fps.loc[fid].values.tobytes()).hexdigest())
    #     list_md5.append(md5)
    # filtered_fps['label'] = list_md5
    # filtered_table['label'] = list_md5
    # filtered_canopus['label'] = list_md5
    feature_data = pd.DataFrame(columns=['label', '#featureID', 'csi_smiles',
                                         'ms2_smiles', 'ms2_library_match',
                                         'parent_mass', 'retention_time',
                                         'NPC#pathway', 'NPC#superclass', 'NPC#class',
                                         'ClassyFire#most specific class', 'ClassyFire#level 5',
                                         'ClassyFire#subclass', 'ClassyFire#class', 'ClassyFire#superclass',
                                         'ClassyFire#all classifications'])
    feature_data['label'] = list_md5
    feature_data['#featureID'] = overlap
    feature_data['csi_smiles'] = list(smiles.loc[overlap, 'csi_smiles'])
    feature_data['ms2_smiles'] = list(smiles.loc[overlap, 'ms2_smiles'])
    feature_data['ms2_library_match'] = list(
        smiles.loc[overlap, 'ms2_library_match'])
    feature_data['parent_mass'] = list(smiles.loc[overlap, 'parent_mass'])
    feature_data['retention_time'] = list(smiles.loc[overlap,
                                                     'retention_time'])
    feature_data['NPC#pathway'] = list(filtered_canopus.loc[overlap, 'NPC#pathway'])
    feature_data['NPC#superclass'] = list(filtered_canopus.loc[overlap, 'NPC#superclass'])
    feature_data['NPC#class'] = list(filtered_canopus.loc[overlap, 'NPC#class'])
    feature_data['ClassyFire#most specific class'] = list(
        filtered_canopus.loc[overlap, 'ClassyFire#most specific class'])
    feature_data['ClassyFire#level 5'] = list(filtered_canopus.loc[overlap, 'ClassyFire#level 5'])
    feature_data['ClassyFire#subclass'] = list(filtered_canopus.loc[overlap, 'ClassyFire#subclass'])
    feature_data['ClassyFire#class'] = list(filtered_canopus.loc[overlap, 'ClassyFire#class'])
    feature_data['ClassyFire#superclass'] = list(filtered_canopus.loc[overlap, 'ClassyFire#superclass'])
    feature_data['ClassyFire#all classifications'] = list(
        filtered_canopus.loc[overlap, 'ClassyFire#all classifications'])

    feature_data.set_index('label', inplace=True)
    relabel_fps = filtered_fps.groupby('label').first()
    matched_table = filtered_table.groupby('label').sum()
    # biom requires that ids be strings
    npfeatures = matched_table.values
    matched_table = biom.table.Table(
        data=npfeatures, observation_ids=matched_table.index.astype(str),
        sample_ids=matched_table.columns.astype(str))

    return relabel_fps, matched_table, feature_data
