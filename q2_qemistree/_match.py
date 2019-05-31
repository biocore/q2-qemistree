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


def match_tables(collated_fingerprints: pd.DataFrame,
                 feature_table: biom.Table,
                 feature_data: pd.DataFrame):
    '''
    This function filters the feature table to retain only features with
    fingerprints.

    Parameters
    ----------
    collated_fingerprints : pd.DataFrame
        table containing mass-spec molecular substructures (columns) for each
        mass-spec feature (index)
    feature_table : biom.Table
        feature tables with mass-spec feature intensity per sample
    feature_data : pd.DataFrame
        metadata (row ID, row m/z) about features in the feature table

    Raises
    ------
    ValueError
        If features in collated fingerprint table are not a subset of
        features in ``feature_table``
        If feature identifiers in feature table and feature data do not match
        If 'row m/z' missing in feature data columms
    UserWarning
        If features in collated fingerprint table are not a subset of
        features in corresponding feature table

    Returns
    -------
    pd.DataFrame
        fingerprint table with binary fingerprint vector
    biom.Table
        feature table that is filtered to contain only the features with
        predicted fingerprints.
    pd.DataFrame
        feature data table that is filtered to contain only the features with
        predicted fingerprints.
    '''

    fps = collated_fingerprints.copy()
    allfps = list(fps.index)
    ftable = feature_table.to_dataframe(dense=True)
    allfeatrs = set(ftable.index)
    if 'row m/z' not in feature_data.columns:
        raise ValueError("Feature data does not contain 'row m/z'"")
    if not set(feature_data.index) == allfeatrs:
        raise ValueError('The identifiers in feature data and feature table '
                         ' do not match')
    if not set(allfps).issubset(allfeatrs):
        extra_tips = set(allfps) - set(allfps).intersection(allfeatrs)
        raise ValueError('The following tips were not '
                         'found in the feature table:\n' +
                         ', '.join([str(i) for i in extra_tips]))
    binary_fps = (fps > 0.5).astype(int)
    filtered_fdata = feature_data.reindex(allfps)
    filtered_ftable = ftable.reindex(allfps)

    return binary_fps, filtered_ftable, filtered_fdata
