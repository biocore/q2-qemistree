# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import hashlib


def match_label(collated_fingerprints: biom.Table, feature_table: biom.Table):
    '''
    This function filters the feature table to retain only features with
    fingerprints. It also relabels features with MD5 hash of its
    binary fingerprint vector.
    '''
    fps = collated_fingerprints.to_dataframe(dense=True)
    allfps = set(fps.index)
    if fps.empty:
        raise ValueError("Cannot have empty fingerprint table")
    table = feature_table.to_dataframe(dense=True)
    allfeatrs = set(table.index)
    if not allfps.issubset(allfeatrs):
        extra_tips = allfps-allfps.intersection(allfeatrs)
        raise ValueError('The following tips were not '
                         'found in the feature table:\n' +
                         ', '.join([str(i) for i in extra_tips]))
    filtered_table = table.reindex(allfps)
    fps = (fps > 0.5).astype(int)
    for fid in fps.index:
        md5 = str(hashlib.sha256(fps.loc[fid].values.tobytes()).hexdigest())
        fps.loc[fid, 'label'] = md5
        filtered_table.loc[fid, 'label'] = md5
    relabel_fps = fps.groupby('label').first()
    matched_table = filtered_table.groupby('label').sum()
    # biom requires that ids be strings
    npfeatures = matched_table.values
    matched_table = biom.table.Table(data=npfeatures,
                                     observation_ids=
                                     matched_table.index.astype(str),
                                     sample_ids=
                                     matched_table.columns.astype(str))

    return relabel_fps, matched_table
