# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import warnings
import hashlib
import numpy as np
from skbio import TreeNode


def match_label(collated_fingerprints: biom.Table, feature_table: biom.Table):
    '''
    This function filters the feature table to retain only features with
    fingerprints. It also relabels with features with md5 hash of its
    binary fingerprint vector.
    '''
    fingerprints = collated_fingerprints.to_dataframe()
    if fingerprints.shape == (0, 0):
        raise ValueError("Cannot have empty fingerprint table")
    table = feature_table.to_dataframe()
    fps = set(fingerprints.index)
    allfeatrs = set(table.index)
    if not fps.issubset(allfeatrs):
        extra_tips = fps-fps.intersection(allfeatrs)
        warnings.warn(UserWarning('The following tips were not '
                                  'found in the feature table:\n' +
                                  ', '.join([str(i) for i in extra_tips])))
    filtered_table = table.reindex(fps)
    fingerprints = (fingerprints > 0.5).astype(int)
    for fid in fingerprints.index:
        md5 = hashlib.sha256(fingerprints[fid].values.tobytes()).hexdigest()
        fingerprints.loc[fid, 'label'] = md5
        filtered_table.loc[fid, 'label'] = md5
    relabeled_fingerprints = fingerprints.groupby(fingerprints.label).first()
    matched_table = filtered_table.groupby(filtered_table.label).sum()
    npfid = np.asarray(matched_table)
    # biom requires that ids be strings
    matched_table = biom.table.Table(data=npfid,
                                     observation_ids=\
                                     matched_table.index.astype(str),
                                     sample_ids=\
                                     matched_table.columns.astype(str))

    return relabeled_fingerprints, matched_table
