# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import biom
import pandas as pd
import numpy as np

from ._semantics import CSIDirFmt


def collate_fingerprint(csi_result: CSIDirFmt) -> biom.Table:
    '''
    This function collates chemical fingerprints for mass-spec
    features in an experiment.

    Parameters
    ----------
    csi_result : CSI:FingerID output folder

    Raises
    ------
    ValueError
        If ``fptable`` (collated fingerprint table) is empty

    Returns
    -------
    biom.Table
        biom table containing mass-spec feature IDs (in rows) and molecular
        substructure IDs (in columns). Values are presence (1) or absence (0)
        of a particular substructure.
    '''
    if isinstance(csi_result, CSIDirFmt):
        csi_result = str(csi_result.path)

    fpfoldrs = os.listdir(csi_result)
    molfp = dict()
    for foldr in fpfoldrs:
        if os.path.isdir(os.path.join(csi_result, foldr)):
            fidpath = os.path.join(csi_result, foldr)
            fid = foldr.split('_')[-1]
            if 'fingerprints' in os.listdir(fidpath):
                fname = os.listdir(os.path.join(fidpath, 'fingerprints'))[0]
                with open(os.path.join(fidpath, 'fingerprints', fname)) as f:
                    fp = f.read().strip().split('\n')
                molfp[fid] = fp

    fingerids = pd.DataFrame.from_dict(molfp, orient='index')
    if fingerids.shape == (0, 0):
        raise ValueError('Fingerprint file is empty!')
    fingerids.index.name = '#featureID'
    npfid = np.asarray(fingerids)

    # biom requires that ids be strings
    fptable = biom.table.Table(data=npfid,
                               observation_ids=fingerids.index.astype(str),
                               sample_ids=fingerids.columns.astype(str))
    return fptable
