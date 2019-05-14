# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pandas as pd
import pkg_resources

from ._semantics import CSIDirFmt


data = pkg_resources.resource_filename('q2_qemistree', 'data')


def collate_fingerprint(csi_result: CSIDirFmt, qc_properties: bool = True):
    '''
    This function collates predicted chemical fingerprints for mass-spec
    features in an experiment.
    '''
    if isinstance(csi_result, CSIDirFmt):
        csi_result = str(csi_result.get_path())

    fpfoldrs = os.listdir(csi_result)
    molfp = dict()
    for foldr in fpfoldrs:
        if os.path.isdir(os.path.join(csi_result, foldr)):
            fid = foldr.split('_')[-1]
            fidpath = os.path.join(csi_result, foldr)
            if 'fingerprints' in os.listdir(fidpath):
                fname = os.listdir(os.path.join(fidpath, 'fingerprints'))[0]
                with open(os.path.join(fidpath, 'fingerprints', fname)) as f:
                    fp = f.read().strip().split('\n')
                molfp[fid] = [float(val) for val in fp]
    fingerids = pd.DataFrame.from_dict(molfp, orient='index')
    if fingerids.shape == (0, 0):
        raise ValueError('Fingerprint file is empty!')
    substructrs = pd.read_table(os.path.join(csi_result, 'fingerprints.csv'),
                                index_col='relativeIndex', dtype=str)
    fingerids.index.name = '#featureID'
    fingerids.columns = substructrs.loc[fingerids.columns, 'absoluteIndex']
    if qc_properties is True:
        properties = pd.read_table(os.path.join(data,
                                   'molecular_properties.csv'),
                                   index_col='absoluteIndex')
        pubchem_indx = list(properties.loc[properties.type == 'PUBCHEM'].index)
        pubchem_indx = list(map(str, pubchem_indx))
        fingerids = fingerids[pubchem_indx]
    return fingerids
