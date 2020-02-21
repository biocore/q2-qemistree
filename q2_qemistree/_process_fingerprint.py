# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pandas as pd
import numpy as np
import pkg_resources

from ._semantics import CSIDirFmt


data = pkg_resources.resource_filename('q2_qemistree', 'data')


def collate_fingerprint(csi_result: CSIDirFmt, qc_properties: bool = False,
                        metric: str = 'euclidean'):
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
    collated_fps = pd.DataFrame.from_dict(molfp, orient='index')
    if metric == 'jaccard':
        collated_fps = (collated_fps > 0.5).astype(int)
    if collated_fps.shape == (0, 0):
        raise ValueError('Fingerprint file is empty!')
    substructrs = pd.read_csv(os.path.join(csi_result, 'fingerprints.csv'),
                              index_col='relativeIndex', dtype=str, sep='\t')
    collated_fps.index.name = '#featureID'
    collated_fps.columns = substructrs.loc[collated_fps.columns,
                                           'absoluteIndex']
    if qc_properties:
        properties = os.path.join(data, 'molecular_properties.csv')
        properties = pd.read_csv(properties, index_col='absoluteIndex',
                                 sep='\t')
        pubchem_indx = list(properties.loc[properties.type == 'PUBCHEM'].index)
        pubchem_indx = list(map(str, pubchem_indx))
        collated_fps = collated_fps[pubchem_indx]
    return collated_fps


def get_feature_smiles(csi_result: CSIDirFmt, collated_fps: pd.DataFrame,
                       ms2_match: pd.DataFrame = None):
    '''This function gets the SMILES of mass-spec features from
    CSI:FingerID and optionally, MS/MS library match results
    '''
    if isinstance(csi_result, CSIDirFmt):
        csi_result = str(csi_result.get_path())
    csi_summary = os.path.join(csi_result, 'summary_csi_fingerid.csv')
    csi_summary = pd.read_csv(csi_summary, dtype=str,
                              sep='\t').set_index('experimentName')
    smiles = pd.DataFrame(index=collated_fps.index)
    smiles['csi_smiles'] = csi_summary.loc[smiles.index, 'smiles']
    smiles['ms2_smiles'] = np.nan
    smiles['ms2_library_match'] = np.nan
    smiles['ms2_adduct'] = np.nan
    if ms2_match is not None:
        ms2_match.index = ms2_match.index.astype(str)
        ms2_ids = ms2_match.index.intersection(smiles.index)
        smiles['ms2_smiles'] = ms2_match.loc[ms2_ids, 'Smiles']
        smiles['ms2_library_match'] = ms2_match.loc[ms2_ids, 'Compound_Name']
        smiles['ms2_adduct'] = ms2_match.loc[ms2_ids, 'Adduct']
    smiles = smiles.fillna('missing').apply(
        lambda x : x.replace(' ', 'missing'))
    return smiles


def process_csi_results(csi_result: CSIDirFmt, ms2_match: pd.DataFrame = None,
                        qc_properties: bool = False,
                        metric: str = 'euclidean') -> (pd.DataFrame,
                                                       pd.DataFrame):
    '''This function parses CSI:FingerID result to generate tables
    of collated molecular fingerprints and SMILES for mass-spec features
    '''
    collated_fps = collate_fingerprint(csi_result, qc_properties, metric)
    feature_smiles = get_feature_smiles(csi_result, collated_fps, ms2_match)
    return collated_fps, feature_smiles
