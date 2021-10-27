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


def collate_fingerprint(csi_result: CSIDirFmt,
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
    substructrs = pd.read_csv(os.path.join(csi_result, 'csi_fingerid.tsv'),
                              index_col='relativeIndex', dtype=str, sep='\t')
    collated_fps.index.name = '#featureID'
    collated_fps.columns = substructrs.loc[collated_fps.columns,
                                           'absoluteIndex']
    return collated_fps


def get_feature_smiles(csi_result: CSIDirFmt, collated_fps: pd.DataFrame,
                       library_match: pd.DataFrame = None):
    '''This function gets the SMILES of mass-spec features from
    CSI:FingerID and optionally, MS/MS library match results
    '''
    if isinstance(csi_result, CSIDirFmt):
        csi_result = str(csi_result.get_path())
    csi_summary = os.path.join(csi_result, 'compound_identifications.tsv')
    csi_summary = pd.read_csv(csi_summary, dtype=str,
                              sep='\t')
    csi_summary.index = csi_summary['id'].str.split('_', 2, expand=True)[2]
    smiles = pd.DataFrame(index=collated_fps.index)
    smiles['csi_smiles'] = csi_summary.reindex(smiles.index)['smiles'].str.strip()
    smiles['ms2_smiles'] = np.nan
    smiles['ms2_library_match'] = np.nan
    smiles['parent_mass'] = np.nan
    smiles['retention_time'] = np.nan
    if library_match is not None:
        library_match.index = library_match.index.astype(str)
        ms2_ids = library_match.index.intersection(smiles.index)
        smiles['ms2_smiles'] = library_match.loc[ms2_ids, 'Smiles'].str.strip()
        smiles['ms2_library_match'] = library_match.loc[ms2_ids, 'LibraryID']
        smiles['parent_mass'] = library_match.loc[ms2_ids, 'parent mass']
        smiles['retention_time'] = library_match.loc[ms2_ids, 'RTConsensus']
    smiles = smiles.fillna('missing').apply(
        lambda x: x.replace({' ': 'missing', '': 'missing'}))
    return smiles


def process_csi_results(csi_result: CSIDirFmt,
                        library_match: pd.DataFrame = None,
                        metric: str = 'euclidean') -> (pd.DataFrame,
                                                       pd.DataFrame):
    '''This function parses CSI:FingerID result to generate tables
    of collated molecular fingerprints and SMILES for mass-spec features
    '''
    collated_fps = collate_fingerprint(csi_result, metric)
    feature_smiles = get_feature_smiles(csi_result, collated_fps,
                                        library_match)
    return collated_fps, feature_smiles
