# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import requests
import pandas as pd


def get_classyfire_taxonomy(feature_data: pd.DataFrame,
                            classyfire_levels: list = None) -> pd.DataFrame:
    '''This function uses structural predictions of molecules (SMILES)
    to run Classyfire and obtain chemical taxonomy for each mass-spec feature.
    It appends chemical taxonomy of features to feature data table.
    Parameters
    ----------
    feature_data : pd.DataFrame
        table that maps MD5 hash of mass-spec features to their structural
        annotations (smiles)
    feature_table : list
        list of chemical taxonomic levels to obtain from Classyfire.
        Default: ['kingdom', 'superclass', 'class','subclass',
                  'direct_parent']
    Raises
    ------
    ValueError
        If `classyfire_levels` are not a subset of ['kingdom', 'superclass',
        'class', 'subclass', 'direct_parent']
    ValueError
        If feature data does not contain the column 'smiles'
    Returns
    -------
    pd.DataFrame
        feature data table that contains Classyfire annotations per
        mass-spec feature
    '''
    if classyfire_levels is None:
        classyfire_levels = ['kingdom', 'superclass', 'class',
                             'subclass', 'direct_parent']
    if not set(classyfire_levels).issubset({'kingdom', 'superclass', 'class',
                                            'subclass', 'direct_parent'}):
        raise ValueError('One or more Classyfire taxonomic levels are invalid')
    if 'smiles' not in feature_data.columns:
        raise ValueError('Feature data table must contain the column `smiles` '
                         'to run Classyfire')
    classyfire = {}
    for idx in feature_data.index:
        smiles = feature_data.loc[idx, 'smiles']
        if pd.notna(smiles):
            url = 'https://gnps-structure.ucsd.edu/inchikey?smiles=' + smiles
            inchikey = requests.get(url).text
            url = 'https://gnps-classyfire.ucsd.edu/entities/' + str(inchikey) + '.json'
            response = requests.get(url)
            if response.status_code !=  404:
                response = response.json()
                taxonomy = [response[level]['name']
                            if bool(response) and response[level] is not None
                            else 'unclassified'
                            for level in classyfire_levels]
                classyfire[idx] = taxonomy
            else:
                classyfire[idx] = 'unclassified'
        else:
            classyfire[idx] = 'unclassified'
    classyfire = pd.DataFrame.from_dict(classyfire).T
    classyfire.columns = classyfire_levels
    classified_feature_data = pd.concat([feature_data, classyfire],
                                        sort='False', axis =1)
    return classified_feature_data
