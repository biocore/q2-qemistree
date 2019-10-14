# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import requests
import pandas as pd
import warnings


def get_classyfire_taxonomy(feature_data: pd.DataFrame) -> pd.DataFrame:
    '''This function uses structural predictions of molecules (SMILES)
    to run Classyfire and obtain chemical taxonomy for each mass-spec feature.
    It appends chemical taxonomy of features to feature data table.

    Parameters
    ----------
    feature_data : pd.DataFrame
        table that maps MD5 hash of mass-spec features to their structural
        annotations (smiles)

    Raises
    ------
    ValueError
        If feature data does not contain the column 'smiles'

    Returns
    -------
    pd.DataFrame
        feature data table with Classyfire annotations ('kingdom',
        'superclass', 'class','subclass', 'direct_parent') per mass-spec
        feature
    '''
    classyfire_levels = ['kingdom', 'superclass', 'class', 'subclass',
                         'direct_parent']
    if 'smiles' not in feature_data.columns:
        raise ValueError('Feature data table must contain the column `smiles` '
                         'to run Classyfire')
    classyfire = {}
    unexpected = {}
    for idx in feature_data.index:
        smiles = feature_data.loc[idx, 'smiles']
        if pd.notna(smiles):
            url_smiles = 'https://gnps-structure.ucsd.edu/inchikey?smiles='
            inchikey = requests.get(url_smiles+smiles).text
            url_inchi = 'https://gnps-classyfire.ucsd.edu/entities/'
            response = requests.get(url_inchi+str(inchikey)+'.json')
            if response.status_code == 200:
                response = response.json()
                taxonomy = [response[level]['name']
                            if bool(response) and response[level] is not None
                            else 'unclassified'
                            for level in classyfire_levels]
                classyfire[idx] = taxonomy
            elif response.status_code == 404:
                classyfire[idx] = 'unclassified'
            else:
                classyfire[idx] = 'unexpected server response'
                unexpected[smiles] = response.status_code
        else:
            classyfire[idx] = 'unclassified'
    if bool(unexpected):
        warnings.warn('The following structures (keys) produced unexpected '
                      'server response codes (values):\n'+
                      unexpected + '\n' + 'More information about '
                      'response status codes can be found here:\n' +
                      'https://www.ietf.org/assignments/http-status-codes/'
                      'http-status-codes.txt')
    classyfire = pd.DataFrame(classyfire).T
    classyfire.columns = classyfire_levels
    classified_feature_data = pd.concat([feature_data, classyfire],
                                        sort=False, axis=1)
    return classified_feature_data
