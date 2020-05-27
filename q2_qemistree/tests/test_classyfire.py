# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import pandas as pd
from q2_qemistree import get_classyfire_taxonomy


class TestClassyfire(TestCase):
    def setUp(self):
        self.no_smiles = pd.DataFrame(index=['a', 'b', 'c'], data=[1, 2, 3],
                                      columns=['#featureID'])
        self.smiles = pd.DataFrame(index=['a', 'b', 'c'], data=[
            ['missing', 'missing'],
            [' O=C(O)[C@@H](N)Cc1ccccc1', 'missing'],
            ['missing', 'CC(=NC(=O)CC(=NC(=O)C)OOC(=O)C)O']],
            columns=['csi_smiles', 'ms2_smiles'])
        self.nan_smiles = pd.DataFrame(index=['a', 'b', 'c'],
                                       data=[['missing', 'missing'],
                                             ['missing', 'missing'],
                                             ['missing', 'missing']],
                                       columns=['csi_smiles', 'ms2_smiles'])
        self.mal_smiles = pd.DataFrame(index=['a', 'b'],
                                       data=[['missing', 'foo'],
                                             ['bar', 'missing']],
                                       columns=['csi_smiles', 'ms2_smiles'])
        self.levels = set(['kingdom', 'superclass', 'class', 'subclass',
                           'direct_parent', 'structure_source'])

    def test_no_smiles(self):
        msg = ('Feature data table must contain the columns `csi_smiles` '
               'and `ms2_smiles` to run Classyfire')
        with self.assertRaisesRegex(ValueError, msg):
            get_classyfire_taxonomy(self.no_smiles)

    def test_nan_smiles(self):
        msg = ("The feature data table should have at least one structural "
               "annotation to run Classyfire")
        with self.assertRaisesRegex(ValueError, msg):
            get_classyfire_taxonomy(self.nan_smiles)

    def test_malformed_smiles(self):
        msg = "The following structures"
        with self.assertWarnsRegex(UserWarning, msg):
            get_classyfire_taxonomy(self.mal_smiles)

    def test_classyfire_output(self):
        classified = get_classyfire_taxonomy(self.smiles)
        classified_mols = classified[classified['kingdom'] != 'unclassified']
        self.assertTrue(pd.isna(classified_mols).shape, 0)
        self.assertTrue(classified_mols.loc['b',
                        'kingdom'] == 'Organic compounds')
        self.assertTrue((self.levels.issubset(set(classified.columns))))


if __name__ == '__main__':
    main()
