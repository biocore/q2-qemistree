# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import os
import pandas as pd
import numpy as np
from q2_qemistree import get_classyfire_taxonomy


class TestClassyfire(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        no_smiles = os.path.join(THIS_DIR, 'data/feature_data_no_smiles.txt')
        self.no_smiles = pd.read_csv(no_smiles, sep='\t')
        self.no_smiles.set_index('label')
        smiles = os.path.join(THIS_DIR, 'data/feature_data_smiles.txt')
        self.smiles = pd.read_csv(smiles, sep='\t')
        self.smiles.set_index('label')
        self.nan_smiles = pd.DataFrame(index=['a', 'b', 'c'],
                                       data=[[np.nan, np.nan],
                                             [np.nan, np.nan],
                                             [np.nan, np.nan]],
                                       columns=['csi_smiles', 'ms2_smiles'])
        self.mal_smiles = pd.DataFrame(index=['a', 'b'],
                                       data=[[np.nan, 'foo'],
                                             ['bar', np.nan]],
                                       columns=['csi_smiles', 'ms2_smiles'])
        self.levels = set(['kingdom', 'superclass', 'class', 'subclass',
                          'direct_parent', 'annotation_type'])

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
        classified_feature_data = get_classyfire_taxonomy(self.smiles)
        self.assertTrue((self.levels.issubset(set(
            classified_feature_data.columns))))


if __name__ == '__main__':
    main()
