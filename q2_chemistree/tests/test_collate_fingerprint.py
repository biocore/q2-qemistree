# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from biom import load_table
import pandas as pd
import os
import pkg_resources
from q2_chemistree import collate_fingerprint

data = pkg_resources.resource_filename('q2_chemistree', 'data')

class FingerprintTests(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        self.featureTable = os.path.join(THIS_DIR,
                                         'data/features_formated.biom')
        self.emptycsi = os.path.join(THIS_DIR, 'data/emptycsi')
        self.goodcsi = os.path.join(THIS_DIR, 'data/goodcsi')
        self.properties = pd.read_table(os.path.join(data,
                                        'molecular_properties.csv'), dtype=str)
        self.properties.set_index('absoluteIndex', inplace=True)

    def test_fingerprintOut(self):
        with self.assertRaises(ValueError):
            collate_fingerprint(self.emptycsi)

    def test_featureMatch(self):
        tablefp = collate_fingerprint(self.goodcsi)
        features = load_table(self.featureTable)
        allfeatrs = set(features.ids(axis='observation'))
        fpfeatrs = set(tablefp.ids(axis='observation'))
        self.assertEqual(fpfeatrs <= allfeatrs, True)

    def test_pubchemTrue(self):
        tablefp = collate_fingerprint(self.goodcsi, qc_properties=True)
        indx = self.properties.loc[self.properties.type == 'PUBCHEM'].index
        self.assertEqual(set(tablefp.ids(axis='sample')) == set(indx), True)

    def test_pubchemFalse(self):
        tablefp = collate_fingerprint(self.goodcsi, qc_properties=False)
        indx = self.properties.index
        self.assertEqual(set(tablefp.ids(axis='sample')) == set(indx), True)


if __name__ == '__main__':
    main()
