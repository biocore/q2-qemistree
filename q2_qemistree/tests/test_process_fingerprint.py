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
import qiime2

from q2_qemistree import CSIDirFmt
from q2_qemistree._process_fingerprint import (collate_fingerprint,
                                               get_feature_smiles)

data = pkg_resources.resource_filename('q2_qemistree', 'data')


class FingerprintTests(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        self.featureTable = os.path.join(THIS_DIR,
                                         'data/features_formated.biom')
        self.emptycsi = os.path.join(os.path.join(THIS_DIR,
                                                  'data/emptycsi'))
        self.goodcsi = qiime2.Artifact.load(os.path.join(THIS_DIR,
                                                         'data/csiFolder.qza'))
        properties_path = os.path.join(data, 'molecular_properties.tsv')
        self.properties = pd.read_csv(properties_path, dtype=str, sep='\t')
        self.properties.set_index('absoluteIndex', inplace=True)

    def test_fingerprintOut(self):
        msg = "Fingerprint file is empty!"
        with self.assertRaisesRegex(ValueError, msg):
            collate_fingerprint(self.emptycsi)

    def test_featureMatch1(self):
        goodcsi = self.goodcsi.view(CSIDirFmt)
        tablefp = collate_fingerprint(goodcsi)
        features = load_table(self.featureTable)
        allfeatrs = set(features.ids(axis='observation'))
        fpfeatrs = set(tablefp.index)
        self.assertEqual(fpfeatrs <= allfeatrs, True)

    def test_featureMatch2(self):
        goodcsi = self.goodcsi.view(CSIDirFmt)
        tablefp = collate_fingerprint(goodcsi)
        smiles = get_feature_smiles(goodcsi, tablefp)
        fpfeatrs = set(tablefp.index)
        smlfeatrs = set(smiles.index)
        self.assertEqual(fpfeatrs == smlfeatrs, True)

    def test_pubchemFalse(self):
        goodcsi = self.goodcsi.view(CSIDirFmt)
        tablefp = collate_fingerprint(goodcsi)
        indx = self.properties.index
        self.assertEqual(set(tablefp.columns) == set(indx), True)


if __name__ == '__main__':
    main()
