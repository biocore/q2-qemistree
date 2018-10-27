# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from biom import load_table
import os
from q2_chemistree import collate_fingerprint


class fingerprintTests(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        self.featureTable = os.path.join(THIS_DIR,
                                         'data/features_formated.biom')
        self.emptycsi = os.path.join(THIS_DIR, 'data/emptycsi')
        self.goodcsi = os.path.join(THIS_DIR, 'data/goodcsi')


    def test_fingerprintOut(self):
        with self.assertRaises(ValueError):
            collate_fingerprint(self.emptycsi)

    def test_featureMatch(self):
        tablefp = collate_fingerprint(self.goodcsi)
        features = load_table(self.featureTable)
        allfeatrs = set(features.ids(axis='observation'))
        fpfeatrs = set(tablefp.ids(axis='observation'))
        self.assertEqual(fpfeatrs <= allfeatrs, True)


if __name__ == '__main__':
    main()
