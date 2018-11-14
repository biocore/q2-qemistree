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
from q2_chemistree import fingerprint


class fingerprintTests(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        self.badsirpath = os.path.join(THIS_DIR, 'data/foo/bin')
        self.goodsirpath = os.path.join(THIS_DIR, 'data/'
                                        'sirius-linux64-headless-4.0.1/bin')
        self.badionsfp = os.path.join(THIS_DIR, 'data/foo.mgf')
        self.goodionsfp = os.path.join(THIS_DIR, 'data/sirius.mgf')
        self.featureTable = os.path.join(THIS_DIR,
                                         'data/features_formated.biom')
        self.emptycsi = os.path.join(THIS_DIR, 'data/emptycsi')
        self.goodcsi = os.path.join(THIS_DIR, 'data/goodcsi')

    def test_siriusBinaryPath(self):
        with self.assertRaises(OSError):
            fingerprint(self.badsirpath, self.goodionsfp, ppm_max=15,
                        profile='orbitrap', n_jobs=1)

    def test_pipeline(self):
        obs = os.environ['_JAVA_OPTIONS']
        table = fingerprint(self.goodsirpath, self.goodionsfp, ppm_max=15,
                            profile='orbitrap', n_jobs=1)
        features = load_table(self.featureTable)
        allfeatrs = set(features.ids(axis='observation'))
        fpfeatrs = set(table.ids(axis='observation'))
        self.assertEqual(fpfeatrs <= allfeatrs, True)
        self.assertEqual(obs, os.environ['_JAVA_OPTIONS'])

    def test_java_flags(self):
        obs = os.environ['_JAVA_OPTIONS']
        table = fingerprint(self.goodsirpath, self.goodionsfp, ppm_max=15,
                            profile='orbitrap', n_jobs=1,
                            java_flags='-Xms16G')

        features = load_table(self.featureTable)
        allfeatrs = set(features.ids(axis='observation'))
        fpfeatrs = set(table.ids(axis='observation'))
        self.assertEqual(fpfeatrs <= allfeatrs, True)
        self.assertEqual(obs, os.environ['_JAVA_OPTIONS'])


if __name__ == '__main__':
    main()
