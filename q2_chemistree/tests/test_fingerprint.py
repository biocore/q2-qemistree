import functools
from unittest import TestCase, main
import pandas as pd
import os
from q2_chemistree import fingerprint, collatefp


class fingerprintTests(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        self.badsirpath = os.path.join(THIS_DIR, 'data/foo/bin')
        self.goodsirpath = os.path.join(THIS_DIR, 'data/sirius-osx64-4.0.1/bin')
        self.badionsfp = os.path.join(THIS_DIR, 'data/foo.mgf')
        self.goodionsfp = os.path.join(THIS_DIR, 'data/sirius.mgf')
        self.featureTable = os.path.join(THIS_DIR, 'data/features.csv')
        self.emptycsi = os.path.join(THIS_DIR, 'data/emptycsi')
        self.goodcsi = os.path.join(THIS_DIR, 'data/goodcsi')

    def test_siriusBinaryPath(self):
        with self.assertRaises(OSError):
            fingerprint(self.badsirpath, self.goodionsfp, ppm_max=15,
                        profile='orbitrap', n_jobs=1)

    def test_mgfPath(self):
        with self.assertRaises(OSError):
            fingerprint(self.goodsirpath, self.badionsfp, ppm_max=15,
                        profile='orbitrap', n_jobs=1)

    def test_fingerprintOut(self):
        with self.assertRaises(RuntimeError):
            collatefp(self.emptycsi)

    def test_featureMatch(self):
        tablefp = collatefp(self.goodcsi)
        features = pd.read_table(self.featureTable, sep=',', dtype=str)
        allfeatrs = set(features['row ID'])
        fpfeatrs = set(tablefp.to_dataframe().index)
        self.assertEqual(fpfeatrs <= allfeatrs, True)

    def test_pipeline(self):
        tablefp = fingerprint(self.goodsirpath, self.goodionsfp, ppm_max=15,
                              profile='orbitrap', n_jobs=1)
        features = pd.read_table(self.featureTable, sep=',', dtype=str)
        allfeatrs = set(features['row ID'])
        fpfeatrs = set(tablefp.to_dataframe().index)
        self.assertEqual(fpfeatrs <= allfeatrs, True)


if __name__ == '__main__':
    main()
