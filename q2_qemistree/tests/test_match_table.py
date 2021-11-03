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
from biom import load_table
import qiime2

from q2_qemistree import CSIDirFmt
from q2_qemistree._process_fingerprint import collate_fingerprint
from q2_qemistree._match import get_matched_tables


class TestMatch(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        table = pd.DataFrame()
        self.emptyfps = table
        table = pd.DataFrame(index=['a', 'b', 'c'], data=['a', 'b', 'c'])
        self.wrongtips = table
        goodtable = os.path.join(THIS_DIR, 'data/features_formated.biom')
        self.features = load_table(goodtable)
        goodsmiles = os.path.join(THIS_DIR, 'data/features_smiles.txt')
        self.smiles = pd.read_csv(goodsmiles, dtype=str, sep='\t')
        self.smiles = self.smiles.set_index('#featureID')
        self.goodcsi = qiime2.Artifact.load(os.path.join(THIS_DIR, 'data/csiFolder.qza'))
        goodcsi = self.goodcsi.view(CSIDirFmt)
        self.tablefp = collate_fingerprint(goodcsi)

    def test_emptyTable(self):
        msg = "Cannot have empty fingerprint table"
        with self.assertRaisesRegex(ValueError, msg):
            get_matched_tables(self.emptyfps, self.smiles, self.features)

    def test_tipMismatch(self):
        msg = "^The following fingerprints were not found"
        with self.assertWarnsRegex(UserWarning, msg):
            get_matched_tables(self.wrongtips, self.smiles, self.features)

    def test_matchFdata(self):
        relabeled_fps, matched_ft, matched_fdata = get_matched_tables(
            self.tablefp, self.smiles, self.features)
        fdata_featrs = sorted(list(matched_fdata.index))
        fdata_cols = sorted(list(matched_fdata.columns))
        featrs = sorted(list(matched_ft.ids(axis='observation')))
        self.assertEqual(fdata_featrs, featrs)
        self.assertEqual(fdata_cols, sorted(['#featureID', 'csi_smiles',
                                             'ms2_smiles', 'ms2_library_match',
                                             'parent_mass', 'retention_time']))

    def test_matchFps(self):
        relabeled_fps, matched_ft, matched_fdata = get_matched_tables(
            self.tablefp, self.smiles, self.features)
        featrs = sorted(list(matched_ft.ids(axis='observation')))
        fps = sorted(list(relabeled_fps.index))
        self.assertEqual(fps, featrs)


if __name__ == '__main__':
    main()
