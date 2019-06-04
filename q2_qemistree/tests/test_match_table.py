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
import qiime2
from biom import load_table

from q2_qemistree._collate_fingerprint import collate_fingerprint
from q2_qemistree._match import match_tables
from q2_qemistree import CSIDirFmt
from q2_qemistree._transformer import _read_dataframe


class TestMatch(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        table = pd.DataFrame()
        self.emptyfps = table
        table = pd.DataFrame(index=['a', 'b', 'c'], data=['a', 'b', 'c'])
        self.wrongtips = table
        goodtable = os.path.join(THIS_DIR, 'data/feature_table1.biom')
        self.features = load_table(goodtable)
        goodfdata = os.path.join(THIS_DIR, 'data/feature_data1.txt')
        self.fdata = _read_dataframe(goodfdata)
        wrongfdata = os.path.join(THIS_DIR, 'data/feature_data2.txt')
        self.wrongfdata = _read_dataframe(wrongfdata)
        missingfdata = os.path.join(THIS_DIR, 'data/fdata1_missing.txt')
        self.missingfdata = _read_dataframe(missingfdata)
        goodcsi = qiime2.Artifact.load(os.path.join(THIS_DIR,
                                                    'data/goodcsi1.qza'))
        self.tablefp = collate_fingerprint(goodcsi.view(CSIDirFmt))

    def test_missingFdata(self):
        msg = "Feature data does not contain 'row m/z'"
        with self.assertRaisesRegex(ValueError, msg):
            match_tables(self.tablefp, self.features, self.missingfdata)

    def test_matchFtableFdata(self):
        msg = "The identifiers in feature data and feature table do not match"
        with self.assertRaisesRegex(ValueError, msg):
            match_tables(self.tablefp, self.features, self.wrongfdata)

    def test_tipMismatch(self):
        msg = "^The following tips were not found in the feature table:"
        with self.assertRaisesRegex(ValueError, msg):
            match_tables(self.wrongtips, self.features, self.fdata)

    def test_matchFdata(self):
        relabeled_fps, matched_ft, matched_fdata = match_tables(self.tablefp,
                                                                self.features,
                                                                self.fdata)
        fdata_featrs = sorted(list(matched_fdata.index))
        featrs = sorted(list(matched_ft.index))
        self.assertEqual(fdata_featrs, featrs)

    def test_matchFps(self):
        relabeled_fps, matched_ft, matched_fdata = match_tables(self.tablefp,
                                                                self.features,
                                                                self.fdata)
        featrs = sorted(list(matched_ft.index))
        fps = sorted(list(relabeled_fps.index))
        self.assertEqual(fps, featrs)


if __name__ == '__main__':
    main()
