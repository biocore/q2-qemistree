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

from q2_chemistree._collate_fingerprint import collate_fingerprint
from q2_chemistree._match import match_label


class TestMatch(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        table = pd.DataFrame()
        self.emptyfps = table
        table = pd.DataFrame(index=['a', 'b', 'c'], data=['a', 'b', 'c'])
        self.wrongtips = table
        goodtable = os.path.join(THIS_DIR, 'data/features_formated.biom')
        self.features = load_table(goodtable)
        goodcsi = os.path.join(THIS_DIR, 'data/goodcsi')
        self.tablefp = collate_fingerprint(goodcsi)

    def test_emptyTable(self):
        msg = "Cannot have empty fingerprint table"
        with self.assertRaisesRegex(ValueError, msg):
            match_label(self.emptyfps, self.features)

    def test_tipMismatch(self):
        msg = "^The following tips were not found in the feature table:"
        with self.assertRaisesRegex(ValueError, msg):
            match_label(self.wrongtips, self.features)

    def test_matchFdata(self):
        relabeled_fps, matched_ft, matched_fdata = match_label(self.tablefp,
                                                               self.features)
        fdata_featrs = sorted(list(matched_fdata.index))
        featrs = sorted(list(matched_ft.ids(axis='observation')))
        self.assertEqual(fdata_featrs, featrs)

    def test_matchFps(self):
        relabeled_fps, matched_ft, matched_fdata = match_label(self.tablefp,
                                                               self.features)
        featrs = sorted(list(matched_ft.ids(axis='observation')))
        fps = sorted(list(relabeled_fps.index))
        self.assertEqual(fps, featrs)


if __name__ == '__main__':
    main()
