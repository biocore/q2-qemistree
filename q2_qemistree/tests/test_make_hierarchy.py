# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import os
import qiime2
from biom.table import Table
from biom import load_table
from q2_qemistree import make_hierarchy
from q2_qemistree import CSIDirFmt

from q2_qemistree._collate_fingerprint import collate_fingerprint
from q2_qemistree._transformer import _read_dataframe


class TestHierarchy(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        tablefp = Table({}, [], [])
        self.emptyfeatures = tablefp
        goodtable = os.path.join(THIS_DIR, 'data/feature_table1.biom')
        self.features = load_table(goodtable)
        goodfdata = os.path.join(THIS_DIR, 'data/feature_data1.txt')
        self.fdata = _read_dataframe(goodfdata)
        goodtable = os.path.join(THIS_DIR, 'data/feature_table2.biom')
        self.features2 = load_table(goodtable)
        goodfdata = os.path.join(THIS_DIR, 'data/feature_data2.txt')
        self.fdata2 = _read_dataframe(goodfdata)
        self.goodcsi = qiime2.Artifact.load(os.path.join(THIS_DIR,
                                                         'data/goodcsi1.qza'))
        goodcsi = self.goodcsi.view(CSIDirFmt)
        self.collated = collate_fingerprint(goodcsi)
        self.goodcsi2 = qiime2.Artifact.load(os.path.join(
                                            THIS_DIR, 'data/goodcsi2.qza'))
        goodcsi = self.goodcsi2.view(CSIDirFmt)
        self.collated2 = collate_fingerprint(goodcsi)

    def test_unequal_inputs1(self):
        goodcsi = self.goodcsi.view(CSIDirFmt)
        msg = ("The feature tables and CSI results should have a one-to-one"
               " correspondance.")
        with self.assertRaisesRegex(ValueError, msg):
            make_hierarchy([goodcsi], [self.features, self.features2],
                           [self.fdata, self.fdata2])

    def test_unequal_inputs2(self):
        goodcsi = self.goodcsi.view(CSIDirFmt)
        goodcsi2 = self.goodcsi2.view(CSIDirFmt)
        msg = ("The feature tables and feature data should have a one-to-one"
               " correspondance.")
        with self.assertRaisesRegex(ValueError, msg):
            make_hierarchy([goodcsi, goodcsi2],
                           [self.features, self.features2], [self.fdata])

    def test_mergeFeatureDataSingle(self):
        goodcsi1 = self.goodcsi.view(CSIDirFmt)
        treeout, merged_fts, merged_fdata = make_hierarchy(
            [goodcsi1], [self.features], [self.fdata])
        featrs = sorted(list(merged_fts.ids(axis='observation')))
        fdata_featrs = sorted(list(merged_fdata.index))
        self.assertEqual(len(featrs) == 5, True)
        self.assertEqual(fdata_featrs, featrs)

    def test_mergeFeatureDataMultiple(self):
        goodcsi1 = self.goodcsi.view(CSIDirFmt)
        goodcsi2 = self.goodcsi2.view(CSIDirFmt)
        treeout, merged_fts, merged_fdata = make_hierarchy(
            [goodcsi1, goodcsi2], [self.features, self.features2],
            [self.fdata, self.fdata2])
        featrs = sorted(list(merged_fts.ids(axis='observation')))
        fdata_featrs = sorted(list(merged_fdata.index))
        self.assertEqual(len(featrs) == 11, True)
        self.assertEqual(fdata_featrs, featrs)

    def test_emptyFeatures(self):
        goodcsi = self.goodcsi.view(CSIDirFmt)
        with self.assertRaises(ValueError):
            make_hierarchy([goodcsi], [self.emptyfeatures], [self.fdata])

    def test_tipMatchSingle(self):
        goodcsi = self.goodcsi.view(CSIDirFmt)
        treeout, feature_table, merged_fdata = make_hierarchy(
            [goodcsi], [self.features], [self.fdata])
        tip_names = {node.name for node in treeout.tips()}
        self.assertEqual(tip_names, set(feature_table._observation_ids))

    def test_Pipeline(self):
        goodcsi1 = self.goodcsi.view(CSIDirFmt)
        goodcsi2 = self.goodcsi2.view(CSIDirFmt)
        treeout, merged_fts, merged_fdata = make_hierarchy(
            [goodcsi1, goodcsi2], [self.features, self.features2],
            [self.fdata, self.fdata2])
        tip_names = {node.name for node in treeout.tips()}
        self.assertEqual(tip_names, set(merged_fts._observation_ids))


if __name__ == '__main__':
    main()
