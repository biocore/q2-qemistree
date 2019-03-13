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
from q2_chemistree import make_hierarchy
from q2_chemistree import CSIDirFmt


class TestHierarchy(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        tablefp = Table({}, [], [])
        self.emptyfeatures = tablefp
        goodtable = os.path.join(THIS_DIR, 'data/features_formated.biom')
        self.features = load_table(goodtable)
        self.goodcsi = qiime2.Artifact.load(os.path.join(THIS_DIR,
                                                         'data/csiFolder.qza'))

    def test_emptyFeatures(self):
        goodcsi = self.goodcsi.view(CSIDirFmt)
        with self.assertRaises(ValueError):
            make_hierarchy(goodcsi, self.emptyfeatures)

    def test_tipMatch(self):
        goodcsi = self.goodcsi.view(CSIDirFmt)
        treeout, feature_table = make_hierarchy(goodcsi, self.features)
        tip_names = {node.name for node in treeout.tips()}
        self.assertEqual(tip_names, set(feature_table._observation_ids))


if __name__ == '__main__':
    main()
