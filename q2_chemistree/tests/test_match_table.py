# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import os
from biom import load_table
from biom.table import Table
from q2_chemistree import match_table, collatefp, make_hierarchy


class test_match(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        table = Table({}, [], [])
        self.emptyfeatures = table
        table = Table({}, ['a', 'b', 'c'], [])
        self.wrongtips = table
        self.goodtable = os.path.join(THIS_DIR, 'data/features_formated.biom')
        self.goodcsi = os.path.join(THIS_DIR, 'data/goodcsi')
        self.goodthresh = 0.5
        tablefp = collatefp(self.goodcsi)
        treeout = make_hierarchy(tablefp, prob_threshold=self.goodthresh)
        self.goodtree = treeout

    def test_emptyTable(self):
        with self.assertRaises(ValueError):
            match_table(self.goodtree, self.emptyfeatures)

    def test_tipMismatch(self):
        with self.assertRaises(ValueError):
            match_table(self.goodtree, self.wrongtips)

    def test_matchPipeline(self):
        tips = {node.name for node in self.goodtree.tips()}
        features = load_table(self.goodtable)
        tableout = match_table(self.goodtree, features)
        tableids = set(tableout.ids(axis='observation'))
        self.assertEqual(tips, tableids)


if __name__ == '__main__':
    main()
