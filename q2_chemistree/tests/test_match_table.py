import functools
from unittest import TestCase, main
from biom.table import Table
import pandas as pd
import os
from q2_chemistree import match_table, collatefp, make_hierarchy


class test_match(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        table = Table({}, [], [])
        self.emptyfeatures = table
        table = pd.DataFrame({'#row ID': ['X', 'Y', 'Z'], 'A': [0, 1, 1],
                             'B': [0, 0, 1], 'C': [0, 0, 0]})
        self.wrongtips = table
        self.wrongformat = os.path.join(THIS_DIR, 'data/features.csv')
        self.goodtable = os.path.join(THIS_DIR, 'data/features_formated.csv')
        self.goodcsi = os.path.join(THIS_DIR, 'data/goodcsi')
        self.goodthresh = 0.5
        tablefp = collatefp(self.goodcsi)
        treeout = make_hierarchy(tablefp, prob_threshold=self.goodthresh)
        self.goodtree = treeout

    def test_emptyTable(self):
        with self.assertRaises(ValueError):
            match_table(self.goodtree, self.emptyfeatures)

    def test_rowID(self):
        features = pd.read_table(self.wrongformat, sep=',', dtype=str)
        with self.assertRaises(RuntimeError):
            match_table(self.goodtree, features)

    def test_tipMismatch(self):
        with self.assertRaises(RuntimeError):
            match_table(self.goodtree, self.wrongtips)

    def test_matchPipeline(self):
        tips = {node.name for node in self.goodtree.tips()}
        features = pd.read_table(self.goodtable, sep=',', dtype=str)
        tableout = match_table(self.goodtree, features)
        tableids =  set(tableout.index)
        self.assertEqual(tips, tableids)


if __name__ == '__main__':
    main()
