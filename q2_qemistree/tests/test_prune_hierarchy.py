# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import pandas as pd
import os
import numpy as np
from io import StringIO
from skbio import TreeNode
from q2_qemistree import prune_hierarchy


class TestPruning(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        self.tree = TreeNode.read(["((A, B)C, D)root;"])
        self.no_column = pd.DataFrame(index=['A', 'B', 'D'],
                                      data=[np.nan, np.nan, np.nan])
        self.no_annotation = pd.DataFrame(index=['A', 'B', 'D'],
                                          data=['unclassified', 'unclassified',
                                                'unclassified'],
                                          columns=['class'])
        self.one_annotation = pd.DataFrame(index=['A', 'B', 'D'],
                                           data=['unclassified', 'unclassified',
                                                 'foo'],
                                           columns=['class'])
        self.no_overlap = pd.DataFrame(index=['X', 'Y', 'Z'],
                                       data=['foo', 'bar', 'unclassified'],
                                       columns=['class'])
        self.one_overlap = pd.DataFrame(index=['A', 'Y', 'Z'],
                                        data=['foo', 'bar', 'unclassified'],
                                        columns=['class'])

    def test_no_smiles(self):
        msg = "Feature data does not contain the column 'smiles'. Molecular "
        "hierarchy could not be pruned."
        with self.assertRaisesRegex(ValueError, msg):
            prune_hierarchy(self.no_column, self.tree, prune_type='smiles')

    def test_no_classyfire(self):
        msg = 'Classyfire level = class not present in feature data. '
        'Molecular hierarchy could not be pruned.'
        with self.assertRaisesRegex(ValueError, msg):
            prune_hierarchy(self.no_column, self.tree, prune_type='classyfire')

    def test_no_annotation(self):
        msg = 'Tree pruning aborted! There are less than two tree tips with '
        'annotations. Please check if correct feature data table was provided.'
        with self.assertRaisesRegex(ValueError, msg):
            prune_hierarchy(self.no_annotation, self.tree,
                            prune_type='classyfire')

    def test_one_annotation(self):
        msg = 'Tree pruning aborted! There are less than two tree tips with '
        'annotations. Please check if correct feature data table was provided.'
        with self.assertRaisesRegex(ValueError, msg):
            prune_hierarchy(self.one_annotation, self.tree,
                            prune_type='classyfire')

    def test_no_overlap(self):
        msg = 'Tree pruning aborted! There are less than two tree tips with '
        'annotations. Please check if correct feature data table was provided.'
        with self.assertRaisesRegex(ValueError, msg):
            prune_hierarchy(self.no_overlap, self.tree,
                            prune_type='classyfire')

    def test_one_overlap(self):
        msg = 'Tree pruning aborted! There are less than two tree tips with '
        'annotations. Please check if correct feature data table was provided.'
        with self.assertRaisesRegex(ValueError, msg):
            prune_hierarchy(self.one_overlap, self.tree,
                            prune_type='classyfire')


if __name__ == '__main__':
    main()
