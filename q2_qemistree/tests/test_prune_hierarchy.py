# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import pandas as pd
from q2_qemistree import prune_hierarchy


class TestPruning(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        no_smiles = os.path.join(THIS_DIR, 'data/feature_data_no_smiles.txt')
        self.no_smiles = pd.read_csv(no_smiles, sep='\t')
        self.no_smiles.set_index('label')
        # MAKE UP DATA

        # smiles = os.path.join(THIS_DIR, 'data/feature_data_smiles.txt')
        # self.smiles = pd.read_csv(smiles, sep='\t')
        # self.smiles.set_index('label')
        # self.nan_smiles = pd.DataFrame(index=['a', 'b', 'c'],
        #                                data=[np.nan, np.nan, np.nan],
        #                                columns=['smiles'])
        # self.mal_smiles = pd.DataFrame(index=['a', 'b', 'c'],
        #                                data=[np.nan, 'foo', 'bar'],
        #                                columns=['smiles'])
        # self.levels = set(['kingdom', 'superclass', 'class', 'subclass',
        #                   'direct_parent'])

    def test_no_smiles(self):
        # merged fps from make_hierarchy
        # load with feature data no smiles
        pass

    def test_no_classyfire(self):
        # run make_hierarchy
        # don't run classyfire
        pass

    def test_more_fingerprints(self):
        pass

    def test_more_feature_data(self):
        pass

    def test_one_overlap(self):
        pass

    def no_overlap(self)
        pass


if __name__ == '__main__':
    main()
