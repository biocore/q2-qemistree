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
from q2_qemistree import get_classyfire_taxonomy


class TestHierarchy(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        badfdata = os.path.join(THIS_DIR, 'data/feature_data_no_smiles.txt')
        self.badfdata = pd.read_csv(badfdata, sep='\t')
        self.badfdata.set_index('label')

    def test_no_smiles(self):
        msg = 'Feature data table must contain the column `smiles` '
        'to run Classyfire'
        with self.assertRaisesRegex(ValueError, msg):
            get_classyfire_taxonomy(self.badfdata)


if __name__ == '__main__':
    main()
