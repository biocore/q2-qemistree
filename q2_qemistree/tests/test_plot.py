# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import os
import pandas as pd

from q2_qemistree._plot import values_to_colors, format_colors, format_labels


EXP_COLORS = """TREE_COLORS
SEPARATOR TAB
DATA
0	clade	#1b9e77	dashed	4
1	clade	#1b9e77	normal	6
2	clade	#1b9e77	normal	6
3	clade	#1b9e77	normal	6
4	clade	#1b9e77	normal	6
5	clade	#d95f02	dashed	4
6	clade	#d95f02	dashed	4
7	clade	#d95f02	dashed	4
8	clade	#d95f02	dashed	4
9	clade	#d95f02	dashed	4
10	clade	#d95f02	dashed	4"""

EXP_LABELS = """LABELS
SEPARATOR TAB
DATA
0	x
1	x
2	x
3	x
4	x
5	y
6	y
7	y
8	y
9	y
10	y"""


class TestPlot(TestCase):
    def setUp(self):
        data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'data')
        fp = os.path.join(data_dir, 'feature_data_itol.txt')

        self.feature_data = pd.read_csv(fp, sep='\t')

    def test_values_to_colors_repeated(self):
        with self.assertWarns(UserWarning):
            res = values_to_colors(self.feature_data['important'],
                                   color_palette='Dark2')

        # do not check for the exact colors, since these might change
        self.assertEqual(len(res), 10)
        self.assertEqual(set(res.keys()),
                         {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'})
        self.assertEqual(res['a'], res['i'])
        self.assertEqual(res['b'], res['j'])

    def test_values_to_colors(self):
        res = values_to_colors(self.feature_data['not_so_important'],
                               color_palette='Dark2')

        # do not check for the exact colors, since these might change
        self.assertEqual(len(res), 2)
        self.assertEqual(set(res.keys()), {'x', 'y'})

    def test_format_colors(self):
        res = format_colors(self.feature_data, 'not_so_important', 'Dark2')
        self.assertEqual(res, EXP_COLORS)

    def test_format_labels_no_ms2(self):
        res = format_labels(self.feature_data, 'not_so_important', False,
                            None)
        self.assertEqual(res, EXP_LABELS)


if __name__ == '__main__':
    main()
