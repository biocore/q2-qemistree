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
import biom

from q2_qemistree._plot import (values_to_colors, format_colors,
                                format_labels, format_barplots)


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

EXP_LABELS_A = """LABELS
SEPARATOR TAB
DATA
0	Spectral Match to Bleepbloop
1	Caffeine
2	Fakeiine
3	x
4	x
5	Spectral Match to Glu-Val from NIST14
6	Spectral Match to Glu-Val from NIST14
7	y
8	y
9	y
10	y"""

EXP_LABELS_B = """LABELS
SEPARATOR TAB
DATA
0	Spectral Match to Bleepbloop
1	Caffeine
2	Fakeiine
3	5
4	6
5	Spectral Match to Glu-Val from NIST14
6	Spectral Match to Glu-Val from NIST14
7	10
8	22
9	23
10	24"""

EXP_LABELS_C = """LABELS
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

EXP_BARS = """DATASET_MULTIBAR
SEPARATOR TAB
DATASET_LABEL	Relative Abundance
FIELD_COLORS	#f77189	#36ada4
FIELD_LABELS	group dogs	group cats
LEGEND_TITLE	Relative Abundance
LEGEND_SHAPES	1	1
LEGEND_COLORS	#f77189	#36ada4
LEGEND_LABELS	group dogs	group cats
WIDTH	100
DATA
57135347ed549fa52376d5cca207d57c	0.25	0.75
a443a84cf2c49a6b758208e113ad1fd3	0.3125	0.6875
6b527eb72120dda52f9b3952d20fc128	0.35	0.65"""

class TestPlot(TestCase):
    def setUp(self):
        data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'data')
        fp = os.path.join(data_dir, 'feature_data_itol.txt')
        self.feature_data = pd.read_csv(fp, sep='\t')
        fp = os.path.join(data_dir, 'grouped_feature_table.biom')
        self.grouped_table = biom.load_table(fp)

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

    def test_format_labels(self):
        print(self.feature_data)
        res = format_labels(self.feature_data, 'not_so_important', True,
                            None)
        self.assertEqual(res, EXP_LABELS_A)

    def test_format_labels_parent_mz(self):
        res = format_labels(self.feature_data, 'ms2_library_match', True,
                            '#featureID')
        self.assertEqual(res, EXP_LABELS_B)

    def test_format_labels_no_ms2(self):
        res = format_labels(self.feature_data, 'not_so_important', False,
                            None)
        self.assertEqual(res, EXP_LABELS_C)

    def test_format_barplots(self):
        res = format_barplots(self.grouped_table)
        self.assertEqual(res, EXP_BARS)


if __name__ == '__main__':
    main()
