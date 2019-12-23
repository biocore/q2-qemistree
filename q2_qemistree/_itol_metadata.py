#!/usr/bin/env python
# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import warnings
import pandas as pd
import seaborn as sns
import click
from qiime2 import Artifact


def classyfire_to_colors(classified_feature_data: pd.DataFrame,
                         feature_data_column: str, color_palette: str):
    '''This function generates a color map (dict) for unique Classyfire
    annotations in a user-specified Classyfire level.'''
    color_map = {}
    annotations = classified_feature_data[feature_data_column].unique()
    colors = sns.color_palette(color_palette,
                               n_colors=len(annotations)).as_hex()
    # give a heads up to the user
    if len(set(colors)) < len(annotations):
        warnings.warn("The mapping between colors and annotations"
                      " is not unique, some colors have been repeated",
                      UserWarning)
    for i, value in enumerate(annotations):
        color_map[value] = colors[i]
    return color_map


@click.command()
@click.option('--classified-feature-data', required=True, type=str,
              help='Path to feature data with Classyfire taxonomy.')
@click.option('--feature-data-column', default='class', type=str,
              help="One of the Classyfire levels in ['kingdom', "
              "'superclass', 'class', 'subclass', 'direct_parent']")
@click.option('--color-file-path', default='./itol_colors.txt', type=str,
              help='Path to file with colors specifications for tree clades')
@click.option('--label-file-path', default='./itol_labels.txt', type=str,
              help='Path to file with label specifications for tree tips')
@click.option('--color-palette', default='bright', type=str,
              help='Color palette for tree clades. One of the options'
              ' allowed by seaborn.color_palette()')
def get_itol_visualization(classified_feature_data: str,
                           feature_data_column: str = 'class',
                           ms2_label: bool = True,
                           color_file_path: str = './itol_colors.txt',
                           label_file_path: str = './itol_labels.txt',
                           color_palette: str = 'husl'):
    '''This function creates iTOL metadata files to specify clade colors and
    tip labels based on Classyfire annotations.'''
    fdata = Artifact.load(classified_feature_data).view(pd.DataFrame)
    color_map = classyfire_to_colors(fdata, feature_data_column, color_palette)
    with open(color_file_path, 'w+') as fh:
        fh.write('TREE_COLORS\n'
                 'SEPARATOR TAB\n'
                 'DATA\n')
        for idx in fdata.index:
            color = color_map[fdata.loc[idx, feature_data_column]]
            if fdata.loc[idx, 'annotation_type'] == 'MS2':
                fh.write(idx + '\t' + 'clade\t' +
                         color + '\tnormal\t6\n')
            if fdata.loc[idx, 'annotation_type'] == 'CSIFingerID':
                fh.write(idx + '\t' + 'clade\t' +
                         color + '\tdashed\t4\n')
    with open(label_file_path, 'w+') as fh:
        fh.write('LABELS\n'
                 'SEPARATOR TAB\n'
                 'DATA\n')
        if ms2_label:
            for idx in fdata.index:
                ms2_compound = fdata.loc[idx, 'ms2_compound']
                if pd.notna(ms2_compound) and not ms2_compound.isspace():
                    label = ms2_compound
                else:
                    label = fdata.loc[idx, feature_data_column]
                    fh.write(idx + '\t' + label + '\n')
        else:
            for idx in fdata.index:
                label = fdata.loc[idx, feature_data_column]
                fh.write(idx + '\t' + label + '\n')


if __name__ == '__main__':
    get_itol_visualization()
