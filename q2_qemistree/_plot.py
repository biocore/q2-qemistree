#!/usr/bin/env python
# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import warnings
import seaborn as sns
import pandas as pd
import numpy as np
import q2templates
import pkg_resources

from shutil import copyfile
from os.path import join
from q2_types.tree import NewickFormat
from itolapi import Itol

TEMPLATES = pkg_resources.resource_filename('q2_qemistree', 'assets')


def values_to_colors(coloring_category: str, color_palette: str):
    '''This function generates a color map (dict) for unique values in a
    user-specified feature metadata column.'''
    color_map = {}
    annotations = coloring_category.unique()
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


def format_colors(feature_metadata, category, color_palette):
    colors = []

    color_map = values_to_colors(feature_metadata[category], color_palette)

    colors.append('TREE_COLORS')
    colors.append('SEPARATOR TAB')
    colors.append('DATA')

    for idx in feature_metadata.index:
        color = color_map[feature_metadata.loc[idx, category]]

        if feature_metadata.loc[idx, 'annotation_type'] == 'MS2':
            style, width = 'normal', 6
        else:
            style, width = 'dashed', 4

        colors.append('%s\tclade\t%s\t%s\t%s' % (idx, color, style, width))

    return '\n'.join(colors)


def format_labels(feature_metadata, category, ms2_label, parent_mz):
    labels = []

    missing_values = {'unclassified', 'unexpected server response',
                      'SMILE parse error', np.nan}

    labels.append('LABELS')
    labels.append('SEPARATOR TAB')
    labels.append('DATA')

    if ms2_label:
        for idx in feature_metadata.index:
            ms2_compound = feature_metadata.loc[idx, 'ms2_compound']
            if pd.notna(ms2_compound) and not ms2_compound.isspace():
                label = ms2_compound
            else:
                label = feature_metadata.loc[idx, category]

            if parent_mz and label in missing_values:
                label = feature_metadata.loc[idx, parent_mz]

            labels.append('%s\t%s' % (idx, label))
    else:
        for idx in feature_metadata.index:
            label = feature_metadata.loc[idx, category]

            if parent_mz and label in missing_values:
                label = feature_metadata.loc[idx, parent_mz]

            labels.append('%s\t%s' % (idx, label))

    return '\n'.join(labels)


def plot(output_dir: str, table: biom.Table, tree: NewickFormat,
         feature_metadata: pd.DataFrame, category: str,
         color_palette: str = 'Dark2', ms2_label: bool = False,
         parent_mz: str = None) -> None:

    if category not in feature_metadata.columns:
        raise ValueError('Could not find %s in the feature data, the available'
                         ' columns are: %s.' %
                         ', '.join(feature_metadata.columns.astype()))

    color_fp = join(output_dir, 'colors.tsv')
    with open(color_fp, 'w') as fh:
        fh.write(format_colors(feature_metadata, category, color_palette))

    label_fp = join(output_dir, 'labels.tsv')
    with open(label_fp, 'w') as fh:
        fh.write(format_labels(feature_metadata, category, ms2_label,
                               parent_mz))

    # itol won't accept a file unless it has a .tree or .txt extension
    target = join(output_dir, 'qemistree.tree')
    copyfile(str(tree), target)

    # upload the tree, labels and tip coloring files
    itol_uploader = Itol()
    itol_uploader.add_file(target)
    itol_uploader.add_file(label_fp)
    itol_uploader.add_file(color_fp)

    status = itol_uploader.upload()

    if not status:
        raise ValueError(itol_uploader.comm.upload_output)

    url = itol_uploader.get_webpage()
    print(url)

    index_path = join(TEMPLATES, 'index.html')
    q2templates.render(index_path, output_dir, context={'url': url})
