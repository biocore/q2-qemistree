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


def values_to_colors(categories, color_palette: str):
    '''This function generates a color map (dict) for unique values in a
    user-specified feature metadata column.'''
    color_map = {}
    colors = sns.color_palette(color_palette,
                               n_colors=len(categories)).as_hex()
    # give a heads up to the user
    if len(set(colors)) < len(categories):
        warnings.warn("The mapping between colors and categories"
                      " is not unique, some colors have been repeated",
                      UserWarning)
    for i, value in enumerate(categories):
        color_map[value] = colors[i]
    return color_map


def format_colors(feature_metadata, category, color_palette):
    colors = []
    annotations = feature_metadata[category].unique()
    color_map = values_to_colors(annotations, color_palette)

    colors.append('TREE_COLORS')
    colors.append('SEPARATOR TAB')
    colors.append('DATA')

    for idx in feature_metadata.index:
        color = color_map[feature_metadata.loc[idx, category]]

        if feature_metadata.loc[idx, 'structure_source'] == 'MS2':
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
            ms2_compound = feature_metadata.loc[idx, 'ms2_library_match']
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


def format_barplots(table: biom.Table):
    barplots = []
    barplots.append('DATASET_MULTIBAR')
    barplots.append('SEPARATOR TAB')
    barplots.append('DATASET_LABEL\tRelative Abundance')

    table = table.norm(axis='observation', inplace=False)
    table = table.to_dataframe(dense=True)

    field_labels = list(table.columns)
    field_colors = values_to_colors(field_labels, 'husl').values()

    barplots.append('FIELD_COLORS\t'+'\t'.join(field_colors))
    barplots.append('FIELD_LABELS\t'+'\t'.join(field_labels))

    barplots.append('LEGEND_TITLE\tRelative Abundance')
    barplots.append('LEGEND_SHAPES\t'+'\t'.join(['1']*len(field_colors)))
    barplots.append('LEGEND_COLORS\t'+'\t'.join(field_colors))
    barplots.append('LEGEND_LABELS\t'+'\t'.join(field_labels))
    barplots.append('WIDTH\t100')

    barplots.append('DATA')
    table = table.reset_index()
    for idx in table.index:
        barplots.append('\t'.join(table.loc[idx].apply(str)))

    return '\n'.join(barplots)


def plot(output_dir: str, tree: NewickFormat, feature_metadata: pd.DataFrame,
         category: str = 'class', ms2_label: bool = True,
         color_palette: str = 'Dark2', parent_mz: str = None,
         grouped_table: biom.Table = None) -> None:
    '''This function plots the phenetic tree in iTOL with clade colors,
    feature labels and relative abundance per sample group.

    Parameters
    ----------
    tree : NewickFormat
        Phenetic tree
    feature_metadata : pd.DataFrame
        Feature metadata
    grouped_table : biom.Table, optional
        Feature table of samples grouped by categories
    category : str, optional
        The feature data column used to color and label the tips.
        Default 'class'
    color_palette : str, optional
        The color palette to use for coloring tips
    ms2_label : bool, optional
        Whether to label the tips with the MS2 value
    parent_mz : str, optional
        If the feature is unclassified, label the tips using this
        column's value

    Raises
    ------
    ValueError
        If `category` is not a column in `feature_metadata`
    UserWarning
        If the number of unique values in `category` is greater than the number
        of unique colors in `color_palette`

    Returns
    -------
    None
    '''

    if category not in feature_metadata.columns:
        raise ValueError('Could not find %s in the feature data, the available'
                         ' columns are: %s.' %
                         (category,
                          ', '.join(feature_metadata.columns.astype())))

    color_fp = join(output_dir, 'colors.txt')
    with open(color_fp, 'w') as fh:
        fh.write(format_colors(feature_metadata, category, color_palette))

    label_fp = join(output_dir, 'labels.txt')
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
    if grouped_table is not None:
        barplot_fp = join(output_dir, 'barplots.txt')
        with open(barplot_fp, 'w') as fh:
            fh.write(format_barplots(grouped_table))
        itol_uploader.add_file(barplot_fp)

    status = itol_uploader.upload()

    if not status:
        raise ValueError(itol_uploader.comm.upload_output)

    url = itol_uploader.get_webpage()
    print(url)

    index_path = join(TEMPLATES, 'index.html')
    q2templates.render(index_path, output_dir,
                       context={'url': url,
                                'has_barplots': grouped_table is not None})
