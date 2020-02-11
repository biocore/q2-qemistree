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


def classyfire_to_colors(coloring_category: str,
                         color_palette: str):
    '''This function generates a color map (dict) for unique Classyfire
    annotations in a user-specified Classyfire level.'''
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


def plot(output_dir: str, table: biom.Table, tree: NewickFormat,
         feature_metadata: pd.DataFrame, category: str,
         color_palette: str = 'Dark2', ms2_label: bool = False,
         parent_mz: str = None) -> None:

    missing_values = {'unclassified', 'unexpected server response',
                      'SMILE parse error', np.nan}

    if category not in feature_metadata.columns:
        raise ValueError('Could not find %s in the feature data')

    color_map = classyfire_to_colors(feature_metadata[category],
                                     color_palette)

    color_fp = join(output_dir, 'colors.tsv')
    with open(color_fp, 'w') as fh:
        fh.write('TREE_COLORS\n'
                 'SEPARATOR TAB\n'
                 'DATA\n')

        for idx in feature_metadata.index:
            color = color_map[feature_metadata.loc[idx, category]]
            if feature_metadata.loc[idx, 'annotation_type'] == 'MS2':
                fh.write(idx + '\t' + 'clade\t' + color + '\tnormal\t6\n')
            else:
                fh.write(idx + '\t' + 'clade\t' + color + '\tdashed\t4\n')

    label_fp = join(output_dir, 'labels.tsv')
    with open(label_fp, 'w') as fh:
        fh.write('LABELS\n'
                 'SEPARATOR TAB\n'
                 'DATA\n')

        if ms2_label:
            for idx in feature_metadata.index:
                ms2_compound = feature_metadata.loc[idx, 'ms2_compound']
                if pd.notna(ms2_compound) and not ms2_compound.isspace():
                    label = ms2_compound
                else:
                    label = feature_metadata.loc[idx, category]
                if parent_mz and label in missing_values:
                    label = feature_metadata.loc[idx, parent_mz]
                fh.write(str(idx) + '\t' + str(label) + '\n')
        else:
            for idx in feature_metadata.index:
                label = feature_metadata.loc[idx, category]
                if parent_mz and label in missing_values:
                    label = feature_metadata.loc[idx, parent_mz]
                fh.write(idx + '\t' + label + '\n')

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
