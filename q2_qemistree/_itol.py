#!/usr/bin/env python
# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
from qiime2 import CategoricalMetadataColumn


def get_itol_barchart(table: biom.Table,
                      metadata: CategoricalMetadataColumn) -> biom.Table:
    '''Generate a table in QIIME 2 artifact format which can be directly
    parsed by iTOL and yield a multi-value bar chart.

    Parameters
    ----------
    table : biom.Table
        Table of feature frequencies by sample.
    metadata : qiime2.CategoricalMetadataColumn
        Categorical sample metadata column.

    Returns
    -------
    biom.Table
        Table of mean feature frequencies per category per sample.
    '''
    column = metadata.drop_missing_values()
    catmap = column.to_series().to_dict()
    return table.collapse(lambda i, _: catmap[i], norm=True, axis='sample')
