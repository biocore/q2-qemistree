# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from ._version import get_versions

from ._fingerprint import (compute_fragmentation_trees,
                           rerank_molecular_formulas,
                           predict_fingerprints)
from ._hierarchy import make_hierarchy
from ._match import match_table
from ._semantics import (MassSpectrometryFeatures, MGFDirFmt,
                         CSIFolder, CSIDirFmt, ZodiacFolder, ZodiacDirFmt,
                         SiriusFolder, SiriusDirFmt, OutputDirs)
from ._collate_fingerprint import collate_fingerprint

__all__ = ['compute_fragmentation_trees', 'rerank_molecular_formulas',
           'predict_fingerprints', 'collate_fingerprint', 'make_hierarchy',
           'match_table', 'MassSpectrometryFeatures', 'MGFDirFmt',
           'CSIFolder', 'CSIDirFmt', 'ZodiacFolder', 'ZodiacDirFmt',
           'SiriusFolder', 'SiriusDirFmt', 'OutputDirs']

__version__ = get_versions()['version']
