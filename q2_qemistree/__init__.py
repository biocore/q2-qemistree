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
from ._classyfire import get_classyfire_taxonomy
from ._hierarchy import make_hierarchy
from ._prune_hierarchy import prune_hierarchy
from ._semantics import (MassSpectrometryFeatures, MGFDirFmt,
                         CSIFolder, CSIDirFmt, ZodiacFolder, ZodiacDirFmt,
                         SiriusFolder, SiriusDirFmt, OutputDirs)

__all__ = ['compute_fragmentation_trees', 'rerank_molecular_formulas',
           'predict_fingerprints', 'make_hierarchy', 'get_classyfire_taxonomy',
           'prune_hierarchy', 'plot', 'MassSpectrometryFeatures', 'MGFDirFmt',
           'CSIFolder', 'CSIDirFmt', 'ZodiacFolder', 'ZodiacDirFmt',
           'SiriusFolder', 'SiriusDirFmt', 'OutputDirs']

__version__ = get_versions()['version']
