# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from ._version import get_versions

from ._fingerprint import fingerprint, collatefp
from ._hierarchy import make_hierarchy
from ._semantics import MassSpectrometryFeatures

__all__ = ['fingerprint', 'make_hierarchy', 'MassSpectrometryFeatures',
           'collatefp']

__version__ = get_versions()['version']
del get_versions
