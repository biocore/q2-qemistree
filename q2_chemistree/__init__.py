# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from ._version import get_versions

__version__ = get_versions()['version']
del get_versions

from .make_hierarchy import make_hierarchy
from .fingerprint import fingerprint, collatefp, run_command

__all__ = ['make_hierarchy', 'fingerprint', 'collatefp', 'run_command']
