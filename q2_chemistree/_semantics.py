# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin.model as model
from qiime2.plugin import SemanticType


class MGFFile(model.TextFileFormat):
    def sniff(self):
        # we don't really parse this file
        return True


MGFDirFmt = model.SingleFileDirectoryFormat('MGFFile', 'features.mgf', MGFFile)
MassSpectrometryFeatures = SemanticType('MassSpectrometryFeatures')


class CSIDirFmt(model.DirectoryFormat):
    def validate(self, level=None):
        return True  # :D


CSIFingerprintFolder = SemanticType('CSIFingerprintFolder')
