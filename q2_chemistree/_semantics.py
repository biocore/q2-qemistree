# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin.model as model
from qiime2.plugin import SemanticType
import os


class MGFFile(model.TextFileFormat):
    def sniff(self):
        # we don't really parse this file
        return True


MGFDirFmt = model.SingleFileDirectoryFormat('MGFFile', 'features.mgf', MGFFile)
MassSpectrometryFeatures = SemanticType('MassSpectrometryFeatures')

class OutputDirs(model.DirectoryFormat):

    def get_folder_name(self):
        return 'base-output'

    def get_path(self):
        return os.path.join(str(self.path), self.get_folder_name())

    def validate(self, level=None):
        return True  # :D

class CSIDirFmt(OutputDirs):
    def get_folder_name(self):
        return 'csi-output'


CSIFolder = SemanticType('CSIFolder')


class SiriusDirFmt(OutputDirs):
    def get_folder_name(self):
        return 'sirius-output'


SiriusFolder = SemanticType('SiriusFolder')


class ZodiacDirFmt(OutputDirs):
    def get_folder_name(self):
        return 'zodiac-output'


ZodiacFolder = SemanticType('ZodiacFolder')
