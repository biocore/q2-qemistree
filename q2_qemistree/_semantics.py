# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin.model as model
from qiime2.plugin import SemanticType
from q2_types.feature_data import FeatureData
import os
import warnings


def validate_mgf(iterable):

    # ms1 and ms2 count the number of records for the MS1 and MS2 spectra
    read_id, new_feature_id, ms1, ms2 = None, None, 0, 0

    for line in iterable:
        if line.startswith('FEATURE_ID='):
            # get the feature identifier without the new line
            new_feature_id = line.split('=')[1].strip()

            if read_id is None:
                read_id = new_feature_id
            else:
                if new_feature_id != read_id:
                    if ms2 < 1:
                        warnings.warn('At least one feature (Feature ID = '
                                      '"%s") does not have MS2 '
                                      'information' % read_id, UserWarning)

                    # if the previous record is good, then we'll reset
                    # the current state variables and start over
                    read_id = new_feature_id

                    ms1, ms2 = 0, 0

        elif read_id is not None:
            if line.startswith('MSLEVEL=1'):
                ms1 += 1

                if ms1 > 1:
                    raise ValueError('Feature "%s" has more than one MSLEVEL=1'
                                     ' record' % read_id)
            elif line.startswith('MSLEVEL=2'):
                ms2 += 1
            elif line.startswith('END IONS'):
                # after reading a whole record, we should have an MS1, this is
                # assuming that MS1 records always come before any MS2 records
                if ms1 < 1:
                    raise ValueError('Feature "%s" does not have an '
                                     'MSLEVEL=1 record' % read_id)
    return True


class MGFFile(model.TextFileFormat):
    def sniff(self):
        # checks for unique MS1s and at least one MS2 per MS1
        with open(str(self)) as f:
            return validate_mgf(f)


MGFDirFmt = model.SingleFileDirectoryFormat('MGFFile', 'features.mgf', MGFFile)
MassSpectrometryFeatures = SemanticType('MassSpectrometryFeatures')


class TSVMolecules(model.TextFileFormat):
    def sniff(self):
        return True


TSVMoleculesFormat = model.SingleFileDirectoryFormat('TSVMoleculesFormat',
                                                     'feature_data.tsv',
                                                     TSVMolecules)
Molecules = SemanticType('Molecules', variant_of=FeatureData.field['type'])


class OutputDirs(model.DirectoryFormat):

    def get_folder_name(self):
        return 'base-output'

    def get_path(self):
        """Get the path to the directory where the outputs are saved"""
        return os.path.join(str(self.path), self.get_folder_name())

    def validate(self, level=None):
        return True  # :D


class SiriusDirFmt(OutputDirs):
    def get_folder_name(self):
        return 'sirius-output'


SiriusFolder = SemanticType('SiriusFolder')
