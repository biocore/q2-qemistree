# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from setuptools import setup, find_packages
import versioneer

setup(
    name='q2-qemistree',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    author='The q2-qemistree development team',
    author_email='a3tripat@ucsd.edu',
    description='Molecular tree inference for metabolomics analysis.',
    license='BSD-2-Clause',
    url="https://github.com/biocore/qemistree",
    entry_points={
        'qiime2.plugins': ['q2-qemistree=q2_qemistree.plugin_setup:plugin']
    },
    package_data={'q2_qemistree': ['data/molecular_properties.csv',
                                   'assets/index.html', 'citations.bib']},
    install_requires=['itolapi']
)
