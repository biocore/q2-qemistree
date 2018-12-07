# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from biom import load_table
import qiime2
import os
from q2_chemistree import (MGFDirFmt, SiriusDirFmt, ZodiacDirFmt,
                           CSIDirFmt, OutputDirs)
from q2_chemistree import (compute_fragmentation_trees,
                          rerank_molecular_formulas,
                          predict_fingerprints)
from q2_chemistree._fingerprint import artifactory
from qiime2.plugin.model import DirectoryFormat


class FingerprintTests(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        self.badsirpath = os.path.join(THIS_DIR, 'data/foo/bin')
        self.goodsirpath = os.path.join(THIS_DIR, 'data/'
                                        'sirius-osx64-4.0.1/bin')
        # MassSpectrometryFeatures
        self.ions = qiime2.Artifact.load(os.path.join(THIS_DIR,
                                            'data/sirius.mgf.qza'))
        # self.ions = ions.view(MGFDirFmt)
        # SiriusFolder
        self.sirout = qiime2.Artifact.load(os.path.join(THIS_DIR,
                                            'data/SiriusFolder.qza'))
        # self.sirout = sirout.view(SiriusDirFmt)
        #ZodiacFolder
        self.zodout = qiime2.Artifact.load(os.path.join(THIS_DIR,
                                            'data/zodiacFolder.qza'))
        # self.zodout = zodout.view(ZodiacDirFmt)
        #CSIFolder
        csiout = qiime2.Artifact.load(os.path.join(THIS_DIR,
                                            'data/csiFolder.qza'))
        # self.csiout = csiout.view(CSIDirFmt)

    def test_artifactory(self):
        # everything is working fine
        obs = os.environ.get('_JAVA_OPTIONS', '')
        res = artifactory(self.goodsirpath, ['--help'],
                          constructor=OutputDirs, java_flags='-Xms2G')
        self.assertEqual(obs, os.environ.get('_JAVA_OPTIONS'))
        self.assertTrue(isinstance(res, OutputDirs))
        # exceptions are raised
        with self.assertRaises(OSError):
            res = artifactory(self.badsirpath, ['--help'],
                              constructor=OutputDirs)

    def test_fragmentation_trees(self):
        ions = self.ions.view(MGFDirFmt)
        result = compute_fragmentation_trees(sirius_path=self.goodsirpath,
                                             features=ions,
                                             ppm_max=15, profile='orbitrap')
        contents = os.listdir(result.get_path())
        self.assertTrue(('version.txt' in contents))

    def test_reranking(self):
        ions = self.ions.view(MGFDirFmt)
        sirout = self.sirout.view(SiriusDirFmt)
        result = rerank_molecular_formulas(sirius_path=self.goodsirpath,
                                          fragmentation_trees=sirout,
                                          features=ions)
        contents = os.listdir(result.get_path())
        self.assertTrue(('zodiac_summary.csv' in contents))


    def test_fingerid(self):
        zodout = self.zodout.view(ZodiacDirFmt)
        result = predict_fingerprints(sirius_path=self.goodsirpath,
                                     molecular_formulas=zodout,
                                     ppm_max=15)
        contents = os.listdir(result.get_path())
        self.assertTrue(('summary_csi_fingerid.csv' in contents))


if __name__ == '__main__':
    main()
