# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import qiime2
import os
from q2_qemistree import MGFDirFmt, SiriusDirFmt, ZodiacDirFmt, OutputDirs
from q2_qemistree import (compute_fragmentation_trees,
                          rerank_molecular_formulas,
                          predict_fingerprints)
from q2_qemistree._fingerprint import artifactory


class FingerprintTests(TestCase):
    def setUp(self):
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        self.badsirpath = os.path.join(THIS_DIR, 'data/foo/bin')
        self.goodsirpath = os.path.join(THIS_DIR, 'sirius.app/Contents/MacOS')
        # MassSpectrometryFeatures
        self.ions = qiime2.Artifact.load(os.path.join(THIS_DIR,
                                                      'data/sirius.mgf.qza'))
        # SiriusFolder
        self.sirout = qiime2.Artifact.load(os.path.join(THIS_DIR,
                                                        'data/sirFolder.qza'))
        # ZodiacFolder
        self.zodout = qiime2.Artifact.load(os.path.join(THIS_DIR,
                                                        'data/zodFolder.qza'))

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
        self.assertTrue(('formula_identifications.tsv' in contents))

        contents = os.listdir(result.path)
        self.assertTrue(('stderr.txt' in contents))
        self.assertTrue(('stdout.txt' in contents))

    def test_fragmentation_trees_negative_ionization(self):
        ions = self.ions.view(MGFDirFmt)
        result = compute_fragmentation_trees(sirius_path=self.goodsirpath,
                                             features=ions,
                                             ppm_max=15, profile='orbitrap',
                                             ions_considered='[M-H]-')
        contents = os.listdir(result.get_path())
        self.assertTrue(('formula_identifications.tsv' in contents))

        contents = os.listdir(result.path)
        self.assertTrue(('stderr.txt' in contents))
        self.assertTrue(('stdout.txt' in contents))

    def test_fragmentation_trees_exception(self):
        ions = self.ions.view(MGFDirFmt)
        pass

    def test_reranking(self):
        ions = self.ions.view(MGFDirFmt)
        sirout = self.sirout.view(SiriusDirFmt)
        result = rerank_molecular_formulas(sirius_path=self.goodsirpath,
                                           fragmentation_trees=sirout,
                                           features=ions)
        contents = os.listdir(result.get_path())
        self.assertTrue(('formula_identifications.tsv' in contents))

        contents = os.listdir(result.path)
        self.assertTrue(('stderr.txt' in contents))
        self.assertTrue(('stdout.txt' in contents))

    def test_fingerid(self):
        zodout = self.zodout.view(ZodiacDirFmt)
        result = predict_fingerprints(sirius_path=self.goodsirpath,
                                      molecular_formulas=zodout, ppm_max=15)
        contents = os.listdir(result.get_path())
        self.assertTrue(('compound_identifications.tsv' in contents))

        contents = os.listdir(result.path)
        self.assertTrue(('stderr.txt' in contents))
        self.assertTrue(('stdout.txt' in contents))


if __name__ == '__main__':
    main()

