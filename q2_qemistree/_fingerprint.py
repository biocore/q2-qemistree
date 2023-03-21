# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import subprocess
import os

from ._semantics import MGFDirFmt, SiriusDirFmt
from qiime2.plugin import Str, List


def run_command(cmd, output_fp, error_fp, verbose=True):
    if verbose:
        print("Running external command line application. This may print "
              "messages to stdout and/or stderr.")
        print("The command being run is below. This command cannot "
              "be manually re-run as it will depend on temporary files that "
              "no longer exist.")
        print("\nCommand:", end=' ')
        print(" ".join(cmd), end='\n\n')

    with open(output_fp, 'w') as output_f, open(error_fp, 'w') as error_f:
        subprocess.run(cmd, stdout=output_f, stderr=error_f, check=True)


def artifactory(sirius_path: str, parameters: list, java_flags: str = None,
                constructor=None):
    artifact = constructor()
    if not os.path.exists(sirius_path):
        raise OSError("SIRIUS could not be located")
    sirius = os.path.join(sirius_path, 'sirius')

    initial_flags = os.environ.get('_JAVA_OPTIONS', '')
    if java_flags is not None:
        # append the flags to any existing options
        os.environ['_JAVA_OPTIONS'] = initial_flags + ' ' + java_flags
    cmdsir = ([sirius, '-o', artifact.get_path()] + parameters)

    stdout = os.path.join(str(artifact.path), 'stdout.txt')
    stderr = os.path.join(str(artifact.path), 'stderr.txt')

    run_command(cmdsir, stdout, stderr)

    if java_flags is not None:
        os.environ['_JAVA_OPTIONS'] = initial_flags

    return artifact


def compute_fingerprint_classes(sirius_path: str, features: MGFDirFmt,
                                ppm_max: int, profile: str,
                                tree_timeout: int = 1600,
                                maxmz: int = 600, n_jobs: int = 1,
                                num_candidates: int = 50,
                                database: List[Str] = ['all'],
                                ions_considered: List[Str] = ['[M+H]+'],
                                java_flags: str = None,
                                zodiac_threshold: float = 0.98,
                                fingerid_db: str = 'bio') -> SiriusDirFmt:
    '''Compute fragmentation trees for candidate molecular formulas using SIRIUS; Reranks molecular formula
    candidates generated using ZODIAC; Predict molecular fingerprints using CSI:FingerID; Predict compound
    classes using CANOPUS (ClassyFire and NPClassifier).

    Parameters
    ----------
    sirius_path : str
        Path to Sirius executable (without including the word sirius).
    features : MGFDirFmt
        MGF file for Sirius
    ppm_max : int
        allowed parts per million tolerance for decomposing masses
    profile: str
        configuration profile for mass-spec platform used
    tree_timeout : int, optional
        time for computation per fragmentation tree in seconds. 0 for an
        infinite amount of time
    maxmz : int, optional
        considers compounds with a precursor mz lower or equal to this
        value (int)
    n_jobs : int, optional
        Number of cpu cores to use. If not specified Sirius uses all available
        cores
    num_candidates: int, optional
        number of fragmentation trees to compute per feature
    database: str, optional
        search formulas in given database
    ions_considered : str, optional
        The ion type/adduct of the MS/MS data. You can also provide a
        comma separated list of adducts (default: '[M+H]+')
    java_flags : str, optional
        Setup additional flags for the Java virtual machine.
    zodiac_threshold: float, optional
        threshold filter for molecular formula re-ranking. Higher value
        recommended for less false positives (default: 0.98)
    fingerid_db: str, optional
        Search structure in given database. You can also provide a
        comma-separated list of databases (default: 'bio')

    Returns
    -------
    SiriusDirFmt
        Directory with computed fragmentation trees
    '''

    # qiime2 will check that the only possible modes are positive, negative or
    # auto

    params = ['-i', os.path.join(str(features.path), 'features.mgf'),
              '--maxmz', str(maxmz),
              '--processors', str(n_jobs),
              '--initial-compound-buffer', str(1),
              'formula',
              '--profile', str(profile),
              '--database', str(database),
              '--candidates', str(num_candidates),
              '--ions-considered', str(', '.join(ions_considered)),
              '--tree-timeout', str(tree_timeout),
              '--ppm-max', str(ppm_max),
              'zodiac',
              '--thresholdFilter', str(zodiac_threshold),
              'fingerprint',
              'structure',
              '--database', str(fingerid_db),
              'canopus',
              'write-summaries']

    return artifactory(sirius_path, params, java_flags, SiriusDirFmt)
