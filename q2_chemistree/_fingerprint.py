# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import subprocess
import os
import biom
import pandas as pd
import numpy as np
import shutil
import tempfile

from ._semantics import MGFDirFmt


def run_command(cmd, output_fp, verbose=True):
    if verbose:
        print("Running external command line application. This may print "
              "messages to stdout and/or stderr.")
        print("The command being run is below. This command cannot "
              "be manually re-run as it will depend on temporary files that "
              "no longer exist.")
        print("\nCommand:", end=' ')
        print(" ".join(cmd), end='\n\n')

    with open(output_fp, 'w') as output_f:
        subprocess.run(cmd, stdout=output_f, check=True)


def collatefp(csiout):
    '''
    This function collates chemical fingerprints for mass-spec
    features in an experiment.

    Parameters
    ----------
    csiout : path to CSI:FingerID output folder

    Raises
    ------
    ValueError
        If ``fptable`` (collated fingerprint table) is empty

    Returns
    -------
    biom.Table
        biom table containing mass-spec feature IDs (in rows) and molecular
        substructure IDs (in columns). Values are presence (1) or absence (0)
        of a particular substructure.
    '''

    fpfoldrs = os.listdir(csiout)
    molfp = dict()
    for foldr in fpfoldrs:
        if os.path.isdir(os.path.join(csiout, foldr)):
            fidpath = os.path.join(csiout, foldr)
            fid = foldr.split('_')[-1]
            if 'fingerprints' in os.listdir(fidpath):
                fname = os.listdir(os.path.join(fidpath, 'fingerprints'))[0]
                with open(os.path.join(fidpath, 'fingerprints', fname)) as f:
                    fp = f.read().strip().split('\n')
                molfp[fid] = fp

    fingerids = pd.DataFrame.from_dict(molfp, orient='index')
    if fingerids.shape == (0, 0):
        raise ValueError('Fingerprint file is empty!')
    fingerids.index.name = '#featureID'
    npfid = np.asarray(fingerids)

    # biom requires that ids be strings
    fptable = biom.table.Table(data=npfid,
                               observation_ids=fingerids.index.astype(str),
                               sample_ids=fingerids.columns.astype(str))
    return fptable


def fingerprint(sirius_path: str, features: MGFDirFmt, ppm_max: int,
                profile: str, n_jobs: int = 1,
                num_candidates: int = 75, tree_timeout: int = 1600,
                database: str = 'all', fingerid_db: str = 'pubchem',
                maxmz: int = 600, zodiac_threshold: float = 0.95) -> biom.Table:
    '''
    This function generates and collates chemical fingerprints for mass-spec
    features in an experiment.

    Parameters
    ----------
    sirius_path : path to Sirius executable (str)
    features : MGF file for SIRIUS (str)
    ppm_max : allowed parts per million tolerance for decomposing masses (int)
    profile : configuration profile for mass-spec platform used (str)
    n_jobs : Number of cpu cores to use. If not specified Sirius uses
                 all available cores (int)
    num_candidates : number of fragmentation trees to compute per feature (int)
    tree_timeout : time for computation per fragmentation tree in seconds.
                   0 for an infinite amount of time (int)
    database : search formulas in given database (str)
    fingerid_db : search structure in given database (str)
    maxmz : considers compounds with a precursor mz lower or equal to
            this value (int)
    zodiac_threshold : threshold filter for molecular formula re-ranking.
                       Higher value recommended for
                       less false positives (float)

    Raises
    ------
    OSError:
        If ``sirius_path`` not found

    Returns
    -------
    biom.Table
        biom table containing mass-spec feature IDs (in rows) and molecular
        substructure IDs (in columns). Values are presence (1) or absence (0)
        of a particular substructure.
    '''
    if isinstance(features, MGFDirFmt):
        features = str(features.path) + '/features.mgf'

    tmpdir = tempfile.mkdtemp()

    if not os.path.exists(sirius_path):
        raise OSError("SIRIUS could not be located")
    sirius = os.path.join(sirius_path, 'sirius')

    tmpsir = os.path.join(tmpdir, 'tmpsir')
    cmdsir = [str(sirius), '--quiet',
              '--initial-compound-buffer', str(1),
              '--max-compound-buffer', str(32), '--profile', str(profile),
              '--database', str(database),
              '--candidates', str(num_candidates),
              '--processors', str(n_jobs),
              '--auto-charge', '--trust-ion-prediction',
              '--maxmz', str(maxmz),
              '--tree-timeout', str(tree_timeout),
              '--ppm-max', str(ppm_max),
              '-o', str(tmpsir), str(features)]

    tmpzod = os.path.join(tmpdir, 'tmpzod')
    cmdzod = [str(sirius), '--zodiac', '--sirius', str(tmpsir),
              '-o', str(tmpzod),
              '--thresholdfilter', str(zodiac_threshold),
              '--processors', str(n_jobs),
              '--spectra', str(features)]

    tmpcsi = os.path.join(tmpdir, 'tmpcsi')
    cmdfid = [str(sirius), '--processors', str(n_jobs), '--fingerid',
              '--fingerid-db', str(fingerid_db), '--ppm-max', str(ppm_max),
              '-o', str(tmpcsi), str(tmpzod)]

    run_command(cmdsir, os.path.join(tmpdir, 'sirout'))
    run_command(cmdzod, os.path.join(tmpdir, 'zodout'))
    run_command(cmdfid, os.path.join(tmpdir, 'csiout'))

    fptable = collatefp(tmpcsi)
    shutil.rmtree(tmpdir)
    return fptable
