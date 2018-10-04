import subprocess
import os
import biom
import pandas as pd
import numpy as np
import shutil
import tempfile


def run_command(cmd, output_fp, sirpath, verbose=True):
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
    sirout : path to SIRIUS output folder

    Returns:
    ----------
    table : biom object
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
        raise RuntimeError('Fingerprint file is empty!')
    fingerids.index.name = '#featureID'
    npfid = np.asarray(fingerids)
    fptable = biom.table.Table(data=npfid, observation_ids=fingerids.index,
                               sample_ids=fingerids.columns)
    return fptable


def fingerprint(sirpath: str, features: str, ppmlim: int, instrument: str,
                nproc: int=1, nft: int=75, ftsec: int=1600,
                database: str='all', dbcsi: str='bio', mzlim: int=800,
                zodthresh: float=0.8, minconloc: int=15) -> biom.Table:
    '''
    This function generates and collates chemical fingerprints for mass-spec
    features in an experiment.

    Parameters
    ----------
    sirpath : path to SIRIUS binaries on user's computer
    features : path to MGF file for SIRIUS (str)
    ppmlim : parts per million tolerance (int)
    nft : number of fragmentation trees to compute per feature (int)
    ftsec : time for computation per fragmentation tree in seconds (int)
    database : database for SIRIUS (str)
    dbcsi : database for CSIFingerID (str)
    instrument : mass-spec platform used (str)
    mzlim : Maximum precursor to search (int)
    nproc : Number of processors used for computation (int)
    zodthresh : threshold filter for zodiac (float)
    minconloc : minimum local connections for zodiac (int)

    Returns:
    ----------
    table : biom object
        biom table containing mass-spec feature IDs (in rows) and molecular
        substructure IDs (in columns). Values are presence (1) or absence (0)
        of a particular substructure.
    '''

    tmpdir = tempfile.mkdtemp()

    if not os.path.exists(sirpath):
        raise OSError("SIRIUS could not be located")
    sirius = os.path.join(sirpath, 'sirius')
    if not os.path.exists(features):
        raise OSError("MGF file could not be located")

    tmpsir = os.path.join(tmpdir, 'tmpsir')
    cmdsir = [str(sirius), '--quiet',
              '--initial-compound-buffer', str(1),
              '--max-compound-buffer', str(32), '--profile', str(instrument),
              '--database', str(database), '--candidates',  str(nft),
              '--processors', str(nproc),
              '--auto-charge', '--trust-ion-prediction',
              '--maxmz', str(mzlim),
              '--tree-timeout', str(ftsec),
              '--ppm-max', str(ppmlim),
              '-o', str(tmpsir), str(features)]

    tmpzod = os.path.join(tmpdir, 'tmpzod')
    cmdzod = [str(sirius), '--zodiac', '--sirius', str(tmpsir),
              '-o', str(tmpzod),
              '--thresholdfilter', str(zodthresh),
              '--processors', str(nproc),
              '--minLocalConnections', str(minconloc),
              '--spectra', str(features)]

    tmpcsi = os.path.join(tmpdir, 'tmpcsi')
    cmdfid = [str(sirius), '--fingerid', '--fingerid-db', str(dbcsi),
              '--ppm-max', str(ppmlim), '-o', str(tmpcsi), str(tmpzod)]

    run_command(cmdsir, os.path.join(tmpdir, 'sirout'), sirpath)
    run_command(cmdzod, os.path.join(tmpdir, 'zodout'), sirpath)
    run_command(cmdfid, os.path.join(tmpdir, 'csiout'), sirpath)

    fptable = collatefp(tmpcsi)
    shutil.rmtree(tmpdir)

    return fptable
