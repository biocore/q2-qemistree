# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from ._fingerprint import fingerprint
from ._hierarchy import make_hierarchy
from ._semantics import MassSpectrometryFeatures, MGFDirFmt

from qiime2.plugin import Plugin, Str, Range, Choices, Float, Int
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.tree import Phylogeny, Rooted

plugin = Plugin(
    name='chemistree',
    version='0.0.0',
    website='https://github.com/biocore/q2-chemistree',
    package='q2_chemistree',
    description='Hierarchical orderings for mass spectrometry data',
    short_description='Plugin for exploring chemical diversity.',
)

# type registration
plugin.register_views(MGFDirFmt)
plugin.register_semantic_types(MassSpectrometryFeatures)
plugin.register_semantic_type_to_format(MassSpectrometryFeatures,
                                        artifact_format=MGFDirFmt)

PARAMS = {
    'database': Str % Choices(['all', 'PubChem', 'bio']),
    'sirpath': Str,
    'ppmlim': Int,
    'instrument': Str,
    'nproc': Int,
    'nft': Int,
    'ftsec': Int,
    'dbcsi': Str,
    'mzlim': Int,
    'zodthresh': Float,
    'minconloc': Int,
}

PARAMS_DESC = {
    'database': 'Database for Sirius',
    'sirpath': 'Path to the Sirius executable',
    'ppmlim': 'Parts per million tolerance for ions',
    'instrument': 'Mass-spec platform used',
    'nproc': 'Number of processors used by Sirius',
    'nft': 'Number of fragmentation trees to compute per feature',
    'ftsec': 'Time for computation per fragmentation tree in seconds',
    'dbcsi': 'Database for CSIFingerID',
    'mzlim': 'Maximum precursor to search',
    'zodthresh': 'Threshold filter for zodiac',
    'minconloc': 'Minimum local connections for zodiac',
}

# method registration
plugin.methods.register_function(
    function=fingerprint,
    name='Fingerprint MS data',
    description='Create a contingency table of molecular substructures for a '
                'set of mass spectrometry features',
    inputs={'features': MassSpectrometryFeatures},
    parameters=PARAMS,
    input_descriptions={'features': 'List of MS1 ions and corresponding '
                                    'MS2 ions for each MS1.'},
    parameter_descriptions=PARAMS_DESC,
    outputs=[('collated_fingerprints', FeatureTable[Frequency])],
    output_descriptions={'collated_fingerprints': 'Contingency table of the '
                                                  'probabilities of '
                                                  'molecular substructures '
                                                  'within each feature'}
)

plugin.methods.register_function(
    function=make_hierarchy,
    name='Create a molecular tree',
    description='Build a phylogeny based on molecular substructures',
    inputs={'collated_fingerprints': FeatureTable[Frequency]},
    parameters={'threshold': Float % Range(0, 1, inclusive_end=True)},
    input_descriptions={'collated_fingerprints': 'Contingency table of the '
                                                 'probabilities of '
                                                 'molecular substructures '
                                                 'within each feature'},
    parameter_descriptions={'threshold': 'Probability threshold below which a'
                                         ' substructure is considered absent'},
    outputs=[('tree', Phylogeny[Rooted])],
    output_descriptions={'tree': 'Tree of relatedness between mass '
                                 'spectrometry features based on the chemical '
                                 'substructures within those features'}
)
