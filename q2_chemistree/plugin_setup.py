# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from ._fingerprint import (compute_fragmentation_trees,
                           rerank_molecular_formulas,
                           predict_fingerprints)
from ._hierarchy import make_hierarchy
from ._match import match_table
from ._collate_fingerprint import collate_fingerprint
from ._semantics import (MassSpectrometryFeatures, MGFDirFmt,
                         SiriusFolder, SiriusDirFmt,
                         ZodiacFolder, ZodiacDirFmt,
                         CSIFolder, CSIDirFmt)

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

plugin.register_views(SiriusDirFmt)
plugin.register_semantic_types(SiriusFolder)
plugin.register_semantic_type_to_format(SiriusFolder,
                                        artifact_format=SiriusDirFmt)

plugin.register_views(ZodiacDirFmt)
plugin.register_semantic_types(ZodiacFolder)
plugin.register_semantic_type_to_format(ZodiacFolder,
                                        artifact_format=ZodiacDirFmt)

plugin.register_views(CSIDirFmt)
plugin.register_semantic_types(CSIFolder)
plugin.register_semantic_type_to_format(CSIFolder,
                                        artifact_format=CSIDirFmt)

PARAMS = {
    'ionization_mode': Str % Choices(['positive', 'negative', 'auto']),
    'database': Str % Choices(['all', 'pubchem']),
    'sirius_path': Str,
    'profile': Str % Choices(['qtof', 'orbitrap', 'fticr']),
    'fingerid_db': Str % Choices(['all', 'pubchem', 'bio', 'kegg', 'hmdb']),
    'ppm_max': Int % Range(0, 30, inclusive_end=True),
    'n_jobs': Int % Range(1, None),
    'num_candidates': Int % Range(5, 100, inclusive_end=True),
    'tree_timeout': Int % Range(600, 3000, inclusive_end=True),
    'maxmz': Int % Range(100, 850, inclusive_end=True),
    'zodiac_threshold': Float % Range(0, 1, inclusive_end=True),
    'java_flags': Str
}

PARAMS_DESC = {
    'ionization_mode': 'Ionization mode for mass spectrometry',
    'database': 'search formulas in given database',
    'sirius_path': 'path to Sirius executable',
    'ppm_max': 'allowed parts per million tolerance for decomposing masses',
    'profile': 'configuration profile for mass-spec platform used',
    'n_jobs': 'Number of cpu cores to use',
    'num_candidates': 'number of fragmentation trees to compute per feature',
    'tree_timeout': 'time for computation per fragmentation tree in seconds',
    'fingerid_db': 'search structure in given database',
    'maxmz': 'consider compounds with a precursor mz lower or equal to this',
    'zodiac_threshold': 'threshold filter for molecular formula re-ranking',
    'java_flags': 'Setup additional flags for the Java virtual machine. '
                  'For Sirius it is recommended that you modify the initial '
                  'and maximum heap size. For example to set an initial and '
                  'maximum heap size of 16GB and 64GB (respectively) specify '
                  '"-Xms16G -Xmx64G". Note that the quotes are important.'
}

# method registration
keys = ['sirius_path', 'features', 'ppm_max', 'tree_timeout', 'maxmz',
        'n_jobs', 'num_candidates', 'database', 'profile', 'java_flags']
plugin.methods.register_function(
    function=compute_fragmentation_trees,
    name='Compute fragmentation trees for candidate molecular formulas',
    description='Use Sirius to compute fragmentation trees',
    inputs={'features': MassSpectrometryFeatures},
    parameters={k: v for k, v in PARAMS.items() if k in keys},
    input_descriptions={'features': 'List of MS1 ions and corresponding '
                                    'MS2 ions for each MS1.'},
    parameter_descriptions={k: v
                            for k, v in PARAMS_DESC.items() if k in keys},
    outputs=[('fragmentation_trees', SiriusFolder)],
    output_descriptions={'fragmentation_trees': 'fragmentation trees '
                                                'computed per feature '
                                                'by Sirius'}
)

keys = ['sirius_path', 'zodiac_threshold', 'n_jobs', 'java_flags']
plugin.methods.register_function(
    function=rerank_molecular_formulas,
    name='Reranks candidate molecular formulas',
    description='Use Zodiac to rerank candidate molecular formulas',
    inputs={'features': MassSpectrometryFeatures,
            'fragmentation_trees': SiriusFolder},
    parameters={k: v for k, v in PARAMS.items() if k in keys},
    input_descriptions={'features': 'List of MS1 ions and corresponding '
                                    'MS2 ions for each MS1.'},
    parameter_descriptions={k: v
                            for k, v in PARAMS_DESC.items() if k in keys},
    outputs=[('molecular_formulas', ZodiacFolder)],
    output_descriptions={'molecular_formulas': 'Top scored molecular formula '
                                               'per feature after reranking'
                                               'using Zodiac'}
)

keys = ['sirius_path', 'ppm_max', 'n_jobs', 'fingerid_db', 'java_flags']
plugin.methods.register_function(
    function=predict_fingerprints,
    name='Predict fingerprints for molecular formulas',
    description='Use CSI:FingerID to predict molecular formulas',
    inputs={'molecular_formulas': ZodiacFolder},
    parameters={k: v for k, v in PARAMS.items() if k in keys},
    input_descriptions={
        'molecular_formulas': 'The output from running Zodiac'},
    parameter_descriptions={k: v
                            for k, v in PARAMS_DESC.items() if k in keys},
    outputs=[('predicted_fingerprints', CSIFolder)],
    output_descriptions={'predicted_fingerprints': 'Predicted substructures '
                                                   'per feature using '
                                                   'CSI:FingerID'}
)

plugin.methods.register_function(
    function=make_hierarchy,
    name='Create a molecular tree',
    description='Build a phylogeny based on molecular substructures',
    inputs={'collated_fingerprints': FeatureTable[Frequency]},
    parameters={'prob_threshold': Float % Range(0, 1, inclusive_end=True)},
    input_descriptions={'collated_fingerprints': 'Contingency table of the '
                                                 'probabilities of '
                                                 'molecular substructures '
                                                 'within each feature'},
    parameter_descriptions={'prob_threshold': 'Probability threshold below '
                                              'which a substructure is '
                                              'considered absent.'},
    outputs=[('tree', Phylogeny[Rooted])],
    output_descriptions={'tree': 'Tree of relatedness between mass '
                                 'spectrometry features based on the chemical '
                                 'substructures within those features'}
)

plugin.methods.register_function(
    function=match_table,
    name='Match feature table to tree tips',
    description='Filters feature table to to match tree tips',
    inputs={'tree': Phylogeny[Rooted],
            'feature_table': FeatureTable[Frequency]},
    input_descriptions={'tree': 'Phylogenetic tree with the features that will'
                                ' be retained on the feature table',
                        'feature_table': 'Feature table that will be filtered '
                                         'based on the features of the '
                                         'phylogenetic tree'},
    parameters={},
    outputs=[('filtered_feature_table', FeatureTable[Frequency])],
    output_descriptions={'filtered_feature_table': 'filtered feature table '
                                                   'that contains only the '
                                                   'features present in '
                                                   'the tree'}
)

plugin.methods.register_function(
    function=collate_fingerprint,
    name='Collate fingerprints into a table',
    description='Collate fingerprints predicted by CSI:FingerID',
    inputs={'csi_result': CSIFolder},
    input_descriptions={'csi_result': 'CSI:FingerID output folder'},
    parameters={},
    parameter_descriptions={},
    outputs=[('collated_fingerprints', FeatureTable[Frequency])],
    output_descriptions={'collated_fingerprints': 'Contingency table of the '
                                                  'probabilities of '
                                                  'molecular substructures '
                                                  'within each feature'}
)
