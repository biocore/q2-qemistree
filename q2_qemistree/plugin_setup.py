# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib
from ._fingerprint import (compute_fragmentation_trees,
                           rerank_molecular_formulas,
                           predict_fingerprints)
from ._hierarchy import make_hierarchy
from ._prune_hierarchy import prune_hierarchy
from ._classyfire import get_classyfire_taxonomy
from ._semantics import (MassSpectrometryFeatures, MGFDirFmt,
                         SiriusFolder, SiriusDirFmt,
                         ZodiacFolder, ZodiacDirFmt,
                         CSIFolder, CSIDirFmt,
                         FeatureData, TSVMoleculesFormat, Molecules)

from qiime2.plugin import Plugin, Str, Range, Choices, Float, Int, Bool, List
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.tree import Phylogeny, Rooted

plugin = Plugin(
    name='qemistree',
    version='0.0.0',
    website='https://github.com/biocore/q2-qemistree',
    package='q2_qemistree',
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

plugin.register_views(TSVMoleculesFormat)
plugin.register_semantic_types(Molecules)
plugin.register_semantic_type_to_format(FeatureData[Molecules],
                                        artifact_format=TSVMoleculesFormat)

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
        'n_jobs', 'num_candidates', 'database', 'profile', 'java_flags',
        'ionization_mode']
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
    inputs={'csi_results': List[CSIFolder],
            'feature_tables': List[FeatureTable[Frequency]]},
    parameters={'qc_properties': Bool},
    input_descriptions={'csi_results': 'one or more CSI:FingerID '
                                       'output folders',
                        'feature_tables': 'one or more feature tables with '
                                          'mass-spec feature intensity '
                                          'per sample'},
    parameter_descriptions={'qc_properties': 'filters molecular properties to '
                                             'retain PUBCHEM fingerprints'},
    outputs=[('tree', Phylogeny[Rooted]),
             ('merged_feature_table', FeatureTable[Frequency]),
             ('merged_fingerprints', FeatureData[Molecules]),
             ('merged_feature_data', FeatureData[Molecules])],
    output_descriptions={'tree': 'Tree of relatedness between mass '
                                 'spectrometry features based on the chemical '
                                 'substructures within those features',
                         'merged_feature_table': 'filtered feature table '
                                                 'that contains only the '
                                                 'features present in '
                                                 'the tree',
                         'merged_fingerprints': 'table of molecular '
                                                'fingerprints indexed by the '
                                                'MD5 hash of the fingerprint '
                                                'vectors predicted '
                                                'by CSI:FingerID',
                         'merged_feature_data': 'mapping of unique feature '
                                                'identifiers in input '
                                                'feature tables to MD5 hash '
                                                'of feature fingerprints'}
)

plugin.methods.register_function(
    function=get_classyfire_taxonomy,
    name='Generate Classyfire annotations',
    description='Predicts chemical taxonomy based on molecule structures',
    inputs={'feature_data': FeatureData[Molecules]},
    parameters={},
    input_descriptions={'feature_data': 'Feature data table that maps MD5 '
                                        'hash of mass-spec features to their '
                                        'structural annotations (SMILES)'},
    parameter_descriptions={},
    outputs=[('classified_feature_data', FeatureData[Molecules])],
    output_descriptions={'classified_feature_data': 'Feature data table that '
                                                    'contains Classyfire '
                                                    'annotations per mass-'
                                                    'spec feature'}
)

plugin.methods.register_function(
    function=prune_hierarchy,
    name='Prune hierarchy of molecules',
    description='Removes non-annotated tree tips based on feature data',
    inputs={'feature_data': FeatureData[Molecules],
            'tree': Phylogeny[Rooted]},
    parameters={'column': Str % Choices(['smiles', 'kingdom', 'superclass',
                                         'class', 'subclass',
                                         'direct_parent'])},
    input_descriptions={'feature_data': 'Feature data table with Classyfire '
                                        'annotations and/or SMILES.',
                        'tree': 'Tree of relatedness of molecules.'},
    parameter_descriptions={'column': 'Features with missing values in this '
                                      'column will be removed from the tree'},
    outputs=[('pruned_tree', Phylogeny[Rooted])],
    output_descriptions={'pruned_tree': 'Pruned tree of molecules with '
                                        'tips that are in feature data'}
)

importlib.import_module('q2_qemistree._transformer')
