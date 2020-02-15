# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import q2_qemistree
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
from ._plot import plot

from qiime2.plugin import (Plugin, Str, Range, Choices, Float, Int, Bool, List,
                           Citations)
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.tree import Phylogeny, Rooted


citations = Citations.load('citations.bib', package='q2_qemistree')

plugin = Plugin(
    name='qemistree',
    version=q2_qemistree.__version__,
    website='https://github.com/biocore/q2-qemistree',
    package='q2_qemistree',
    description='Hierarchical orderings for mass spectrometry data',
    short_description='Plugin for exploring chemical diversity.',
    citations=citations,
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
                                                'by Sirius'},
    citations=[citations['duhrkop2015sirius']]
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
                                               'using Zodiac'},
    citations=[citations['duhrkop2015sirius']]
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
                                                   'CSI:FingerID'},
    citations=[citations['duhrkop2015sirius']]
)

plugin.methods.register_function(
    function=make_hierarchy,
    name='Create a molecular tree',
    description='Build a phylogeny based on molecular substructures',
    inputs={'csi_results': List[CSIFolder],
            'feature_tables': List[FeatureTable[Frequency]],
            'ms2_matches': List[FeatureData[Molecules]]},
    parameters={'qc_properties': Bool,
                'metric': Str % Choices(['euclidean', 'jaccard'])},
    input_descriptions={'csi_results': 'one or more CSI:FingerID '
                                       'output folders',
                        'feature_tables': 'one or more feature tables with '
                                          'mass-spec feature intensity '
                                          'per sample',
                        'ms2_matches': 'one or more tables with MS/MS library '
                                       'match for mass-spec features'},
    parameter_descriptions={'qc_properties': 'filters molecular properties to '
                                             'retain PUBCHEM fingerprints',
                            'metric': 'metric for hierarchical clustering of '
                                      'fingerprints. If the Jaccard metric is'
                                      ' selected, molecular fingerprints are '
                                      'first binarized (probabilities above '
                                      '0.5 are True, and False otherwise).'},
    outputs=[('tree', Phylogeny[Rooted]),
             ('feature_table', FeatureTable[Frequency]),
             ('feature_data', FeatureData[Molecules])],
    output_descriptions={'tree': 'Tree of relatedness between mass '
                                 'spectrometry features based on the chemical '
                                 'substructures within those features',
                         'feature_table': 'filtered feature table '
                                          'that contains only the '
                                          'features present in '
                                          'the tree',
                         'feature_data': 'mapping of unique feature '
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
                                                    'spec feature'},
    citations=[citations['djoumbou2016classyfire']]
)

plugin.methods.register_function(
    function=prune_hierarchy,
    name='Prune hierarchy of molecules',
    description='Removes the tips of the tree based on feature data',
    inputs={'feature_data': FeatureData[Molecules],
            'tree': Phylogeny[Rooted]},
    parameters={'column': Str},
    input_descriptions={'feature_data': 'Feature data table with '
                                        'molecules to keep',
                        'tree': 'Tree of relatedness of molecules.'},
    parameter_descriptions={'column': 'A column in feature data table. '
                                      'Features with missing values in this '
                                      'column will be removed from the tree. '
                                      'If no column name is specified then '
                                      'the tree will be pruned to only '
                                      'contain features in the feature data.'},
    outputs=[('pruned_tree', Phylogeny[Rooted])],
    output_descriptions={'pruned_tree': 'Pruned tree of molecules with '
                                        'tips that are in feature data'}
)

plugin.visualizers.register_function(
    function=plot,
    name='Generate an iTOL bar chart annotation file',
    description=('Calculate mean frequency per category per sample and render '
                 'as bar height.'),
    inputs={'grouped_table': FeatureTable[Frequency],
            'tree': Phylogeny[Rooted],
            'feature_metadata': FeatureData[Molecules]
            },
    parameters={
        'category': Str,
        'parent_mz': Str,
        'color_palette': Str % Choices(['Pastel1', 'Pastel2', 'Paired',
                                        'Accent', 'Dark2', 'Set1', 'Set2',
                                        'Set3', 'tab10', 'tab20', 'tab20b',
                                        'tab20c', 'Greys', 'Purples', 'Blues',
                                        'Greens', 'Oranges', 'Reds', 'YlOrBr',
                                        'YlOrRd', 'OrRd', 'PuRd', 'RdPu',
                                        'BuPu', 'GnBu', 'PuBu', 'YlGnBu',
                                        'PuBuGn', 'BuGn', 'YlGn']),
        'ms2_label': Bool,
                },
    input_descriptions={'grouped_table': 'Feature table of samples '
                                         'grouped by categories',
                        'tree': 'Phenetic tree',
                        'feature_metadata': 'Feature metadata'
                        },
    parameter_descriptions={
        'category': 'The feature data column used to color and label the tips',
        'color_palette': 'The color palette to use for coloring tips. '
                         'For examples, see: https://matplotlib.org/'
                         'tutorials/colors/colormaps.html',
        'ms2_label': 'Whether to label the tips with the MS2 value',
        'parent_mz': 'If the feature is unclassified, label the tips using '
                     'this column\'s value'
    },
    citations=[citations['letunic2019itol']])

importlib.import_module('q2_qemistree._transformer')
