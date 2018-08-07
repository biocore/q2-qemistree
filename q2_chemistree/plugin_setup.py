from qiime2.plugin import Plugin

plugin = Plugin(
    name='chemistree',
    version='0.0.0',
    website='https://github.com/biocore/q2-chemistree',
    package='q2_chemistree',
    description='Hierarchical orderings for mass spectrometry data',
    short_description='Plugin for exploring chemical diversity.',
)

def plugin_function():
    pass

plugin.visualizers.register_function(
    function=plugin_function,
    inputs={},
    parameters={},
    input_descriptions={},
    parameter_descriptions={},
    name='test function',
    description='a test to check installation'
)