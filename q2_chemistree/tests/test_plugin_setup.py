# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, q2-chemistree development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from q2_chemistree.plugin_setup import plugin as chemistree_plugin


class PluginSetupTests(unittest.TestCase):

    def test_plugin_setup(self):
        self.assertEqual(chemistree_plugin.name, 'chemistree')
