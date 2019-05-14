# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, q2-qemistree development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from q2_qemistree.plugin_setup import plugin as qemistree_plugin


class PluginSetupTests(unittest.TestCase):

    def test_plugin_setup(self):
        self.assertEqual(qemistree_plugin.name, 'qemistree')
