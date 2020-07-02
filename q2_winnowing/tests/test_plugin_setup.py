
from unittest import TestCase, main as unittest_main

from q2_winnowing.plugin_setup import plugin as winnowing_plugin

class PluginSetupTests( TestCase ):

    def test_plugin_setup(self):
        self.assertEqual( winnowing_plugin.name, "winnowing")


if __name__ == '__main__':
    unittest_main()