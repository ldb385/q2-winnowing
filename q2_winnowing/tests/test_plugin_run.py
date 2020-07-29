
import os

from unittest import TestCase, main as unittest_main
from q2_winnowing.winnow import process

class PluginRunTests( TestCase ):

    package = 'q2_winnowing.tests'
    testing_data = f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/test_run/test_in_data.biom",

    def test_plugin_setup(self):



if __name__ == '__main__':
    unittest_main()