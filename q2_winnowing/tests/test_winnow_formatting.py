
import pkg_resources
import unittest
import uuid

import qiime2.core.archive as archive
from qiime2.sdk import Artifact

from q2_winnowing import Winnowed


class WinnowedFormatTests(unittest.TestCase):
    # test class layout referenced from:
    #   https://github.com/qiime2-graveyard/q2-dummy-types/blob/master/q2_dummy_types/tests/test_mapping.py

    def test_data_import(self):
        fp = pkg_resources.resource_filename(
            'q2_winnowing.tests', 'sample_data/test_in_featureData.tsv')

        # `Artifact.import_data` copies `test_in_featureData.tsv` into the artifact after
        # performing validation on the file.
        artifact = Artifact.import_data(Winnowed, fp)

        self.assertEqual(artifact.type, Winnowed)
        self.assertIsInstance(artifact.uuid, uuid.UUID)

    def test_reader_transformer(self):


    def test_writer_transformer(self):



if __name__ == "__main__":
    unittest.main()