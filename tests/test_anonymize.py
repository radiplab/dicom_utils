import unittest

import dicom-utils as du

class TestAnonymize(unittest.TestCase):

    def test_integration_anonymize(self):
        """
        Integration test of DICOM anonymization method
        """

        test_dicom_path = r".tests/data/test-ct-abdomen"
        test_output_path = r".tests/data/test-output"

        

