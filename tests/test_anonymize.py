import unittest

import dicom_utils as du

class TestAnonymize(unittest.TestCase):

    def test_integration_anonymize(self):
        """
        Integration test of DICOM anonymization method
        """

        test_dicom_path = r"./tests/data/test_ct_abdomen"
        test_output_path = r"./tests/data/test_output"

        preserve_list = ("PatientAge", "PatientBirthDate", "PatientSex")
        du.anonymize(test_dicom_path, test_output_path, preserve_list=preserve_list, MRN='0000001', Accession='0000001', PatientName='LAST^FIRST', SeriesDescription='Anonymized CT Abdomen', verbose=False)

if __name__ == '__main__':
    unittest.main()        

