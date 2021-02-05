import os
import shutil
import unittest

import dicom_utils as du

class TestDicomUtils(unittest.TestCase):  

    def test_anon_data(self):
        dicom_path = r"./tests/data/test_ct_abdomen"

    def test_integration_anonymize(self):
        """
        Integration test of DICOM anonymization method
        """

        test_dicom_path = r"./tests/data/test_ct_abdomen"
        test_output_path = r"./tests/output"

        if os.path.exists(test_output_path):
            shutil.rmtree(test_output_path)
        os.mkdir(test_output_path)

        # Preserve the patient age, birth date, and gender. This may be desired for certain research, etc.
        preserve_list = ("PatientAge", "PatientBirthDate", "PatientSex") 
        anon_dicom_path = os.path.join(test_output_path, 'anon_dicom')

        # Use all available options for anonymize:
        # Specify a new MRN, Accession, patient name - also potentially useful for research
        # Specify a new series description. There are actually 3 series in the test_ct_abdomen folder, 
        # so this will test the zip method's ability to still put them in different zip files.
        du.anonymize(test_dicom_path, anon_dicom_path, preserve_list=preserve_list, MRN='0000001', Accession='0000001', PatientName='LAST^FIRST', SeriesDescription='Anonymized CT Abdomen', verbose=True)

        # Verify values in 1st DICOM file
        ds = du.get_dicom_dataset(anon_dicom_path)
        self.assertEqual(ds.PatientAge, '037Y')
        self.assertEqual(ds.PatientBirthDate, '19791207')
        self.assertEqual(ds.PatientSex, 'M')
        self.assertEqual(ds.PatientID, '0000001')
        self.assertEqual(ds.AccessionNumber, '0000001')
        self.assertEqual(ds.PatientName, 'LAST^FIRST')
        self.assertEqual(ds.SeriesDescription, 'Anonymized CT Abdomen')
        self.assertEqual(ds.ReferringPhysicianName, 'anon')

        # Output metadata for the DICOM - this will require some manual verification
        anon_metadata_path = os.path.join(test_output_path, 'anon_metadata.txt')
        du.output_metadata(anon_dicom_path, anon_metadata_path, verbose=True)

        # Zip the anonymized DICOM. 3 different zip files should be created, one for each series
        zip_output_path = os.path.join(test_output_path, 'zip')
        du.zip_dicom(anon_dicom_path, zip_output_path, verbose=True)
        zipfiles = os.listdir(zip_output_path)
        self.assertEqual(zipfiles[0], 'Anonymized CT Abdomen.zip')
        self.assertEqual(zipfiles[1], 'Anonymized CT Abdomen 2.zip')
        self.assertEqual(zipfiles[2], 'Anonymized CT Abdomen 3.zip')

        shutil.rmtree(test_output_path)

    def test_unit_dcm2jpg(self):
        """
        Unit test of the dcm2jpg method
        """

        test_dicom_path = r"./tests/data/test_ct_abdomen"
        test_output_path = r"./tests/output"

        if os.path.exists(test_output_path):
            shutil.rmtree(test_output_path)
        os.mkdir(test_output_path)

        du.dcm2jpg(test_dicom_path, test_output_path, window_width=2000, window_level=400, verbose=True)

        # In this case, there are 3 series in the test_dicom_path.
        # Only the first series of 102 images should be loaded for conversion
        self.assertEqual(len(os.listdir(test_output_path)), 102)

        # Files should be in jpg format - visually confirm the expected window and level
        self.assertTrue(os.listdir(test_output_path)[0].endswith('.jpg'))

        shutil.rmtree(test_output_path)

if __name__ == '__main__':
    unittest.main()    

