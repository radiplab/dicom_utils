import numpy as np
import os
import pydicom
import shutil
import time

general_keep_list = [
                        'AccessionNumber',
                        'BitsAllocated',
                        'BitsStored',
                        'BodyPartExamined',
                        'Columns',
                        'FrameOfReferenceUID',
                        'HighBit',
                        'HighRRValue',
                        'ImageOrientationPatient',
                        'ImagePositionPatient',
                        'ImageType',
                        'InstanceNumber',
                        'Laterality',
                        'LossyImageCompression',
                        'LowRRValue',
                        'Manufacturer',
                        'ManufacturerModelName'
                        'MediaStorageSOPClassUID',
                        'MediaStorageSOPInstanceUID',
                        'Modality',
                        'PatientAge',
                        'PatientBirthDate',
                        'PatientID',
                        'PatientName',
                        'PatientPosition',
                        'PatientSex',
                        'PhotometricInterpretation',
                        'PixelBandwidth',
                        'PixelIntensityRelationship',
                        'PixelPaddingValue',
                        'PixelRepresentation',
                        'PixelSpacing',
                        'PositionReferenceIndicator',
                        'ReconstructionDiameter',
                        'ReferringPhysicianName',
                        'RescaleIntercept',
                        'RescaleSlope',
                        'RescaleType',
                        'Rows',
                        'ScanOptions',
                        'SeriesDescription',
                        'SeriesInstanceUID',
                        'SeriesNumber',
                        'SequenceVariant',
                        'SliceLocation',
                        'SliceThickness',
                        'SOPClassUID',
                        'SOPInstanceUID',
                        'SpacingBetweenSlices',
                        'SpecificCharacterSet',
                        'StudyDate',
                        'StudyDescription',
                        'StudyID',
                        'StudyInstanceUID',
                        'StudyTime',
                        'TemporalPositionIdentifier',
                        'WindowCenter',
                        'WindowWidth'
                    ]

mr_keep_list = [
                    'AcquisitionDuration',
                    'AcquisitionMatrix',
                    'AcquisitionNumber',
                    'B1rms',
                    'dBdt',
                    'DiffusionBValue',
                    'DiffusionGradientOrientation',
                    'EchoNumbers',
                    'EchoTime',
                    'EchoTrainLength',
                    'FlipAngle',
                    'ImagingFrequency',
                    'InPlanePhaseEncodingDirection',
                    'IntervalsAcquired',
                    'IntervalsRejected',
                    'LowRRValue',
                    'MagneticFieldStrength',
                    'MRAcquisitionType',
                    'NumberOfAverages',
                    'NumberOfPhaseEncodingSteps',
                    'NumberOfTemporalPositions',
                    'PercentPhaseFieldOfView',
                    'PercentSampling',
                    'PixelBandwidth',
                    'PixelPaddingValue',
                    'PixelRepresentation',
                    'PixelSpacing',
                    'RepetitionTime',
                    'Rows',
                    'SamplesPerPixel',
                    'SAR'
                ]

angio_keep_list = [
                    'CineRate',
                    'DerivationDescription',
                    'FrameIncrementPointer',
                    'FrameTimeVector',
                    'ImagerPixelSpacing',
                    'NumberOfFrames',
                    'PositionerMotion',
                    'PositionerPrimaryAngle',
                    'PositionerSecondaryAngle',
                    'RecommendedDisplayFrameRate',
                    'ShutterLeftVerticalEdge',
                    'ShutterLowerHorizontalEdge',
                    'ShutterRightVerticalEdge',
                    'ShutterShape',
                    'ShutterUpperHorizontalEdge',
                    'TableMotion'
                ]

def anonymize(dicom_path, output_path, preserve_list=list(), MRN=None, Accession=None, PatientName=None, SeriesDescription=None, verbose=False):
    """ Anonymizes all DICOM files within the provided path, saving them to the 
    specified output_path. The original DICOM files are NOT modified. All
    contents in the output_path will be deleted.

    Parameters:
        dicom_path (str):The full path to the directory containing DICOM files
         - these will not be modified
        output_path (str): The full path to the directory that will contain the output folder. Any existing files in this directory will not be deleted
        preserve_list (list): List of dataset keywords to not be altered in the anonymization process.
            Can only be applied to a small number of keywords. Otherwise alter lists above.
    """
    start_time = None
    if verbose:
        start_time = time.time()
        print("Beginning DICOM anonymization..", end = '')

    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.mkdir(output_path)

    # Build the output directory
    output_dir = ''
    if MRN is not None:
        output_dir = output_dir + 'MRN' + MRN
    if Accession is not None:
        if output_dir != '':
            output_dir = output_dir + '-'
        output_dir = output_dir + 'Acc' + Accession

    if os.path.exists(os.path.join(output_path, output_dir)):
        shutil.rmtree(os.path.join(output_path, output_dir))
    os.mkdir(os.path.join(output_path, output_dir))
    output_path = os.path.join(output_path, output_dir)

    dicom_files = os.listdir(dicom_path)
    for dicom_file in dicom_files:
        if os.path.isfile(os.path.join(dicom_path, dicom_file)) and not dicom_file.startswith('KO') and not dicom_file.startswith('.DS_Store'):
            shutil.copy(
                os.path.join(dicom_path, dicom_file),
                os.path.join(output_path, dicom_file)
                )

    dicom_files = os.listdir(output_path)
    for dicom_file in dicom_files:
        if os.path.isfile(os.path.join(output_path, dicom_file)) and not dicom_file.startswith('KO') and not dicom_file.startswith('.DS_Store'):
            if not dicom_file.endswith('.dcm'):
                shutil.move(
                    os.path.join(output_path, dicom_file),
                    os.path.join(output_path, dicom_file + '.dcm')
                )

    dicom_files = os.listdir(output_path)
    for dicom_file in dicom_files:
        if dicom_file.endswith('.dcm') and not dicom_file.startswith('KO') and not dicom_file.startswith('.DS_Store'):
            #print(" - Processing " + dicom_file)
            ds = pydicom.filereader.dcmread(
                os.path.join(output_path, dicom_file)
            )

            for data_element in ds:
                if "pixel data" not in str(data_element.name).lower() \
                        and "overlay data" \
                        not in str(data_element.name).lower():
                    delete_de = True

                    # Data elements derived from the following in dwv.js:
                    # dwv.dicom.DicomParser.prototype.parse

                    keep_list = general_keep_list + \
                        mr_keep_list + \
                        angio_keep_list

                    for keep in keep_list:
                        if keep in str(data_element.keyword):
                            delete_de = False

                    if "AccessionNumber" in str(data_element.keyword):
                        if Accession is not None:
                            data_element.value = str(Accession)
                        else:
                            data_element.value = "000000"
                    if "PatientAge" in str(data_element.keyword) and "PatientAge" not in preserve_list:
                        data_element.value = "18"
                    if "PatientBirthDate" in str(data_element.keyword) and "PatientBirthDate" not in preserve_list:
                        data_element.value = "01011900"
                    if "PatientID" in str(data_element.keyword):
                        if MRN is not None:
                            data_element.value = str(MRN)
                        else:
                            data_element.value = "000000"
                    if "PatientName" in str(data_element.keyword):
                        if PatientName is not None:
                            data_element.value = PatientName
                        else:
                            data_element.value = "anon"
                    if "PatientSex" in str(data_element.keyword) and "PatientSex" not in preserve_list:
                        data_element.value = "F"
                    if "ReferringPhysicianName" in str(data_element.keyword):
                        data_element.value = "anon"
                    if "SeriesDescription" in str(data_element.keyword):
                        if SeriesDescription is not None:
                            data_element.value = SeriesDescription                    
                    if "StudyDate" in str(data_element.keyword) and "StudyDate" not in preserve_list:
                        data_element.value = "01011900"
                    if "StudyID" in str(data_element.keyword):
                        data_element.value = "000000"
                    if delete_de:
                        del ds[data_element.tag]

            ds.save_as(os.path.join(output_path, dicom_file))
    if verbose:
        print('.done - time = ' + str(np.around((time.time() - start_time), decimals=1)) + 'sec')
    return output_path