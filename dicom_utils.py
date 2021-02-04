import numpy as np
import os
import pydicom
import shutil
import time
import zipfile
from zipfile import ZipFile

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
        print("Beginning anonymize of DICOM in " + dicom_path + " to " + output_path + "..", end = '')

    # Create the output directory
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.mkdir(output_path)

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

def filter_series_name(series):
    """ Removes special characters from a series name

    Parameters:
        series (str): Series name to be filtered

    Returns:
        (str): Filtered series name
    """    
    series = str(series)
    series = series.replace("/","")
    series = series.replace("\\","")
    series = series.replace(":","")
    series = series.replace("*","")
    series = series.replace("?","")
    series = series.replace("\"","")
    series = series.replace("<","")
    series = series.replace(">","")
    series = series.replace("|","")
    series = series.replace(".","")
    series = series.replace("+","")
    # For some bizarre reason it doesn't work with "fps" in the filename
    series = series.replace("fps", "")
    
    # Remove leading and trailing spaces
    series = series.lstrip()
    series = series.rstrip()  

    return series

def get_dicom_dataset(dicom_path):
    """ Returns a pydicom dataset given a path to a folder containing dicom files (all files should be from the same study). 
    Returns None if there is no valid DICOM file in the provided path

    Parameters:
        dicom_path (str):The full path to the directory containing DICOM files
    """
    dicom_dataset = None
    dicom_files = os.listdir(dicom_path)
    for dicom_file in dicom_files:
        if not dicom_file.startswith('KO') and not dicom_file.startswith('.DS_Store'): 
            dicom_dataset = pydicom.filereader.dcmread(os.path.join(dicom_path, dicom_file))
            break
    return dicom_dataset

def output_metadata(dicom_path, output_path, verbose=False):
    """ Outputs all metadata from a given folder of DICOM files. Any metadata different between DICOM files is output as well.

    Parameters:
        dicom_path (str):The full path to the directory containing DICOM files
        output_path (str): The full path to the desired output text file, including the filename
    """
    start_time = None
    if verbose:
        start_time = time.time()
        print("Beginning output_metadata from DICOM in " + dicom_path + " to " + output_path + "..", end = '')
    
    dicom_folder = os.path.basename(os.path.normpath(output_path))
    
    metadata_file = open(output_path, "w+")
    metadata_file.write("\nFile -- Tag -- Keyword -- Name -- Value\n")
    metadata_file.write("\n***** Case " + dicom_folder + " *****\n")
    
    # Create a dictionary with key of data_element.keyword - output data elements that are different
    metadata_dict = {}
    
    dicom_files = os.listdir(dicom_path)                        
    for dicom_file in dicom_files:
        if dicom_file.endswith('.dcm') and not dicom_file.startswith('KO') and not dicom_file.startswith('.DS_Store'):
            dicom_dataset = pydicom.filereader.dcmread(os.path.join(dicom_path, dicom_file))
            for data_element in dicom_dataset:
                write_data_element = False
                if "pixel data" not in str(data_element.name).lower() and "overlay data" not in str(data_element.name).lower():
                    if data_element.keyword in metadata_dict:
                        # The data element has already been added once
                        if data_element.value not in metadata_dict[data_element.keyword]:
                            # The data element value is different though, so print again
                            metadata_dict[data_element.keyword].append(data_element.value)
                            write_data_element = True
                    else:
                        # The data element has not been added yet - create a list of 1 value and print
                        metadata_dict[data_element.keyword] = [data_element.value]
                        write_data_element = True
                if write_data_element and data_element.keyword not in ['InstanceNumber', 'ImageOrientationPatient', 'ImagePositionPatient', 'SliceLocation', 'WindowCenter', 'WindowWidth', 'SOPInstanceUID']:
                    metadata = str(dicom_file) + " | " + str(data_element.tag) + " | " + str(data_element.keyword) + " | " + str(data_element.name) + " | " + str(data_element.value) + "\n"
                    metadata_file.write(metadata)

    metadata_file.close()

    if verbose:
        print('.done - time = ' + str(np.around((time.time() - start_time), decimals=1)) + 'sec')    

def zip(dicom_path, output_path, compression=zipfile.ZIP_DEFLATED, compresslevel=9, series=None, verbose=False):
    """ Zips all *.dcm files within a given dicom_path. A different zip file is created for each series (SeriesDescription). Individual .dcm files are renamed based on their image number (InstanceNumber)

    Parameters:
        dicom_path (str):The full path to the directory containing DICOM files
        output_path (str): The full path to the desired output zip file (not including filename)
        compression: ZIP compression method to use when writing the archive, and should be ZIP_STORED, ZIP_DEFLATED, ZIP_BZIP2 or ZIP_LZMA; unrecognized values will cause NotImplementedError to be raised. If ZIP_DEFLATED, ZIP_BZIP2 or ZIP_LZMA is specified but the corresponding module (zlib, bz2 or lzma) is not available, RuntimeError is raised. The default is ZIP_STORED.
        compresslevel: Controls the compression level to use when writing files to the archive. When using ZIP_STORED or ZIP_LZMA it has no effect. When using ZIP_DEFLATED integers 0 through 9 are accepted (see zlib for more information). When using ZIP_BZIP2 integers 1 through 9 are accepted (see bz2 for more information).
        series (str): Optional series name. Only works if there is 1 series in the DICOM.
    """
    start_time = None
    if verbose:
        start_time = time.time()
        print("Beginning zip of DICOM in " + dicom_path + " to " + output_path + "..", end = '')

    zipfile_name = None
    # Iterate through .dcm files - add each file to an archive based on its series name
    
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.mkdir(output_path)
    
    # Move DICOM into zip files based on series
    series_dict = {}
    for dcm_file in os.listdir(dicom_path):
        if dcm_file.endswith('.dcm') and not dcm_file.startswith('KO') and not dcm_file.startswith('.DS_Store'):
            dicom_dataset = pydicom.filereader.dcmread(os.path.join(dicom_path, dcm_file)) 
            
            series = ""
            
            series_instance_uid = dicom_dataset.SeriesInstanceUID
            if series_instance_uid in series_dict:
                # The series name has already been assigned
                series = series_dict[series_instance_uid]
            else:
                # The series name has not been assigned
                # Make sure there aren't 2 series with the same name
                if "SeriesDescription" in dicom_dataset:
                    # Make the series name the SeriesDescription
                    series = str(dicom_dataset.SeriesDescription)
                if series == "":
                    if "SeriesNumber" in dicom_dataset:
                        series = str(dicom_dataset.SeriesNumber)
                if series == "":
                    series = "Series"

                # Eliminate special characters: / \ : * ? " < > |
                series = filter_series_name(series)
                    
                if series in series_dict.values():
                    # There are multiple series of the same name
                    # Generate a unique name for the series
                    new_series_name = ""
                    ended = False
                    i = 2
                    while not ended:
                        new_series_name = series + " " + str(i)
                        i = i + 1
                        if new_series_name not in series_dict.values():
                            ended = True
                    series_dict[series_instance_uid] = new_series_name
                    series = new_series_name
                else:
                    # This is a new series - add to dictionary
                    series_dict[series_instance_uid] = series                    

            if "dose report" not in series.lower() and "screen save" not in series.lower() and "summary" not in series.lower():                                
                zipfile_name = series + ".zip"
                
                with ZipFile(os.path.join(output_path, zipfile_name), 'a') as myzip:
                    myzip.write(os.path.join(dicom_path, dcm_file), arcname=dcm_file)
                    myzip.close()
    
    # Rename .dcm files for the image number (InstanceNumber) within each zip file
    # Unzip the .dcm files into the same directory
    for zip_file in os.listdir(output_path):
        if zip_file.endswith('.zip'):
            series = os.path.splitext(zip_file)[0]
            with ZipFile(os.path.join(output_path, zip_file), 'a') as myzip:
                myzip.extractall(output_path)
                myzip.close()
            os.remove(os.path.join(output_path, zip_file))      
            
            # key = filename, value = number with that filename - 1
            dcm_filename_dict = {}
            for dcm_file in os.listdir(output_path):
                if dcm_file.endswith('.dcm') and not dcm_file.startswith('KO') and not dcm_file.startswith('.DS_Store'):
                    dicom_dataset = pydicom.filereader.dcmread(os.path.join(output_path, dcm_file)) 

                    if "dose report" not in series.lower() and "screen save" not in series.lower() and "summary" not in series.lower():                                
                        zipfile_name = series + ".zip"
                        
                        # Name the .dcm file for the slice number (InstanceNumber)
                        new_dcm_file = str(dicom_dataset.InstanceNumber).zfill(7)
                        if new_dcm_file in dcm_filename_dict.keys():
                            new_dcm_file_num = int(dcm_filename_dict[new_dcm_file]) + 1
                            dcm_filename_dict[new_dcm_file] = new_dcm_file_num
                            new_dcm_file = new_dcm_file + "-" + str(new_dcm_file_num)
                        else:
                            dcm_filename_dict[new_dcm_file] = 1
                        
                        new_dcm_file = new_dcm_file + ".dcm"
                        
                        if (dcm_file != new_dcm_file):
                            os.rename(os.path.join(output_path, dcm_file), os.path.join(output_path, new_dcm_file))
                        
                        with ZipFile(os.path.join(output_path, zipfile_name), mode='a', compression=compression, compresslevel=compresslevel) as myzip:
                            myzip.write(os.path.join(output_path, new_dcm_file), arcname=new_dcm_file)
                            myzip.close()
            
            for dcm_file in os.listdir(output_path):
                if dcm_file.endswith('.dcm'):
                    os.remove(os.path.join(output_path, dcm_file))
    if verbose:
        print('.done - time = ' + str(np.around((time.time() - start_time), decimals=1)) + 'sec') 

    return zipfile_name