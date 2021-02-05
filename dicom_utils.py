import cv2
import numpy as np
import os
import pydicom
import shutil
import SimpleITK as sitk
import sys
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
    """ Anonymizes DICOM files. Utilizes the above defined keep lists (modify as desired)

    Parameters:
        dicom_path (str):The full path to the directory containing DICOM files - these will not be modified
        output_path (str): The full path to the directory that will contain the output folder. Any existing files in this directory will be deleted
        preserve_list (list, optional): List of dataset keywords to not be altered in the anonymization process.
            Currently only applied to a small number of keywords (PatientAge, PatientBirthDate, PatientSex, StudyDate). Otherwise alter keep lists and potentially change code below.
        MRN (str, optional): New MRN
        Accession (str, optional): New accession number
        PatientName (str, optional): New patient name. Format = 'LASTNAME^FIRSTNAME'
        SeriesDescription (str, optional): New series description. 
        verbose (bool, optional): Verbose flag
    """
    start_time = None
    if verbose:
        start_time = time.time()
        print("Beginning anonymize of DICOM in " + dicom_path + " to " + output_path + "..", end = '')

    # Concatenate the keep lists
    keep_list = general_keep_list + \
        mr_keep_list + \
        angio_keep_list

    # Create the output directory
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.mkdir(output_path)

    # Copy original files into output_path
    dicom_files = os.listdir(dicom_path)
    for dicom_file in dicom_files:
        if os.path.isfile(os.path.join(dicom_path, dicom_file)) and not dicom_file.startswith('KO') and not dicom_file.startswith('.DS_Store'):
            shutil.copy(
                os.path.join(dicom_path, dicom_file),
                os.path.join(output_path, dicom_file)
                )

    # Add '.dcm' extension if necessary
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

            # Iterate through each metadata element of each DICOM file
            for data_element in ds:
                # Don't delete any actual image data
                if "pixel data" not in str(data_element.name).lower() \
                        and "overlay data" \
                        not in str(data_element.name).lower():
                    delete_de = True

                    # Don't delete the metadata element if it is in a keep list
                    if str(data_element.keyword) in keep_list:
                        delete_de = False

                    # Anonymize metadata elements that are necessary to keep
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

def dcm2jpg(dicom_path, output_path, window_width=None, window_level=None, verbose=False):
    """ Converts DICOM to 8 bit jpg image files.

    Parameters:
        dicom_path (str): Full path to folder containing DICOM files - will not be modified
        output_path (str): Full path to write image files to
        window_width (int, optional): If not specified, the default window width in the DICOM metadata will be used.
        window_level (int, optional): If not specified, the default window level / center in the DICOM metadata will be used.
        verbose (bool, optional): Verbose flag
    """
    start_time = None
    if verbose:
        start_time = time.time()
        print("Beginning dcm2jpg of " + dicom_path + " to " + output_path + "..", end = '')

    image_type = 'jpg'
  
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.mkdir(output_path)
    
    shutil.copytree(dicom_path, os.path.join(output_path, 'dicom_files'))
    dicom_path = os.path.join(output_path, 'dicom_files')
    for dicom_file in os.listdir(dicom_path):
        if not dicom_file.endswith('.dcm'):
            shutil.move(os.path.join(dicom_path, dicom_file), os.path.join(dicom_path, dicom_file + '.dcm'))

    # Load DICOM as sitk for automatic scaling of voxel values (e.g. to Hounsfield units for CT)
    image_reader, image = dcm2sitk(dicom_path)
    image_array = sitk.GetArrayFromImage(image)

    # Load DICOM as pydicom for access to window and level for each slice
    loaded_dicom = dcm2pydicom(dicom_path)

    for i in range(image_array.shape[0]):
        # Load the pydicom dataset for slice i
        ds = loaded_dicom[i]

        # Load the sitk pixel data for slice i
        slice_array = image_array[i,:,:]

        # Determine the window and level for the image - provided or default
        write_window_width = window_width
        write_window_level = window_level
        if window_width is None:
            if type(ds.WindowWidth) == pydicom.multival.MultiValue:
                write_window_width = str(ds.WindowWidth[0])
            else:
                write_window_width = str(ds.WindowWidth)
        if window_level is None:
            if type(ds.WindowCenter) == pydicom.multival.MultiValue:
                write_window_level = str(ds.WindowCenter[0])
            else:        
                write_window_level = str(ds.WindowCenter)

        # Set the window and level - anything above or below the window is white/black
        # This is pretty slick - if you understand the code below, you understand what window width and levels are
        scaled_slice_array = np.array(slice_array)
        min_value = int(write_window_level) - int(int(write_window_width)/2)
        max_value = int(write_window_level) + int(int(write_window_width)/2)
        min_indices = np.where(scaled_slice_array < min_value)
        scaled_slice_array[min_indices] = min_value
        max_indices = np.where(scaled_slice_array > max_value)
        scaled_slice_array[max_indices] = max_value

        # Scale the pixel values between 0 and 255 for an 8 bit image
        scaled_slice_array = np.interp(scaled_slice_array, (min_value, max_value), (0, 255)) 
        
        # Write the pixel data to an image
        image_filename = str(i).zfill(4) + "." + image_type
        cv2.imwrite(os.path.join(output_path, image_filename), scaled_slice_array)

    shutil.rmtree(os.path.join(output_path, 'dicom_files'))

    if verbose:
        print('.done - time = ' + str(np.around((time.time() - start_time), decimals=1)) + 'sec')

def dcm2pydicom(dicom_path, verbose=False):
    """ Loads DICOM files into pydicom format.

    Parameters:
        dicom_path (str): Full path to the directory containing DICOM files
        verbose (bool, optional): Verbose flag
    Returns:
        loaded_dicom (list): all DICOM files in pydicom format
    """
    start_time = None
    if verbose:
        start_time = time.time()
        print("Beginning dcm2pydicom from " + dicom_path + "..", end = '')

    # Initial load of DICOm into pydicom via dcmread
    loaded_dicom = []
    for dicom_file in os.listdir(dicom_path):
        ds = pydicom.dcmread(os.path.join(dicom_path, dicom_file))
        loaded_dicom.append(ds)
    
    # Sort the loaded DICOM by InstanceNumber (this is assumed to be the slice number)
    try:
        loaded_dicom.sort(key = lambda x: int(x.InstanceNumber))
    except:
        # Perfusion sequences don't have instance numbers
        print("Error in load.dicom_as_pydicom: An InstanceNumber was '' - not sorting")
    
    # Calculate and set the true slice thickness. This is often reported incorrectly in DICOM metadata.
    slice_thickness = None
    if hasattr(loaded_dicom[0], 'ImagePositionPatient') and len(loaded_dicom) > 1:
        slice_thickness = np.abs(loaded_dicom[0].ImagePositionPatient[2] - loaded_dicom[1].ImagePositionPatient[2])
    elif hasattr(loaded_dicom[0], 'SliceLocation') and len(loaded_dicom) > 1:
        slice_thickness = np.abs(loaded_dicom[0].SliceLocation - loaded_dicom[1].SliceLocation)
    if slice_thickness is not None:
        for s in loaded_dicom:
            s.SliceThickness = slice_thickness
    
    if verbose:
        print('.done - time = ' + str(np.around((time.time() - start_time), decimals=1)) + 'sec')

    return loaded_dicom    

def dcm2sitk(dicom_path, verbose=False):
    """ Loads DICOM files into SimpleITK format. Loads all tags including private tags.

    Parameters:
        dicom_path: Path to directory containing DICOM files. If multiple series are present, the first series detected will be loaded. (Generally recommended to just have one series in dicom_path)
        verbose (bool): Verbose flag
    Returns:
        sitk_reader, sitk_dicom: DICOM study read into SimpleITK format, reader inclued for writing back to DICOM
    """	
    start_time = None
    if verbose:
        start_time = time.time()
        print("Beginning dcm2sitk from " + dicom_path + "..", end = '')

    series_IDs = sitk.ImageSeriesReader.GetGDCMSeriesIDs(dicom_path)
    if not series_IDs:
        print("ERROR: given directory \""+dicom_path+"\" does not contain a DICOM series.")
        sys.exit(1)
    series_file_names = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(dicom_path, series_IDs[0])

    series_reader = sitk.ImageSeriesReader()
    series_reader.SetFileNames(series_file_names)

    # Configure the reader to load all of the DICOM tags (public+private):
    # By default tags are not loaded (saves time).
    # By default if tags are loaded, the private tags are not loaded.
    # We explicitly configure the reader to load tags, including the private ones.
    series_reader.MetaDataDictionaryArrayUpdateOn()
    series_reader.LoadPrivateTagsOn()
    sitk_dicom = series_reader.Execute()

    if verbose:
        print('.done - time = ' + str(np.around((time.time() - start_time), decimals=1)) + 'sec')

    return series_reader, sitk_dicom

def filter_series_name(series):
    """ Removes special characters from a series name. These are characters that specifically cause errors with the DICOM Web Viewer (DWV)

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
    """ Helper function to return a pydicom dataset given a path to a folder containing dicom files (just returns a dataset for the first DICOM file). 

    Parameters:
        dicom_path (str): The full path to the directory containing DICOM files
    Returns:
        dicom_dataset (pydicom.dataset.Dataset): A pydicom dataset for the first DICOM file. None if no DICOM files are present.
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
        output_path (str): The full path to the output text file including the filename
        verbose (bool, optional): Verbose flag
    """
    start_time = None
    if verbose:
        start_time = time.time()
        print("Beginning output_metadata from DICOM in " + dicom_path + " to " + output_path + "..", end = '')
    
    dicom_folder = os.path.basename(os.path.normpath(output_path))
    
    # Write the header
    metadata_file = open(output_path, "w+")
    metadata_file.write("\nFile -- Tag -- Keyword -- Name -- Value\n")
    metadata_file.write("\n***** Case " + dicom_folder + " *****\n")
    
    # Create a dictionary with key of data_element.keyword
    # This will track which metadata elements have been output. If there is a new element, or a new value for an existing element, either will be written.
    metadata_dict = {}
    
    dicom_files = os.listdir(dicom_path)                        
    for dicom_file in dicom_files:
        if dicom_file.endswith('.dcm') and not dicom_file.startswith('KO') and not dicom_file.startswith('.DS_Store'):
            dicom_dataset = pydicom.filereader.dcmread(os.path.join(dicom_path, dicom_file))
            for data_element in dicom_dataset:
                write_data_element = False
                if "pixel data" not in str(data_element.name).lower() and "overlay data" not in str(data_element.name).lower(): # Don't write out the actual image data
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
                # The following keywords change with each slice, and are virtually guaranteed to not contain PHI, so we don't write those out repeatedlyl
                if write_data_element and data_element.keyword not in ['InstanceNumber', 'ImageOrientationPatient', 'ImagePositionPatient', 'SliceLocation', 'WindowCenter', 'WindowWidth', 'SOPInstanceUID']:
                    metadata = str(dicom_file) + " | " + str(data_element.tag) + " | " + str(data_element.keyword) + " | " + str(data_element.name) + " | " + str(data_element.value) + "\n"
                    metadata_file.write(metadata)

    metadata_file.close()

    if verbose:
        print('.done - time = ' + str(np.around((time.time() - start_time), decimals=1)) + 'sec')    

def zip_dicom(dicom_path, output_path, compression=zipfile.ZIP_DEFLATED, compresslevel=9, series=None, verbose=False):
    """ Zips all *.dcm files within a given dicom_path. A different zip file is created for each series (SeriesDescription). 
        Individual .dcm files are renamed based on their image number (InstanceNumber)
        Default settings will compress at the maximum amount still readable by the DICOM Web Viewer (DWV)

    Parameters:
        dicom_path (str):The full path to the directory containing DICOM files
        output_path (str): The full path to the desired output zip file (not including filename)
        compression (constant, optional): ZIP compression method to use when writing the archive, and should be ZIP_STORED, ZIP_DEFLATED, ZIP_BZIP2 or ZIP_LZMA; unrecognized values will cause NotImplementedError to be raised. If ZIP_DEFLATED, ZIP_BZIP2 or ZIP_LZMA is specified but the corresponding module (zlib, bz2 or lzma) is not available, RuntimeError is raised. The default is ZIP_STORED.
        compresslevel (int, optional): Controls the compression level to use when writing files to the archive. When using ZIP_STORED or ZIP_LZMA it has no effect. When using ZIP_DEFLATED integers 0 through 9 are accepted (see zlib for more information). When using ZIP_BZIP2 integers 1 through 9 are accepted (see bz2 for more information).
        series (str, optional): Manually specify the series name. Only works if there is 1 series in the DICOM.
        verbose (bool, optional): Verbose flag
    """
    start_time = None
    if verbose:
        start_time = time.time()
        print("Beginning zip of DICOM in " + dicom_path + " to " + output_path + "..", end = '')

    zipfile_name = None
    
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.mkdir(output_path)
    
    # Move DICOM into folders based on series
    # Sometimes the same Accession will have multiple series of the same name. SeriesInstanceUID is the true unique ID for a series.
    # If there are multiple series of the same name, add a number to subsequent series for the name. E.g. 'Ax CT' and 'Ax CT 2'
    # series_dict will track the SeriesInstanceUID and series name assigned to that UID
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

            # Do not include series likely to contain burned-in PHI.
            # The if statement below is EXTREMELY site-specific and unreliable. You must manually confirm there is no burned-in PHI.
            if "dose report" not in series.lower() and "screen save" not in series.lower() and "summary" not in series.lower():                                
                '''
                zipfile_name = series + ".zip"
                
                with ZipFile(os.path.join(output_path, zipfile_name), 'a') as myzip:
                    myzip.write(os.path.join(dicom_path, dcm_file), arcname=dcm_file)
                    myzip.close()
                '''
                series_path = os.path.join(output_path, series)
                if not os.path.exists(series_path):
                    os.mkdir(series_path)
                shutil.copy(os.path.join(dicom_path, dcm_file), os.path.join(series_path, dcm_file))
    
    # Iterate through folders and rename .dcm files based on image number (InstanceNumber)
    # Then zip all files in a folder
    for series in os.listdir(output_path):
        series_path = os.path.join(output_path, series)
        if os.path.isdir(series_path):
            # key = filename, value = number of images with that filename - 1
            dcm_filename_dict = {}
            for dcm_file in os.listdir(series_path):
                if dcm_file.endswith('.dcm') and not dcm_file.startswith('KO') and not dcm_file.startswith('.DS_Store'):
                    dicom_dataset = pydicom.filereader.dcmread(os.path.join(series_path, dcm_file)) 

                    if "dose report" not in series.lower() and "screen save" not in series.lower() and "summary" not in series.lower():                                
                        zipfile_name = series + ".zip"
                        
                        # Name the .dcm file for the slice number (InstanceNumber)
                        new_dcm_file = str(dicom_dataset.InstanceNumber).zfill(7)
                        # Occasionally modalities have multiple images with the same InstanceNumber
                        if new_dcm_file in dcm_filename_dict.keys():
                            # There is already an image with that InstanceNumber
                            # Create a new filename. E.g. instead of 0000001.dcm, name is 0000001-2.dcm
                            new_dcm_file_num = int(dcm_filename_dict[new_dcm_file]) + 1 # Number of files with the same InstanceNumber 
                            dcm_filename_dict[new_dcm_file] = new_dcm_file_num # Update the number in the dict
                            new_dcm_file = new_dcm_file + "-" + str(new_dcm_file_num) # Create the new 'dashed' filename
                        else:
                            dcm_filename_dict[new_dcm_file] = 1 # The majority - no duplicate InstanceNumbers, there is just '1' image with the InstanceNumber
                        
                        new_dcm_file = new_dcm_file + ".dcm"
                        
                        if (dcm_file != new_dcm_file):
                            # Rename the .dcm file with the new filename (the 'if' above is necessary because os.rename will fail if they're the same)
                            os.rename(os.path.join(series_path, dcm_file), os.path.join(series_path, new_dcm_file))
                        
                        # Write out the compressed zip
                        with ZipFile(os.path.join(output_path, zipfile_name), mode='a', compression=compression, compresslevel=compresslevel) as myzip:
                            myzip.write(os.path.join(series_path, new_dcm_file), arcname=new_dcm_file)
                            myzip.close()
            
    # Clean up the dicom folders
    for series in os.listdir(output_path):
        series_path = os.path.join(output_path, series)
        if os.path.isdir(series_path):
            shutil.rmtree(series_path)

    if verbose:
        print('.done - time = ' + str(np.around((time.time() - start_time), decimals=1)) + 'sec')