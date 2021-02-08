# dicom-utils
This is a Python library for manipulating DICOM datasets. It contains some of the most useful methods I've developed. These include:
* anonymize: Anonymize DICOM (tested on CT, MRI, and radiographs)
* dcm2jpg: Convert DICOM to jpg
* dcm2pydicom: Load DICOM as pydicom
* dcm2sitk: Load DICOM as SimpleITK
* output_metadata: Output all unique metadata into a text file to confirm anonymization
* zip_dicom: Compress DICOM into a zip file. Helpful if also using the simple DWV in a separate repository of mine. This method does the max compression possible to still have it readable in DWV.

Importantly, I'm sharing these in the sprit of collaboration, hoping they are helpful to someone. This library of course comes with no guarantees, and you assume all responsibility for ensuring they behave as intended. For example, DICOM anonymization is a thorny process, and I would highly recommend you check your metadata to ensure it is actually anonymized as expected.

I'm also relatively new to GitHub. I welcome any suggestions or feedback.

## Getting DICOM
First, you may be wondering: How do I even get the DICOM files to begin with? There are numerous different techniques - the best first step is to talk with your PACS support. This guide assumes you already have the DICOM files. For convenience, this code also comes with 3 small test CT abdomen studies downloaded from TCIA. These are located in tests/data/tcia_ct_abdomens.

## System Requirements
Python 3.7

## Getting Started
I would recommend you start by running the tests in tests/test_dicom_utils.py using the unittest framework. The 2 tests go through usage of all the methods in dicom_utils.py.

## DICOM Anonymization Method
A little more detail on the DICOM anonymization method. This is a surprisingly difficult task, as DICOM is not as standardized as one might think. DICOM studies are littered with custom tags that can cause issues. So, my anonymization method is very explicit: I only keep tags that I need in order to view the images, and there is a long list of them at the beginning of the file.

For research projects, you may want to make changes to the DICOM metadata. For example, maybe you generate a research patient ID that you want to become the MRN of the anonymized study. For that reason, there are several optional arguments to the anonymize method that allow for that.

To confirm everything anonymized completely, I suggest using the output_metadata method. This walks all the metadata in every DICOM file in the provided folder, and outputs all unique data elements to a text file. That way you don't get burnt by PHI contained in the 61st DICOM file that isn't in any other files. If the MRN is different in 1 single file, that will get added to the text file. You can then review the text file to confirm the absence of any PHI. Again, have a high level of paranoia with this and consider verifying independently.

## dcm2jpg
This is a cool little method that converts DICOM to JPG in bulk, useful for creating scrollable experiences in PowerPoint. If you look in the method code, you can see how window width and levelling actually work.

## Tips for the Beginner
While a full walkthrough is beyond the scope of this README, here is a short guide to my entire image processing setup:
* OS: Ubuntu 18.04
* Python: Install Anaconda, then use conda on the command line to create a Python 3.7 virtual environment. This is HIGHLY recommended instead of dealing with the hassle of using the system Python, you can easily just create and destroy Python verions.
* Development Environment: I use Visual Studio Code, which is amazing. Made by Microsoft, but available on every OS and very intuitive. Allows you to easily run and debug code. Also integrates nicely with conda and unittest.
* Getting DICOM Files: I run DCMTK (free) on the Ubuntu machine, and can send DICOM studies from PACS directly to my machine. Definitely the most efficient way to get DICOM studies.
* Storing Studies: DICOM studies can get HUGE. I have a Box unlimited storage account for $35 per month. There's no Box sync client for Ubuntu (arghhh), so I have a Mac Mini, and mount the Mac's Box directory on my Ubuntu machine. Kludgy, but it works.
