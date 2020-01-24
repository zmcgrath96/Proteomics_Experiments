from pyopenms import MSExperiment, MzMLFile
from utils import __file_exists

'''read

DESC:
    read an .mzML file into memory
PARAMS:
    file: str path to the file to import
RETURNS:
    NONE if file is not found

'''
def read(file):
    if not __file_exists(file):
        print('File {} not found. Please make sure that this file exists'.format(file))
        return
    else: 
        spectra = []
        exp = MSExperiment()
        MzMLFile().load(file, exp)

        for s in exp.getSpectra():
            spectra.append({
                'level': s.getMSLevel(),
                'spectrum': list(list(s.get_peaks())[0]), 
                'scan_no': int(str(s.getNativeID()).split('=')[-1].replace("'", ''))
            })

        return spectra