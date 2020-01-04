from pyopenms import MSExperiment, MSSpectrum, MzMLFile, Peak1D, Precursor
from utils import __make_valid_dir_string, __make_dir, __gzip

'''write_mzml

DESC:
    create a mass spectrum file in mzml of sequences
PARAMS:
    file_name: str name to save the file in
    spectra: list of dictionaries of the form [{spectrum: list[floats], precursor_mass: float, sequence: str}]
             This data is written to file (the spectrum and the precursor)
OPTIONAL:
    title_prefix: str name to give as prefix to the name of each spectrum. Default=Spectrum
    output_dir: str name of the directory to save files to. Default=./
    compress: bool whether or not to compress the file. Comrpesses with gzip. Default=True
RETURNS:
    list of strings of file paths
'''
def write_mzml(file_name, spectra, title_prefix='Spectrum ', output_dir='./', compress=True):
    if '.mzml' not in file_name.lower():
        file_name += '.mzML'
    output_dir = __make_valid_dir_string(output_dir)
    __make_dir(output_dir)
    output_file = output_dir + file_name

    exp = MSExperiment()
    sp_count = 0
    for spectrum in spectra:
        spec = MSSpectrum()
        spec.setMSLevel(2)
        name = str.encode(title_prefix + str(sp_count))
        spec.setName(name)
        sp_count += 1
        
        i = [500 for _ in spectrum['spectrum']]
        spec.set_peaks([spectrum['spectrum'], i])
        spec.setMSLevel(2)
        prec = Precursor()
        prec.setCharge(2)
        prec.setMZ(spectrum['precursor_mass'])
        spec.setPrecursors([prec])
        spec.sortByPosition()
        exp.addSpectrum(spec)

    MzMLFile().store(output_file, exp)
    output_file = __gzip(output_file) if compress else output_file

    return output_file