from pyteomics import mgf

def load_mgf(mgf_file_path):
    spectra = []
    mgf_reader = mgf.read(mgf_file_path)
    for spectrum in mgf_reader:
        spectra.append(spectrum)
    print('wait')