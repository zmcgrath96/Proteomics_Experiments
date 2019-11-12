from spectra import load_spectrum

def load_spectra(iter_of_spectra_files):
    spectra = []
    for spectra_file in iter_of_spectra_files:
        spectra.append(load_spectrum.load_mgf(spectra_file))

    return spectra
    