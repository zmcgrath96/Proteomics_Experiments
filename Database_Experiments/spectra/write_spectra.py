from pyopenms import MSExperiment, MSSpectrum, MzMLFile

def write_mzml(output_dir, file_name, spectra, title_prefix='Spectrum '):
    if '.mgf' not in file_name:
        file_name += '.mzml'
    output_file = output_dir + '/' + file_name

    exp = MSExperiment()
    all_spec = []
    for spectrum in spectra:
        spec = MSSpectrum()
        i = [0 for _ in spectrum]
        spec.set_peaks([spectrum, i])
        spec.setMSLevel(2)
        spec.sortByPosition()
        all_spec.append(spec)
    exp.setSpectra(all_spec)
    MzMLFile().store(output_file, exp)

    return output_file