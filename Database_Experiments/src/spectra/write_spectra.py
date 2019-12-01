from pyopenms import MSExperiment, MSSpectrum, MzMLFile, Peak1D, Precursor

def write_mzml(output_dir, file_name, spectra, title_prefix='Spectrum '):
    if '.mgf' not in file_name:
        file_name += '.mzml'
    output_file = output_dir + '/' + file_name

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

    return output_file