from pyteomics import mgf

def write_mgf(output_dir, file_name, spectra, title_prefix='Spectrum '):
    if '.mgf' not in file_name:
        file_name += '.mgf'
    output_file = output_dir + '/' + file_name

    sp_count = 0
    pyteomics_writable = []
    for spectrum in spectra:
        this_spectrum = {'m/z array': [], 'intensity array': [], 'params': {}}
        this_spectrum['m/z array'].append(spectrum)
        this_spectrum['params']['title'] = title_prefix + str(sp_count)
        pyteomics_writable.append(this_spectrum)
        sp_count += 1
    mgf.write(pyteomics_writable, output=output_file)
    return output_file