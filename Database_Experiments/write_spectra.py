def write_mgf(output_dir, file_name, spectra, title_prefix='Spectrum '):
    if '.mgf' not in file_name:
        file_name += '.mgf'
    output_file = output_dir + '/' + file_name

    print('Writing...')
    total_spectra = len(spectra)
    sp_count = 0
    with open(output_file, 'w') as o:
        for spectrum in spectra:
            print('Writing spectra {}/{} [{}%] to {}\r'.format(sp_count, total_spectra, int((sp_count/total_spectra)*100), output_file), end="")
            o.write('BEGIN IONS\n')
            o.write('TITLE={}\n'.format(title_prefix + str(sp_count)))
            for mass in spectrum:
                o.write(str(mass) + '\n')
            o.write('END IONS\n')
            sp_count += 1
    print('Done')