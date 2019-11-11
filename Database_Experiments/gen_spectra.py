from pyopenms import TheoreticalSpectrumGenerator, MSSpectrum, AASequence

def gen_spectra(sequenences):
    spectra = []
    tsg = TheoreticalSpectrumGenerator()
    for sequence in sequenences:
        spec = MSSpectrum()
        pep = AASequence.fromString(sequence)
        tsg.getSpectrum(spec, pep, 1, 1)

        this_spectrum = [peak.getMZ() for peak in spec]
        this_spectrum.sort()
        spectra.append(this_spectrum)

    return spectra
