# --- Begin $MAIN$ program. ---------------
#

# Use Atomic Data Collection from Porter et al (2012) and (2013)
# 2012MNRAS.425L..28P and 2013MNRAS.433L..89P
import atomneb
import numpy as np
import os

# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('..','atomic-data-rc')
atom_rc_file = os.path.join(base_dir,data_dir, 'rc_he_ii_PFSD12.fits')

atom = 'he'
ion = 'ii' # He I
# read Recombination Coefficients (Aeff) of He I
hei_rc_data = atomneb.read_aeff_he_i_pfsd12(atom_rc_file, atom, ion)
hei_rc_data_wave = atomneb.read_aeff_he_i_pfsd12(atom_rc_file, atom, ion, wavelength=True)
# print information needed for Recombination Coefficients (Aeff) of He I
print(hei_rc_data.aeff[0])
n_line = len(hei_rc_data_wave.wavelength)
for i in range(0, n_line):
   print(hei_rc_data_wave.wavelength[i], hei_rc_data_wave.lowerterm[i], hei_rc_data_wave.upperterm[i])

atom = 'he'
ion = 'ii' # He I
# list all Recombination Coefficients (Aeff) data for He I
list_hei_aeff_data = atomneb.search_aeff_he_i_pfsd12(atom_rc_file, atom, ion)
# print all Recombination Coefficients (Aeff) of He I
print(list_hei_aeff_data)

atom = 'he'
ion = 'ii' # He I
# list all Recombination Coefficients (Aeff) references for He I
list_hei_aeff_reference = atomneb.list_aeff_he_i_pfsd12_references(atom_rc_file, atom, ion)
# print all Recombination Coefficients (Aeff) References for He I
print(list_hei_aeff_reference)

atom = 'he'
ion = 'ii' # He I
reference = 'PFSD13'
# get citations for Recombination Coefficients (Aeff) of He I
citation = atomneb.get_aeff_he_i_pfsd12_reference_citation(atom_rc_file, atom, ion, reference=reference)
# print citations for Recombination Coefficients (Aeff) of He I
print(citation)

