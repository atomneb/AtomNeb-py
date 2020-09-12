# --- Begin $MAIN$ program. ---------------
#

# Use Atomic Data Collection from Fang, Storey and Liu (2011) and (2013)
# 2011A&A...530A..18F and 2013A&A...550C...2F
import atomneb
import numpy as np
import os

# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('..','atomic-data-rc')
atom_rc_file = os.path.join(base_dir,data_dir, 'rc_n_iii_FSL13.fits')

atom = 'n'
ion = 'iii' # N II
# read Recombination Coefficients (Aeff) of N II
wavelength_range = [4400.0, 7100.0]
nii_rc_data = atomneb.read_aeff_n_ii_fsl13(atom_rc_file, atom, ion, wavelength_range) # it takes a while!
nii_rc_data_wave = atomneb.read_aeff_n_ii_fsl13(atom_rc_file, atom, ion, wavelength_range, wavelength=True)
# print information needed for Recombination Coefficients (Aeff) of N II
print(nii_rc_data.aeff[0])
n_line = len(nii_rc_data_wave.wavelength)
for i in range(0, n_line):
   print(nii_rc_data_wave.wavelength[i], nii_rc_data_wave.tr[i], nii_rc_data_wave.trans[i])

atom = 'n'
ion = 'iii' # N II
# list all Recombination Coefficients (Aeff) data for N II
wavelength = 5679.56
list_nii_aeff_data = atomneb.search_aeff_n_ii_fsl13(atom_rc_file, atom, ion, wavelength)
# print all Recombination Coefficients (Aeff) of N II
print(np.asarray(list_nii_aeff_data.wavelength))
print(np.asarray(list_nii_aeff_data.aeff))

atom = 'n'
ion = 'iii' # N II
# list all Recombination Coefficients (Aeff) references for N II
list_nii_aeff_reference = atomneb.list_aeff_n_ii_fsl13_references(atom_rc_file, atom, ion)
# print all Recombination Coefficients (Aeff) References for N II
print(list_nii_aeff_reference)

atom = 'n'
ion = 'iii' # N II
reference = 'FSL13'
# get citations for Recombination Coefficients (Aeff) of N II
citation = atomneb.get_aeff_n_ii_fsl13_reference_citation(atom_rc_file, atom, ion)
# print citations for Recombination Coefficients (Aeff) of N II
print(citation)

