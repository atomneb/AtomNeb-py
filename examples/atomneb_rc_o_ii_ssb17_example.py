# --- Begin $MAIN$ program. ---------------
#

# Use Atomic Data Collection from Storey, Sochi and Bastin (2017)
# 2017MNRAS.470..379S
import atomneb
import numpy as np
import os
import tarfile

# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('..','atomic-data-rc')
atom_rc_file = os.path.join(base_dir,data_dir, 'rc_o_iii_SSB17_orl_case_b.fits')

# unpack rc_o_iii_SSB17_orl_case_b.tar.gz
atom_rc_file_tar_gz = os.path.join(base_dir,data_dir, 'rc_o_iii_SSB17_orl_case_b.fits.tar.gz')
atom_rc_path = os.path.join(base_dir,data_dir)
tar = tarfile.open(atom_rc_file_tar_gz, "r:gz")
tar.extractall(path=atom_rc_path)
tar.close()

atom = 'o'
ion = 'iii' # O II
case1 = 'B'
# read Recombination Coefficients (Aeff) of O II
wavelength_range = [5320.0, 5330.0]
oii_rc_data = atomneb.read_aeff_o_ii_ssb17(atom_rc_file, atom, ion, case1, wavelength_range)
oii_rc_data_wave = atomneb.read_aeff_o_ii_ssb17(atom_rc_file, atom, ion, case1, wavelength_range, wavelength=True)
# print information needed for Recombination Coefficients (Aeff) of O II
print(oii_rc_data.aeff[0])
n_line = len(oii_rc_data_wave.wavelength)
for i in range(0, n_line):
   print(oii_rc_data_wave.wavelength[i], oii_rc_data_wave.lower_term[i], oii_rc_data_wave.upper_term[i])

atom = 'o'
ion = 'iii' # O II
case1 = 'B'
# list all Recombination Coefficients (Aeff) data for O II
wavelength = 5325.42
list_oii_aeff_data = atomneb.search_aeff_o_ii_ssb17(atom_rc_file, atom, ion, case1, wavelength)
# print all Recombination Coefficients (Aeff) of O II
print(np.asarray(list_oii_aeff_data.wavelength))
print(np.asarray(list_oii_aeff_data.aeff))

atom = 'o'
ion = 'iii' # O II
# list all Recombination Coefficients (Aeff) references for O II
list_oii_aeff_reference = atomneb.list_aeff_o_ii_ssb17_references(atom_rc_file, atom, ion)
# print all Recombination Coefficients (Aeff) References for O II
print(list_oii_aeff_reference)

atom = 'o'
ion = 'iii' # O II
reference = 'SSB17'
# get citations for Recombination Coefficients (Aeff) of O II
citation = atomneb.get_aeff_o_ii_ssb17_reference_citation(atom_rc_file, atom, ion)
# print citations for Recombination Coefficients (Aeff) of O II
print(citation)

