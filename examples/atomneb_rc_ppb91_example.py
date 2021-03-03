# --- Begin $MAIN$ program. ---------------
#

# Use Atomic Data Collection from the National Institute of Standards and Technology (NIST)
# Atomic Spectra Database, the CHIANTI atomic database, and some improved atomic data from
# Cloudy v13.04 and pyNeb v1.0
import atomneb
import numpy as np
import os

# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('..','atomic-data-rc')
atom_rc_file = os.path.join(base_dir,data_dir, 'rc_PPB91.fits')

atom = 'c'
ion = 'iii' # C II
# read Recombination Coefficients (Aeff) of C II
cii_rc_data = atomneb.read_aeff_ppb91(atom_rc_file, atom, ion)
n_line = len(cii_rc_data['wavelength'])
# print information needed for Recombination Coefficients (Aeff) of C II
for i in range(0, n_line):
   print(cii_rc_data['ion'][i], cii_rc_data['case1'][i], cii_rc_data['wavelength'][i],
         cii_rc_data['a'][i], cii_rc_data['b'][i], cii_rc_data['c'][i],
         cii_rc_data['d'][i], cii_rc_data['br'][i], cii_rc_data['q'][i], cii_rc_data['y'][i])

atom = 'c'
ion = 'iii'
# list all Recombination Coefficients (Aeff) data for C II
list_cii_aeff_data = atomneb.search_aeff_ppb91(atom_rc_file, atom, ion)
# print all Recombination Coefficients (Aeff) of C II
print(list_cii_aeff_data)

atom = 'c'
ion = 'iii'
# list all Recombination Coefficients (Aeff) references for C II
list_cii_aeff_reference = atomneb.list_aeff_ppb91_references(atom_rc_file, atom, ion)
# print all Recombination Coefficients (Aeff) References for C II
print(list_cii_aeff_reference)

atom = 'c'
ion = 'iii'
# get citations for Recombination Coefficients (Aeff) of C II
citation = atomneb.get_aeff_ppb91_reference_citation(atom_rc_file, atom, ion)
# print citations for Recombination Coefficients (Aeff) of C II
print(citation)

