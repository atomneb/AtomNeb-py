# --- Begin $MAIN$ program. ---------------
#

# Use Atomic Data Collection from Storey and Hummer (1995) 1995MNRAS.272...41S
import atomneb
import numpy as np
import os

# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('..','atomic-data-rc')
atom_rc_file = os.path.join(base_dir,data_dir, 'rc_SH95.fits')

atom = 'h'
ion = 'ii' # H I
# read Recombination Coefficients (Aeff) of H I
hi_rc_data = atomneb.read_aeff_sh95(atom_rc_file, atom, ion)
# print information needed for Recombination Coefficients (Aeff) of He I
print(hi_rc_data.aeff[0])

atom = 'h'
ion = 'ii' # H I
# list all Recombination Coefficients (Aeff) data for H I
list_hi_aeff_data = atomneb.search_aeff_sh95(atom_rc_file, atom, ion)
# print all Recombination Coefficients (Aeff) of H I
print(list_hi_aeff_data)


atom = 'h'
ion = 'ii' # H I
# list all Recombination Coefficients (Aeff) references for H I
list_hi_aeff_reference = atomneb.list_aeff_sh95_references(atom_rc_file, atom, ion)
# print all Recombination Coefficients (Aeff) References for H I
print(list_hi_aeff_reference)

atom = 'h'
ion = 'ii' # H I
# get citations for Recombination Coefficients (Aeff) of H I
citation = atomneb.get_aeff_sh95_reference_citation(atom_rc_file, atom, ion)
# print citations for Recombination Coefficients (Aeff) of H I
print(citation)

