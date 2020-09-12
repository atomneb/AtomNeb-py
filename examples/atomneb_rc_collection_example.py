# --- Begin $MAIN$ program. ---------------
#

# Use Atomic Data Collection for C II (Davey et al. 2000 2000A&AS..142...85D),
# N II (Escalante and Victor 1990 1990ApJS...73..513E),
# O II (Storey 1994 1994A&A...282..999S; Liu et al. 1995 1995MNRAS.272..369L),
# and Ne II ions (Kisielius et al. 1998 1998A&AS..133..257K)
import atomneb
import numpy as np
import os

# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('..','atomic-data-rc')
atom_rc_file = os.path.join(base_dir,data_dir, 'rc_collection.fits')

atom = 'c'
ion = 'iii' # C II
# read Recombination Coefficients (Aeff) of C II
cii_rc_data = atomneb.read_aeff_collection(atom_rc_file, atom, ion)
n_line = len(cii_rc_data.wavelength)
# print information needed for Recombination Coefficients (Aeff) of C II
for i in range(0, n_line):
  print(cii_rc_data.wavelength[i], cii_rc_data.a[i], cii_rc_data.b[i], cii_rc_data.c[i], cii_rc_data.d[i], cii_rc_data.f[i])

atom = 'n'
ion = 'iii' # N II
# read Recombination Coefficients (Aeff) of N II
nii_rc_data = atomneb.read_aeff_collection(atom_rc_file, atom, ion)
nii_rc_data_br = atomneb.read_aeff_collection(atom_rc_file, atom, ion, br=True)
n_line = len(nii_rc_data.a)
# print information needed for Recombination Coefficients (Aeff) of N II
for i in range(0, n_line):
  print(nii_rc_data.a[i], nii_rc_data.b[i], nii_rc_data.c[i])
n_line = len(nii_rc_data_br.wavelength)
# print forBranching ratios (Br) of N II
for i in range(0, n_line):
  print(nii_rc_data_br.wavelength[i], nii_rc_data_br.br[i],
        nii_rc_data_br.g1[i], nii_rc_data_br.g2[i],
        nii_rc_data_br.mult1[i], nii_rc_data_br.lowerterm[i], nii_rc_data_br.upperterm[i])

atom = 'o'
ion = 'iii' # O II
# read Recombination Coefficients (Aeff) of O II
oii_rc_data = atomneb.read_aeff_collection(atom_rc_file, atom, ion)
oii_rc_data_br = atomneb.read_aeff_collection(atom_rc_file, atom, ion, br=True)
n_line = len(oii_rc_data.a2)
# print information needed for Recombination Coefficients (Aeff) of O II
for i in range(0, n_line):
  print(oii_rc_data.term[i], oii_rc_data.case1[i],
        oii_rc_data.a2[i], oii_rc_data.a4[i], oii_rc_data.a5[i], oii_rc_data.a6[i],
        oii_rc_data.b[i], oii_rc_data.c[i], oii_rc_data.d[i])
n_line = len(oii_rc_data_br.wavelength)
# print forBranching ratios (Br) of O II
for i in range(0, n_line):
  print(oii_rc_data_br.wavelength[i], oii_rc_data_br.br_a[i], oii_rc_data_br.br_b[i], oii_rc_data_br.br_c[i],
        oii_rc_data_br.g1[i], oii_rc_data_br.g2[i], oii_rc_data_br.mult1[i],
        oii_rc_data_br.lowerterm[i], oii_rc_data_br.upperterm[i])

atom = 'ne'
ion = 'iii' # Ne II
# read Recombination Coefficients (Aeff) of Ne II
neii_rc_data = atomneb.read_aeff_collection(atom_rc_file, atom, ion)
n_line = len(neii_rc_data.wavelength)
# print information needed for Recombination Coefficients (Aeff) of Ne II
for i in range(0, n_line):
  print(neii_rc_data.wavelength[i],
        neii_rc_data.a[i], neii_rc_data.b[i], neii_rc_data.c[i], neii_rc_data.d[i], neii_rc_data.f[i], neii_rc_data.br[i])

atom = 'c'
ion = 'iii' # C III
# list all Recombination Coefficients (Aeff) data for C III
list_cii_aeff_data = atomneb.search_aeff_collection(atom_rc_file, atom, ion)
# print all Recombination Coefficients (Aeff) of C III
print(list_cii_aeff_data)

atom = 'c'
ion = 'iii' # C III
# list all Recombination Coefficients (Aeff) references for C III
list_cii_aeff_reference = atomneb.list_aeff_collection_references(atom_rc_file, atom, ion)
# print all Recombination Coefficients (Aeff) References for C III
print(list_cii_aeff_reference)

atom = 'c'
ion = 'iii' # C III
# get citations for Recombination Coefficients (Aeff) of C III with reference SSB14
citation = atomneb.get_aeff_collection_reference_citation(atom_rc_file, atom, ion)
# print citations for Recombination Coefficients (Aeff) of C III with reference SSB14
print(citation)


