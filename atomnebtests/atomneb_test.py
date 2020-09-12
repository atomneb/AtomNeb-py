"""Tests for atomneb"""

import atomneb
import numpy as np
import os

# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('atomic-data', 'chianti52')
atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')


# read Energy Levels (Ej) list
elj_data_list = atomneb.read_elj_list(atom_elj_file)
# read Collision Strengths (Omegaij) list
omij_data_list = atomneb.read_omij_list(atom_omij_file)
# read Transition Probabilities (Aij) list
aij_data_list = atomneb.read_aij_list(atom_aij_file)

# read Energy Levels (Ej) references
elj_data_reference = atomneb.read_elj_references(atom_elj_file)
# read Collision Strengths (Omegaij) references
omij_data_reference = atomneb.read_omij_references(atom_omij_file)
# read Transition Probabilities (Aij) references
aij_data_reference = atomneb.read_aij_references(atom_aij_file)

atom = 'o'
ion = 'iii'
# read Energy Levels (Ej) of O III upto level number 6
oiii_elj_data = atomneb.read_elj(atom_elj_file, atom, ion, level_num=6)
# print Levels of O III
print(np.asarray(oiii_elj_data.j_v))
# print Energy Levels (cm-1) of O III
print(np.asarray(oiii_elj_data.ej))

# get citations for Energy Levels (Ej) Reference o_iii_elj
citation = atomneb.get_elj_reference_citation(atom_elj_file, 'o_iii_elj')
# print citations for Energy Levels (Ej) Reference o_iii_elj
print(citation)

atom = 'o'
ion = 'iii'
# read Collision Strengths (Omegaij) of O III
oiii_omij_data = atomneb.read_omij(atom_omij_file, atom, ion)
# print Level 1 of Collision Strengths (Omegaij) of O III
print(np.asarray(oiii_omij_data.level1))
# print Level 2 of Collision Strengths (Omegaij) of O III
print(np.asarray(oiii_omij_data.level2))
# print Strength[1] of Collision Strengths (Omegaij) of O III
print(np.asarray(oiii_omij_data.strength)[0])

atom = 'o'
ion = 'iii'
# list all Collision Strengths (Omegaij) data for O III
list_oiii_omij_data = atomneb.search_omij(atom_omij_file, atom, ion)
# print all Collision Strengths (Omegaij) of O III
print(list_oiii_omij_data)

atom = 'o'
ion = 'iii'
# list all Collision Strengths (Omegaij) references for O III
list_oiii_omij_reference = atomneb.list_omij_references(atom_omij_file, atom, ion)
# print all Collision Strengths (Omegaij) References for O III
print(list_oiii_omij_reference)

atom = 'o'
ion = 'iii'
reference = 'CHI52'
# get citations for Collision Strengths (Omegaij) of O III with reference CHI52
citation = atomneb.get_omij_reference_citation(atom_omij_file, atom, ion, reference)
# print citations for Collision Strengths (Omegaij) of O III with reference CHI52
print(citation)

atom = 'o'
ion = 'iii'
reference = 'CHI52'
# read Transition Probabilities (Aij) of O III with reference CHI52
oiii_aij_data = atomneb.read_aij(atom_aij_file, atom, ion)
# print Transition Probabilities (Aij) of O III with reference CHI52
print(oiii_aij_data.aij)

atom = 'o'
ion = 'iii'
reference = 'CHI52'
# get citations for Transition Probabilities (Aij) of O III with reference CHI52
citation = atomneb.get_aij_reference_citation(atom_aij_file, atom, ion, reference)
# print citations for Transition Probabilities (Aij) of O III with reference CHI52
print(citation)

atom = 'o'
ion = 'iii'
# list all Transition Probabilities (Aij) data for O III
list_oiii_aij_data = atomneb.search_aij(atom_aij_file, atom, ion)
# print all Transition Probabilities (Aij) data for O III
print(list_oiii_aij_data)

atom = 'o'
ion = 'iii'
# list all Transition Probabilities (Aij) references for O III
list_oiii_aij_reference = atomneb.list_aij_references(atom_aij_file, atom, ion)
# print all Transition Probabilities (Aij) references for O III
print(list_oiii_aij_reference)


# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('atomic-data', 'chianti60')
atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')


# read Energy Levels (Ej) list
elj_data_list = atomneb.read_elj_list(atom_elj_file)
# read Collision Strengths (Omegaij) list
omij_data_list = atomneb.read_omij_list(atom_omij_file)
# read Transition Probabilities (Aij) list
aij_data_list = atomneb.read_aij_list(atom_aij_file)

# read Energy Levels (Ej) references
elj_data_reference = atomneb.read_elj_references(atom_elj_file)
# read Collision Strengths (Omegaij) references
omij_data_reference = atomneb.read_omij_references(atom_omij_file)
# read Transition Probabilities (Aij) references
aij_data_reference = atomneb.read_aij_references(atom_aij_file)

atom = 'o'
ion = 'iii'
# read Energy Levels (Ej) of O III upto level number 6
oiii_elj_data = atomneb.read_elj(atom_elj_file, atom, ion, level_num=6)
# print Levels of O III
print(np.asarray(oiii_elj_data.j_v))
# print Energy Levels (cm-1) of O III
print(np.asarray(oiii_elj_data.ej))

# get citations for Energy Levels (Ej) Reference o_iii_elj
citation = atomneb.get_elj_reference_citation(atom_elj_file, 'o_iii_elj')
# print citations for Energy Levels (Ej) Reference o_iii_elj
print(citation)

atom = 'o'
ion = 'iii'
# read Collision Strengths (Omegaij) of O III
oiii_omij_data = atomneb.read_omij(atom_omij_file, atom, ion)
# print Level 1 of Collision Strengths (Omegaij) of O III
print(np.asarray(oiii_omij_data.level1))
# print Level 2 of Collision Strengths (Omegaij) of O III
print(np.asarray(oiii_omij_data.level2))
# print Strength[1] of Collision Strengths (Omegaij) of O III
print(np.asarray(oiii_omij_data.strength)[0])

atom = 'o'
ion = 'iii'
# list all Collision Strengths (Omegaij) data for O III
list_oiii_omij_data = atomneb.search_omij(atom_omij_file, atom, ion)
# print all Collision Strengths (Omegaij) of O III
print(list_oiii_omij_data)

atom = 'o'
ion = 'iii'
# list all Collision Strengths (Omegaij) references for O III
list_oiii_omij_reference = atomneb.list_omij_references(atom_omij_file, atom, ion)
# print all Collision Strengths (Omegaij) References for O III
print(list_oiii_omij_reference)

atom = 'o'
ion = 'iii'
reference = 'CHI60'
# get citations for Collision Strengths (Omegaij) of O III with reference CHI60
citation = atomneb.get_omij_reference_citation(atom_omij_file, atom, ion, reference)
# print citations for Collision Strengths (Omegaij) of O III with reference CHI60
print(citation)

atom = 'o'
ion = 'iii'
reference = 'CHI60'
# read Transition Probabilities (Aij) of O III with reference CHI60
oiii_aij_data = atomneb.read_aij(atom_aij_file, atom, ion)
# print Transition Probabilities (Aij) of O III with reference CHI60
print(oiii_aij_data.aij)

atom = 'o'
ion = 'iii'
reference = 'CHI60'
# get citations for Transition Probabilities (Aij) of O III with reference CHI60
citation = atomneb.get_aij_reference_citation(atom_aij_file, atom, ion, reference)
# print citations for Transition Probabilities (Aij) of O III with reference CHI60
print(citation)

atom = 'o'
ion = 'iii'
# list all Transition Probabilities (Aij) data for O III
list_oiii_aij_data = atomneb.search_aij(atom_aij_file, atom, ion)
# print all Transition Probabilities (Aij) data for O III
print(list_oiii_aij_data)

atom = 'o'
ion = 'iii'
# list all Transition Probabilities (Aij) references for O III
list_oiii_aij_reference = atomneb.list_aij_references(atom_aij_file, atom, ion)
# print all Transition Probabilities (Aij) references for O III
print(list_oiii_aij_reference)


# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('atomic-data', 'chianti70')
atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')


# read Energy Levels (Ej) list
elj_data_list = atomneb.read_elj_list(atom_elj_file)
# read Collision Strengths (Omegaij) list
omij_data_list = atomneb.read_omij_list(atom_omij_file)
# read Transition Probabilities (Aij) list
aij_data_list = atomneb.read_aij_list(atom_aij_file)

# read Energy Levels (Ej) references
elj_data_reference = atomneb.read_elj_references(atom_elj_file)
# read Collision Strengths (Omegaij) references
omij_data_reference = atomneb.read_omij_references(atom_omij_file)
# read Transition Probabilities (Aij) references
aij_data_reference = atomneb.read_aij_references(atom_aij_file)

atom = 'o'
ion = 'iii'
# read Energy Levels (Ej) of O III upto level number 6
oiii_elj_data = atomneb.read_elj(atom_elj_file, atom, ion, level_num=6)
# print Levels of O III
print(np.asarray(oiii_elj_data.j_v))
# print Energy Levels (cm-1) of O III
print(np.asarray(oiii_elj_data.ej))

# get citations for Energy Levels (Ej) Reference o_iii_elj
citation = atomneb.get_elj_reference_citation(atom_elj_file, 'o_iii_elj')
# print citations for Energy Levels (Ej) Reference o_iii_elj
print(citation)

atom = 'o'
ion = 'iii'
# read Collision Strengths (Omegaij) of O III
oiii_omij_data = atomneb.read_omij(atom_omij_file, atom, ion)
# print Level 1 of Collision Strengths (Omegaij) of O III
print(np.asarray(oiii_omij_data.level1))
# print Level 2 of Collision Strengths (Omegaij) of O III
print(np.asarray(oiii_omij_data.level2))
# print Strength[1] of Collision Strengths (Omegaij) of O III
print(np.asarray(oiii_omij_data.strength)[0])

atom = 'o'
ion = 'iii'
# list all Collision Strengths (Omegaij) data for O III
list_oiii_omij_data = atomneb.search_omij(atom_omij_file, atom, ion)
# print all Collision Strengths (Omegaij) of O III
print(list_oiii_omij_data)

atom = 'o'
ion = 'iii'
# list all Collision Strengths (Omegaij) references for O III
list_oiii_omij_reference = atomneb.list_omij_references(atom_omij_file, atom, ion)
# print all Collision Strengths (Omegaij) References for O III
print(list_oiii_omij_reference)

atom = 'o'
ion = 'iii'
reference = 'CHI70'
# get citations for Collision Strengths (Omegaij) of O III with reference CHI70
citation = atomneb.get_omij_reference_citation(atom_omij_file, atom, ion, reference)
# print citations for Collision Strengths (Omegaij) of O III with reference CHI70
print(citation)

atom = 'o'
ion = 'iii'
reference = 'CHI70'
# read Transition Probabilities (Aij) of O III with reference CHI70
oiii_aij_data = atomneb.read_aij(atom_aij_file, atom, ion)
# print Transition Probabilities (Aij) of O III with reference CHI70
print(oiii_aij_data.aij)

atom = 'o'
ion = 'iii'
reference = 'CHI70'
# get citations for Transition Probabilities (Aij) of O III with reference CHI70
citation = atomneb.get_aij_reference_citation(atom_aij_file, atom, ion, reference)
# print citations for Transition Probabilities (Aij) of O III with reference CHI70
print(citation)

atom = 'o'
ion = 'iii'
# list all Transition Probabilities (Aij) data for O III
list_oiii_aij_data = atomneb.search_aij(atom_aij_file, atom, ion)
# print all Transition Probabilities (Aij) data for O III
print(list_oiii_aij_data)

atom = 'o'
ion = 'iii'
# list all Transition Probabilities (Aij) references for O III
list_oiii_aij_reference = atomneb.list_aij_references(atom_aij_file, atom, ion)
# print all Transition Probabilities (Aij) references for O III
print(list_oiii_aij_reference)


# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('atomic-data', 'collection')
atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')


# read Energy Levels (Ej) list
elj_data_list = atomneb.read_elj_list(atom_elj_file)
# read Collision Strengths (Omegaij) list
omij_data_list = atomneb.read_omij_list(atom_omij_file)
# read Transition Probabilities (Aij) list
aij_data_list = atomneb.read_aij_list(atom_aij_file)

# read Energy Levels (Ej) references
elj_data_reference = atomneb.read_elj_references(atom_elj_file)
# read Collision Strengths (Omegaij) references
omij_data_reference = atomneb.read_omij_references(atom_omij_file)
# read Transition Probabilities (Aij) references
aij_data_reference = atomneb.read_aij_references(atom_aij_file)

atom = 'o'
ion = 'iii'
# read Energy Levels (Ej) of O III upto level number 6
oiii_elj_data = atomneb.read_elj(atom_elj_file, atom, ion, level_num=6)
# print Levels of O III
print(np.asarray(oiii_elj_data.j_v))
# print Energy Levels (cm-1) of O III
print(np.asarray(oiii_elj_data.ej))

# get citations for Energy Levels (Ej) Reference L7288
citation = atomneb.get_elj_reference_citation(atom_elj_file, 'L7288')
# print citations for Energy Levels (Ej) Reference L7288
print(citation)

atom = 'o'
ion = 'iii'
reference = 'SSB14'
# read Collision Strengths (Omegaij) of O III from Reference SSB14
oiii_omij_data = atomneb.read_omij(atom_omij_file, atom, ion, reference=reference)
# print Level 1 of Collision Strengths (Omegaij) of O III
print(np.asarray(oiii_omij_data.level1))
# print Level 2 of Collision Strengths (Omegaij) of O III
print(np.asarray(oiii_omij_data.level2))
# print Strength[1] of Collision Strengths (Omegaij) of O III
print(np.asarray(oiii_omij_data.strength)[0])

atom = 'o'
ion = 'iii'
# list all Collision Strengths (Omegaij) data for O III
list_oiii_omij_data = atomneb.search_omij(atom_omij_file, atom, ion)
# print all Collision Strengths (Omegaij) of O III
print(list_oiii_omij_data)

atom = 'o'
ion = 'iii'
# list all Collision Strengths (Omegaij) references for O III
list_oiii_omij_reference = atomneb.list_omij_references(atom_omij_file, atom, ion)
# print all Collision Strengths (Omegaij) References for O III
print(list_oiii_omij_reference)

atom = 'o'
ion = 'iii'
reference = 'SSB14'
# get citations for Collision Strengths (Omegaij) of O III with reference SSB14
citation = atomneb.get_omij_reference_citation(atom_omij_file, atom, ion, reference)
# print citations for Collision Strengths (Omegaij) of O III with reference SSB14
print(citation)

atom = 'o'
ion = 'iii'
reference = 'FFT04'
# read Transition Probabilities (Aij) of O III with reference FFT04
oiii_aij_data = atomneb.read_aij(atom_aij_file, atom, ion, reference)
# print Transition Probabilities (Aij) of O III with reference FFT04
print(oiii_aij_data.aij)

atom = 'o'
ion = 'iii'
reference = 'FFT04'
# get citations for Transition Probabilities (Aij) of O III with reference FFT04
citation = atomneb.get_aij_reference_citation(atom_aij_file, atom, ion, reference)
# print citations for Transition Probabilities (Aij) of O III with reference FFT04
print(citation)

atom = 'o'
ion = 'iii'
# list all Transition Probabilities (Aij) data for O III
list_oiii_aij_data = atomneb.search_aij(atom_aij_file, atom, ion)
# print all Transition Probabilities (Aij) data for O III
print(list_oiii_aij_data)

atom = 'o'
ion = 'iii'
# list all Transition Probabilities (Aij) references for O III
list_oiii_aij_reference = atomneb.list_aij_references(atom_aij_file, atom, ion)
# print all Transition Probabilities (Aij) references for O III
print(list_oiii_aij_reference)


# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('atomic-data-rc')
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



# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('atomic-data-rc')
atom_rc_file = os.path.join(base_dir,data_dir, 'rc_PPB91.fits')

atom = 'c'
ion = 'iii' # C II
# read Recombination Coefficients (Aeff) of C II
cii_rc_data = atomneb.read_aeff_ppb91(atom_rc_file, atom, ion)
n_line = len(cii_rc_data.wavelength)
# print information needed for Recombination Coefficients (Aeff) of C II
for i in range(0, n_line):
   print(cii_rc_data.ion[i], cii_rc_data.case1[i], cii_rc_data.wavelength[i],
         cii_rc_data.a[i], cii_rc_data.b[i], cii_rc_data.c[i],
         cii_rc_data.d[i], cii_rc_data.br[i], cii_rc_data.q[i], cii_rc_data.y[i])

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



# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('atomic-data-rc')
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


# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('atomic-data-rc')
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



# Locate datasets
base_dir = os.getcwd()
data_dir = os.path.join('atomic-data-rc')
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

