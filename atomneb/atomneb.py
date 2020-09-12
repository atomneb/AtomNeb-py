# -*- coding: utf-8 -*-

"""
This module contains functions for Atomic Data of Ionized Nebulae
"""

# A. Danehkar
#
# Version 0.2.0, 10/09/2020
# First Release
#

#from astropy.io import fits
#from astropy.table import Table
from astropy.io.fits import getdata
import numpy as np
import pandas as pd

#try:
#    import regex as re
#except ImportError:
#    import re

__all__ = ["read_aij","read_aij_list","search_aij","read_aij_references","list_aij_references",
           "get_aij_reference_citation",
           "read_elj","read_elj_list","read_elj_references","get_elj_reference_citation",
           "read_omij","read_omij_list","search_omij","read_omij_references",
           "list_omij_references","get_omij_reference_citation",
           "read_aeff_collection","read_aeff_collection_list",
           "search_aeff_collection","read_aeff_collection_references",
           "list_aeff_collection_references","get_aeff_collection_reference_citation",
           "read_aeff_sh95","read_aeff_sh95_list","search_aeff_sh95","read_aeff_sh95_references",
           "list_aeff_sh95_references","get_aeff_sh95_reference_citation",
           "read_aeff_ppb91","read_aeff_ppb91_list","search_aeff_ppb91","read_aeff_ppb91_references",
           "list_aeff_ppb91_references","get_aeff_ppb91_reference_citation",
           "read_aeff_he_i_pfsd12","read_aeff_he_i_pfsd12_list","search_aeff_he_i_pfsd12",
           "read_aeff_he_i_pfsd12_references","list_aeff_he_i_pfsd12_references",
           "get_aeff_he_i_pfsd12_reference_citation",
           "read_aeff_n_ii_fsl13","read_aeff_n_ii_fsl13_list","search_aeff_n_ii_fsl13",
           "read_aeff_n_ii_fsl13_references","list_aeff_n_ii_fsl13_references",
           "get_aeff_n_ii_fsl13_reference_citation",
           "read_aeff_o_ii_ssb17","read_aeff_o_ii_ssb17_list","search_aeff_o_ii_ssb17",
           "read_aeff_o_ii_ssb17_references","list_aeff_o_ii_ssb17_references",
           "get_aeff_o_ii_ssb17_reference_citation"]

def read_aij(atom_aij_file, atom, ion, reference=None, level_num=None):
    """
         This function returns the transition probabilities (Aij) from the table extensions
         of the FITS data file ('AtomAij.fits').

     :Returns:
        type=an array of data. This function returns the aij_data:
              { Aij:dblarr(n_level,n_level) }.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('AtoAij.fits')
         atom          : in, required, type=string
                         atom name e.g. 'o'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Keywords:
         reference     : in, type=string
                         set for the reference,  not necessary
         level_num     : in, type=string
                         set for the maximum level number.

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data', 'collection')
         >>> atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
         >>> atom='o'
         >>> ion='iii'
         >>> reference='FFT04'
         >>> oiii_aij_data = atomneb.read_aij(atom_aij_file, atom, ion, reference)
         >>> print(oiii_aij_data.aij)
            0.0000000   2.5960000e-05   3.0300000e-11   2.3220000e-06       0.0000000    0.0021910000
            0.0000000       0.0000000   9.6320000e-05    0.0069510000      0.22550000       230.80000
            0.0000000       0.0000000       0.0000000     0.020290000   0.00069980000       576.50000
            0.0000000       0.0000000       0.0000000       0.0000000       1.6850000    0.0057770000
            0.0000000       0.0000000       0.0000000       0.0000000       0.0000000   3.7600000e-11
            0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000

     :Categories:
       Collisionally Excited Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         24/12/2015, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aij_list(atom_aij_file)
    if (reference is not None) == 1:
        atom_ion_name = atom.lower() + '_' + ion.lower() + '_aij_' + reference.upper()
        aij_list = np.asarray(element_data_list.aij_data)
        ii = np.where(aij_list== atom_ion_name)[0]  #
        if len(ii) == 0:
            print ('could not find the given element or ion')
            return 0
    else:
        atom_ion_name = atom.lower() + '_' + ion.lower() + '_aij'
        aij_list=np.asarray(element_data_list.aij_data)
        #aij_res = [x for x in aij_list if re.search(atom_ion_name, x)]
        aij_res = [x for x in aij_list if atom_ion_name in x]
        iii = np.where(aij_list == aij_res[0])[0]  #
        if len(iii) == 0:
            print('could not find the given element or ion')
            return 0
        ii = min(iii)
    extension = np.asarray(element_data_list.extension)[ii]

    aij, header1 = getdata(atom_aij_file, int(extension), header=True)
    n_level = len(aij)
    #aij_temp=np.zeros((n_level,n_level))
    #aij_template = {'aij': aij_temp }
    aij_template={'aij':np.zeros((n_level,n_level))}
    if (level_num is not None) == 1:
        if level_num <= n_level:
            level_length = level_num
    else:
        level_length = n_level
    aij_data = pd.concat([pd.DataFrame(v) for k, v in aij_template.items()], axis=1, keys=list(aij_template.keys()))
    #aij_data = pd.DataFrame(data=aij_template, index=np.arange(1))
    aij_data.aij = aij
    return aij_data

def read_aij_list(atom_aij_file):
    """
         This function returns the list of transition probabilities (Aij) from the 1st binary table extension
         of the FITS data file ('AtomAij.fits').

     :Private:

     :Returns:
        type=an array of data. This function returns the aij_data_list:
              { Aij_Data:'',
                Extension:0.0}

     :Params:
         Atom_Aij_file  : in, required, type=string
                         the FITS data file name ('AtomAij.fits')

     :Categories:
       Collisionally Excited Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         24/12/2015, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_template = {'aij_data': '', 'extension': 0 }
    data, hdr = getdata(atom_aij_file, 1, header=True)
    aij_data=data.AIJ_DATA
    extension = data.EXTENSION
    element_length = len(aij_data)

    element_data = pd.DataFrame(data=element_template, index=np.arange(element_length))
    #aij_data=aij_data.strip
    aij_data=np.char.strip(aij_data)
    element_data.aij_data=aij_data
    element_data.extension = extension
    return element_data

def search_aij(atom_aij_file, atom, ion):
    """
         This function searches transition probabilities (Aij) for given element
         and ionic levels in the FITS data file ('AtomAij.fits'), and returns the data entry.

     :Returns:
        type=array of data. This function returns the Aij_Data.

     :Params:
         Atom_Aij_file  : in, required, type=string
                         the FITS data file name ('AtoAij.fits')
         atom          : in, required, type=string
                         atom name e.g. 'o'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data', 'collection')
         >>> atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
         >>> atom='o'
         >>> ion='iii'
         >>> list_oiii_aij_data = atomneb.search_aij(atom_aij_file, atom, ion)
         >>> print(list_oiii_aij_data)
            o_iii_aij_FFT04-SZ00 o_iii_aij_FFT04 o_iii_aij_GMZ97-WFD96 o_iii_aij_SZ00-WFD96

     :Categories:
       Collisionally Excited Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         24/12/2015, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aij_list(atom_aij_file)
    atom_ion_name = atom.lower() + '_' + ion.lower() + '_aij_'
    aij_list = np.asarray(element_data_list.aij_data)
    #aij_res = [x for x in aij_list if re.search(atom_ion_name, x)]
    aij_res = [x for x in aij_list if atom_ion_name in x]
    if len(aij_res) == 0:
        print('could not find the given element or ion')
        return 0
    select_aij_data = aij_res
    return select_aij_data

def read_aij_references(atom_aij_file):
    """
         This function returns the reference list of transition probabilities (Aij) from the 1nd binary table extension
         of the FITS data file ('AtomAij.fits').

     :Private:

     :Returns:
        type=an array of data. This function returns the aij_data_reference:
              { Reference:'',
                Citation:''}

     :Params:
         Atom_Aij_file  : in, required, type=string
                         the FITS data file name ('AtomAij.fits')

     :Categories:
       Collisionally Excited Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         24/12/2015, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    reference_template = {'reference': '', 'citation': 0}
    data, hdr = getdata(atom_aij_file, 2, header=True)
    reference = data.ATOMICDATA
    citation = data.REFERENCE
    reference_length = len(reference)

    reference_data = pd.DataFrame(data=reference_template, index=np.arange(reference_length))
    # aij_data=aij_data.strip
    reference = np.char.strip(reference)
    citation = np.char.strip(citation)
    reference_data.reference = reference
    reference_data.citation = citation
    return reference_data

def list_aij_references(atom_aij_file, atom, ion):
    """
         This function returns a list for all references of transition probabilities (Aij)
         for given element and ionic level from the FITS data file ('AtoAij.fits').

     :Returns:
        type=an array of strings. This function returns the references.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('AtoAij.fits')
         atom          : in, required, type=string
                         atom name e.g. 'o'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data', 'collection')
         >>> atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
         >>> atom='o'
         >>> ion='iii'
         >>> list_oiii_aij_reference = atomneb.list_aij_references(atom_aij_file, atom, ion)
         >>> print(list_oiii_aij_reference)
            FFT04-SZ00 FFT04 GMZ97-WFD96 SZ00-WFD96

     :Categories:
       Collisionally Excited Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aij_list(atom_aij_file)

    atom_ion_name = atom.lower() + '_' + ion.lower() + '_aij_'
    aij_list = np.asarray(element_data_list.aij_data)
    #aij_res = [x for x in aij_list if re.search(atom_ion_name, x)]
    aij_res = [x for x in aij_list if atom_ion_name in x]
    ii_length = len(aij_res)
    if ii_length == 0:
        print('could not find the given element or ion')
        return 0
    nlines1 = ii_length
    references = ['' for i in range(nlines1)]
    for i in range(0, nlines1):
        select_aij_data_str1 = aij_res[i].split('_')
        references[i] = select_aij_data_str1[3]
    return references

def get_aij_reference_citation(atom_aij_file, atom, ion, reference):
    """
         This function returns the reference citation for a transition probability (Aij)
         from the 2nd binary table extension of the FITS data file ('AtoAij.fits')

     :Returns:
        type=string. This function returns the Citation.

     :Params:
         Atom_Aij_file : in, required, type=string
                         the FITS data file name ('AtoAij.fits')
         atom          : in, required, type=string
                         atom name e.g. 'o'
         ion           : in, required, type=string
                         ionic level e.g 'iii'
         reference     : in, type=string
                         set for the reference e.g. 'FFT04'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data', 'collection')
         >>> atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
         >>> atom='o'
         >>> ion='iii'
         >>> reference='FFT04'
         >>> citation = atomneb.get_aij_reference_citation(atom_aij_file, atom, ion, reference)
         >>> print(citation)
            Froese Fischer et al 2004, ADNDT 87, 1

     :Categories:
       Collisionally Excited Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         24/12/2015, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_reference = read_aij_references(atom_aij_file)
    reference_list = np.asarray(element_data_reference.reference)
    atom_ion_name = atom.lower() + '_' + ion.lower() + '_aij_' + reference.upper()
    ii = np.where(reference_list == atom_ion_name)[0]  #
    if len(ii) == 0:
        print('could not find the given reference')
        return 0
    citation = np.asarray(element_data_reference.citation[ii])

    return citation

def read_elj(atom_elj_file, atom, ion, level_num=None):
    """
         This function returns the energy levels (Ej) from the table extensions
         of the FITS data file ('AtomElj.fits').

     :Returns:
        type=an array of data. This function returns the elj_data:
              { Configuration:'',
                Term:'',
                J:'',
                J_v:0.0,
                Ej:0.0,
                Reference:''}.

     :Params:
         atom_elj_file  : in, required, type=string
                         the FITS data file name ('AtomElj.fits')
         atom          : in, required, type=string
                         atom name e.g. 'o'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Keywords:
         level_num     : in, type=string
                         set for the maximum level number.

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data', 'collection')
         >>> atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
         >>> atom='o'
         >>> ion='iii'
         >>> oiii_elj_data=atomneb.read_elj(atom_elj_file, atom, ion, level_num=6)
         >>> print(np.asarray(oiii_elj_data.j_v))
            0.00000      1.00000      2.00000      2.00000      0.00000      2.00000
         >>> print(np.asarray(oiii_elj_data.ej))
            0.0000000       113.17800       306.17400       20273.270       43185.740       60324.790

     :Categories:
       Collisionally Excited Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         24/12/2015, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_elj_list(atom_elj_file)

    atom_ion_name = atom.lower() + '_' + ion.lower() + '_elj'
    elj_list = np.asarray(element_data_list.elj_data)
    ii = np.where(elj_list == atom_ion_name)[0]  #
    if len(ii) == 0:
        print('could not find the given element or ion')
        return 0
    extension = np.asarray(element_data_list.extension)[ii]

    # level_template={Configuration:'', Term:'', J:'', J_v: float(0.0), Ej:double(0.0), Reference:''}
    level_template = {'configuration':'', 'term':'', 'j':'', 'j_v': float(0.0), 'ej':float(0.0), 'reference':''}

    data, header1 = getdata(atom_elj_file, int(extension), header=True)
    configuration=data.CONFIGURATION
    term=data.TERM
    j_s=data.J
    j_v=data.J_V
    ej=data.EJ
    reference=data.REFERENCE
    level_length = len(reference)
    if (level_num is not None) == 1:
        if level_num <= level_length:
            level_length = level_num

    level_data = pd.DataFrame(data=level_template, index=np.arange(level_length))
    term=np.char.strip(term)
    j_s=np.char.strip(j_s)
    reference=np.char.strip(reference)

    for i in range(0, level_length):
        level_data.loc[i, 'configuration'] = configuration[i]
        level_data.loc[i, 'term'] = term[i]
        level_data.loc[i, 'j_s'] = j_s[i]
        level_data.loc[i, 'j_v'] = j_v[i]
        level_data.loc[i, 'ej'] = ej[i]
        level_data.loc[i, 'reference'] = reference[i]

    return level_data

def read_elj_list(atom_elj_file):
    """
         This function returns the list of energy levels (Ej) from the 1st binary table extension
         of the FITS data file ('AtomElj.fits')

     :Private:

     :Returns:
        type=an array of data. This function returns the elj_data_list:
              { Elj_Data:'',
                Extension:0.0}

     :Params:
         atom_elj_file  : in, required, type=string
                         the FITS data file name ('AtomElj.fits')

     :Categories:
       Collisionally Excited Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         24/12/2015, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_template = {'elj_data': '', 'extension': 0 }
    data, hdr = getdata(atom_elj_file, 1, header=True)
    elj_data=data.EJ_DATA
    extension = data.EXTENSION
    element_length = len(elj_data)

    element_data = pd.DataFrame(data=element_template, index=np.arange(element_length))
    #aij_data=aij_data.strip
    elj_data=np.char.strip(elj_data)
    element_data.elj_data=elj_data
    element_data.extension = extension
    return element_data

def read_elj_references(atom_elj_file):
    """
         This function returns the reference list of energy levels (Ej) from the 2nd binary table extension
         of the FITS data file ('AtomElj.fits').

     :Private:

     :Returns:
        type=an array of data. This function returns the aij_data_reference:
              { Reference:'',
                Citation:''}

     :Params:
         atom_elj_file  : in, required, type=string
                          the FITS data file name ('AtomElj.fits')

     :Categories:
       Collisionally Excited Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         24/12/2015, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    reference_template = {'reference': '', 'citation': 0}
    data, hdr = getdata(atom_elj_file, 2, header=True)
    reference = data.ATOMICDATA
    citation = data.REFERENCE
    reference_length = len(reference)

    reference_data = pd.DataFrame(data=reference_template, index=np.arange(reference_length))


    reference = np.char.strip(reference)
    citation = np.char.strip(citation)
    reference_data.reference = reference
    reference_data.citation = citation
    return reference_data

def get_elj_reference_citation(atom_elj_file, reference):
    """
         This function returns the reference citation for energy levels (Ej)
         from the 2nd binary table extension of the FITS data file ('AtomElj.fits').

     :Returns:
        type=string. This function returns the Citation.

     :Params:
         atom_elj_file  : in, required, type=string
                         the FITS data file name ('AtomElj.fits')
         reference     : in, type=string
                         set for the reference e.g. 'L7288'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data', 'collection')
         >>> atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
         >>> reference='L7288'
         >>> citation=atomneb.get_elj_reference_citation(atom_elj_file, reference)
         >>> print(citation)
            C. E. Moore, in CRC Series in Evaluated Data in Atomic Physics, 339 pp. (CRC Press, Boca Raton, FL, 1993)

     :Categories:
       Collisionally Excited Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         24/12/2015, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_reference = read_elj_references(atom_elj_file)
    reference_list = np.asarray(element_data_reference.reference)
    ii = np.where(reference_list == reference)[0]  #
    if len(ii) == 0:
        print('could not find the given reference')
        return 0
    citation = np.asarray(element_data_reference.citation[ii])

    return citation

def read_omij(atom_omij_file, atom, ion, reference=None, level_num=None):
    """
         This function returns the collision strengths (omega_ij) from the table extensions
         of the FITS data file ('AtomOmij.fits').

     :Returns:
        type=an array of data. This function returns the omij_data:
              { level1:0,
                level2:0,
                strength:dblarr(temp_steps)}.

     :Params:
         atom_omij_file  : in, required, type=string
                         the FITS data file name ('AtomOmij.fits')
         atom          : in, required, type=string
                         atom name e.g. 'o'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Keywords:
         reference     : in, type=string
                         set for the reference e.g. 'SSB14'
         level_num     : in, type=string
                         set for the maximum level number.

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data', 'collection')
         >>> atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
         >>> atom='o'
         >>> ion='iii'
         >>> reference='SSB14'
         >>> oiii_omij_data=atomneb.read_omij(atom_omij_file, atom, ion, reference=reference)
         >>> print(np.asarray(oiii_omij_data.level1))
            0       1       1       1       1       2       2       2       3       3       4
         >>> print(np.asarray(oiii_omij_data.level2))
            0       2       3       4       5       3       4       5       4       5       5
         >>> print(np.asarray(oiii_omij_data.strength)[0])
            100.00000       125.89254       158.48932       199.52623       251.18864       ...

     :Categories:
       Collisionally Excited Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         24/12/2015, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_omij_list(atom_omij_file)
    if (reference is not None) == 1:
        atom_ion_name = atom.lower() + '_' + ion.lower() + '_omij_' + reference.upper()
        omij_list = np.asarray(element_data_list.omij_data)
        ii = np.where(omij_list== atom_ion_name)[0]  #
        if len(ii)== 0:
            print ('could not find the given element or ion')
            return 0
    else:
        atom_ion_name = atom.lower() + '_' + ion.lower() + '_omij'
        omij_list=np.asarray(element_data_list.omij_data)
        #omij_res = [x for x in omij_list if re.search(atom_ion_name, x)]
        omij_res = [x for x in omij_list if atom_ion_name in x]
        iii = np.where(omij_list == omij_res[0])[0]  #
        if len(iii) == 0:
            print('could not find the given element or ion')
            return 0
        ii = min(iii)
    extension = np.asarray(element_data_list.extension)[ii]

    data, header1 = getdata(atom_omij_file, int(extension), header=True)
    level1=data.LEVEL1
    level2=data.LEVEL2
    strength=data.STRENGTH
    n_line = len(level1)
    temp_steps=len(strength)
    n_level = max([max(level1), max(level2)])
    # omij_template={level1:0, level2:0, strength:dblarr(temp_steps)}
    strength_temp = np.zeros(temp_steps)
    omij_template = {'level1':0, 'level2':0, 'strength': [strength_temp] }
    if (level_num is not None) == 1:
        if level_num <= n_level:
            level_length = level_num
    else:
        level_length = n_level
    #omij_data = pd.concat([pd.DataFrame(v) for k, v in omij_template.items()], axis=3, keys=list(omij_template.keys()))
    omij_data = pd.DataFrame(data=omij_template, index=np.arange(n_line))
    #for i in range(0, n_line):
    #    omij_data.loc[i, 'level1'] = level1[i]
    #    omij_data.loc[i, 'level2'] = level2[i]
    #    omij_data.loc[i, 'strength'] = strength[:,i]
    for i in range(0, n_line):
        omij_data.loc[i, 'level1'] = level1[i]
        omij_data.loc[i, 'level2'] = level2[i]
        #print(strength[i,:])
        #omij_data.strength[i] = strength[i,:]
        omij_data['strength'][:][i] = strength[i,:]
        #omij_data['strength'][i] = strength[i, :]
        #omij_data.loc[i, 'strength'] = strength[i,:]
    #omij_data.level1 = level1
    #omij_data.level2 = level2
    # omij_data.strength = strength
    return omij_data

def read_omij_list(atom_omij_file):
    """
         This function returns the list of collision strengths (omega_ij) from the 1st binary table extension
         of the FITS data file ('AtomOmij.fits').

     :Private:

     :Returns:
        type=an array of data. This function returns the omij_data_list:
              { Omij_Data:'',
                Extension:0.0}

     :Params:
         atom_omij_file  : in, required, type=string
                         the FITS data file name ('AtomOmij.fits')

     :Categories:
       Collisionally Excited Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         24/12/2015, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_template = {'omij_data': '', 'extension': 0 }
    data, hdr = getdata(atom_omij_file, 1, header=True)
    omij_data=data.OMIJ_DATA
    extension = data.EXTENSION
    element_length = len(omij_data)

    element_data = pd.DataFrame(data=element_template, index=np.arange(element_length))
    #omij_data=omij_data.strip
    omij_data=np.char.strip(omij_data)
    element_data.omij_data=omij_data
    element_data.extension = extension
    return element_data

def search_omij(atom_omij_file, atom, ion):
    """
         This function searches collision strengths (omega_ij) for given element
         and ionic levels in the FITS data file ('AtomOmij.fits'), and returns the data entry.

     :Returns:
        type=array of data. This function returns the Omij_Data.

     :Params:
         atom_omij_file  : in, required, type=string
                         the FITS data file name ('AtomOmij.fits')
         atom          : in, required, type=string
                         atom name e.g. 'o'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data', 'collection')
         >>> atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
         >>> atom='o'
         >>> ion='iii'
         >>> list_oiii_omij_data = atomneb.search_omij(atom_omij_file, atom, ion)
         >>> print(list_oiii_omij_data)
            o_iii_omij_AK99 o_iii_omij_LB94 o_iii_omij_Pal12-AK99 o_iii_omij_SSB14

     :Categories:
       Collisionally Excited Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         24/12/2015, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_omij_list(atom_omij_file)
    atom_ion_name = atom.lower() + '_' + ion.lower() + '_omij_'
    omij_list = np.asarray(element_data_list.omij_data)
    #omij_res = [x for x in omij_list if re.search(atom_ion_name, x)]
    omij_res = [x for x in omij_list if atom_ion_name in x]
    if len(omij_res) == 0:
        print('could not find the given element or ion')
        return 0
    select_omij_data = omij_res
    return select_omij_data

def read_omij_references(atom_omij_file):
    """
         This function returns the reference list of collision strengths (omega_ij) from the 2nd binary table extension
         of the FITS data file ('AtomOmij.fits').

     :Private:

     :Returns:
        type=an array of data. This function returns the aij_data_reference:
              { Reference:'',
                Citation:''}

     :Params:
         atom_omij_file  : in, required, type=string
                           the FITS data file name ('AtomOmij.fits')

     :Categories:
       Collisionally Excited Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         24/12/2015, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    reference_template = {'reference': '', 'citation': 0}
    data, hdr = getdata(atom_omij_file, 2, header=True)
    reference = data.ATOMICDATA
    citation = data.REFERENCE
    reference_length = len(reference)

    reference_data = pd.DataFrame(data=reference_template, index=np.arange(reference_length))
    # aij_data=aij_data.strip
    reference = np.char.strip(reference)
    citation = np.char.strip(citation)
    reference_data.reference = reference
    reference_data.citation = citation
    return reference_data

def list_omij_references(atom_omij_file, atom, ion):
    """
         This function returns a list for all references of collision strengths (Omega_ij)
         for given element and ionic level from the FITS data file ('AtomOmij.fits').

     :Returns:
        type=an array of strings. This function returns the references.

     :Params:
         atom_omij_file  : in, required, type=string
                         the FITS data file name ('AtomOmij.fits')
         atom          : in, required, type=string
                         atom name e.g. 'o'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data', 'collection')
         >>> atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
         >>> atom='o'
         >>> ion='iii'
         >>> list_oiii_omij_reference = atomneb.list_omij_references(atom_omij_file, atom, ion)
         >>> print(list_oiii_omij_reference)
            AK99 LB94 Pal12-AK99 SSB14

     :Categories:
       Collisionally Excited Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_omij_list(atom_omij_file)
    atom_ion_name = atom.lower() + '_' + ion.lower() + '_omij_'
    omij_list = np.asarray(element_data_list.omij_data)
    #omij_res = [x for x in omij_list if re.search(atom_ion_name, x)]
    omij_res = [x for x in omij_list if atom_ion_name in x]
    ii_length = len(omij_res)
    if ii_length == 0:
        print('could not find the given element or ion')
        return 0
    nlines1 = ii_length
    references = ['' for i in range(nlines1)]
    for i in range(0, nlines1):
        select_omij_data_str1 = omij_res[i].split('_')
        references[i] = select_omij_data_str1[3]
    return references

def get_omij_reference_citation(atom_omij_file, atom, ion, reference):
    """
        This function returns the reference citation for collision strengths (Omega_ij)
        from the 2nd binary table extension of the FITS data file ('AtomOmij.fits').

    :Returns:
       type=string. This function returns the Citation.

    :Params:
        atom_omij_file  : in, required, type=string
                        the FITS data file name ('AtomOmij.fits')
        atom          : in, required, type=string
                        atom name e.g. 'o'
        ion           : in, required, type=string
                        ionic level e.g 'iii'
        reference     : in, type=string
                        set for the reference e.g. 'SSB14'

    :Examples:
       For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
        >>> base_dir = os.getcwd()
        >>> data_dir = os.path.join('..','atomic-data', 'collection')
        >>> atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        >>> atom='o'
        >>> ion='iii'
        >>> reference='SSB14'
        >>> citation = atomneb.get_omij_reference_citation(atom_omij_file, atom, ion, reference)
        >>> print(citation)
           Storey, P. J., Sochi, T., and Badnell, N. R. 2014, Astron.Astrophys., 441, 3028

    :Categories:
      Collisionally Excited Lines

    :Dirs:
     ./
         Main routines

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.2.0

    :History:
        24/12/2015, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_reference = read_omij_references(atom_omij_file)
    reference_list = np.asarray(element_data_reference.reference)
    atom_ion_name = atom.lower() + '_' + ion.lower() + '_omij_' + reference.upper()
    ii = np.where(reference_list == atom_ion_name)[0]  #
    if len(ii) == 0:
        print('could not find the given reference')
        return 0
    citation = np.asarray(element_data_reference.citation[ii])

    return citation

def read_aeff_collection(atom_rc_file, atom, ion, br=None, reference=None):
    """
         This function returns the effective recombination coefficients (Aeff) from the table extensions
         of the FITS data file ('rc_collection.fits').

     :Returns:
        type=an array of data. This function returns the effective recombination coefficients.
               aeff_data (c_iii_aeff)
               { Wavelength:0.0,
                 a: 0.0,
                 b: 0.0,
                 c: 0.0,
                 d: 0.0,
                 f: 0.0}

               aeff_data (n_iii_aeff)
               { a: 0.0,
                 b: 0.0,
                 c: 0.0}

               aeff_data (n_iii_br)
               {Wavelength: 0.0,
                BR: 0.0, $
                g1:0,
                g2:0,
                Mult1:'',
                LowerTerm:'',
                UpperTerm:'' }

               aeff_data (o_iii_aeff)
               {Term: '',
                Case1: '',
                a2: 0.0,
                a4: 0.0,
                a5: 0.0,
                a6: 0.0,
                b: 0.0,
                c: 0.0,
                d: 0.0}

               aeff_data (o_iii_br)
               {Wavelength:double(0.0),
                Br_A: 0.0,
                Br_B: 0.0,
                Br_C: 0.0,
                g1: 0,
                g2: 0,
                Mult1: '',
                LowerTerm: '',
                UpperTerm: ''}

                aeff_data (ne_iii_aeff)
                {Wavelength:0.0,
                 a: 0.0,
                 b: 0.0,
                 c: 0.0,
                 d: 0.0,
                 f: 0.0,
                 br: 0.0}

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_collection.fits')
         atom          : in, required, type=string
                         atom name e.g. 'c'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Keywords:
         br            : in, type=boolean
                         set for the branching ratios (Br), may not necessary
         reference     : in, type=string
                         set for the reference, not necessary

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data-rc')
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_collection.fits')
         >>> atom='c'
         >>> ion='iii' # C III
         >>> cii_rc_data = atomneb.read_aeff_collection(atom_rc_file, atom, ion)
         >>> n_line = len(cii_rc_data.wavelength)
         >>> for i in range(0, n_line):
         >>>       print(cii_rc_data.wavelength[i], cii_rc_data.a[i],
         >>>             cii_rc_data.b[i], cii_rc_data.c[i],
         >>>             cii_rc_data.d[i], cii_rc_data.f[i])
            914.00000      0.69280000     0.021400000    -0.016300000     -0.24310000     -0.88000000
            962.00000       1.0998000   -0.0042000000    -0.027900000     -0.22940000     -0.96560000
            ...

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aeff_collection_list(atom_rc_file)

    if (br is not None) == 1:
        prefix = '_br'
    else:
        prefix = '_aeff'
    if (reference is not None) == 1:
        atom_ion_name = atom.lower() + '_' + ion.lower() + prefix + reference.upper()
        aeff_list = np.asarray(element_data_list.aeff_data)
        ii = np.where(aeff_list== atom_ion_name)[0]  #
        if len(ii) == 0:
            print ('could not find the given element or ion')
            return 0
    else:
        atom_ion_name = atom.lower() + '_' + ion.lower() + prefix
        aeff_list=np.asarray(element_data_list.aeff_data)
        #aeff_res = [x for x in aeff_list if re.search(atom_ion_name, x)]
        aeff_res = [x for x in aeff_list if atom_ion_name in x]
        iii = np.where(aeff_list == aeff_res[0])[0]  #
        if len(iii) == 0:
            print('could not find the given element or ion')
            return 0
        ii = min(iii)
    extension = np.asarray(element_data_list.extension)[ii]
    atom_ion_name = atom.lower() + '_' + ion.lower() + prefix

    data, header1 = getdata(atom_rc_file, extension, header=True)

    if atom_ion_name == 'c_iii_aeff':
        rc_template={'wavelength': float(0.0), 'a': float(0.0), 'b':float(0.0), 'c': float(0.0), 'd':float(0.0), 'f':float(0.0)}
        wavelength=data.WAVELENGTH
        a = data.A
        b = data.B
        c = data.C
        d = data.D
        f = data.F
        n_line = len(wavelength)
        rc_data = pd.DataFrame(data=rc_template, index=np.arange(n_line))
        rc_data.wavelength = wavelength
        rc_data.a = a
        rc_data.b = b
        rc_data.c = c
        rc_data.d = d
        rc_data.f = f
    elif atom_ion_name == 'n_iii_aeff':
        rc_template={'a': float(0.0), 'b':float(0.0), 'c': float(0.0)}
        a = data.A
        b = data.B
        c = data.C
        n_line = len(a)
        rc_data = pd.DataFrame(data=rc_template, index=np.arange(n_line))
        rc_data.a = a
        rc_data.b = b
        rc_data.c = c
    elif atom_ion_name == 'n_iii_br':
        rc_template={'wavelength':float(0.0), 'br':float(0.0),
                      'g1':int(0), 'g2':int(0),
                      'mult1':'', 'lowerterm':'', 'upperterm':''
                      }
        wavelength=data.WAVELENGTH
        br = data.BR
        g1 = data.G1
        g2 = data.G2
        mult1 = data.MULT1
        lowerterm = data.LOWERTERM
        upperterm = data.UPPERTERM
        mult1 = np.char.strip(mult1)
        lowerterm = np.char.strip(lowerterm)
        upperterm = np.char.strip(upperterm)
        n_line = len(wavelength)
        rc_data = pd.DataFrame(data=rc_template, index=np.arange(n_line))
        rc_data.wavelength = wavelength
        rc_data.br = br
        rc_data.g1 = g1
        rc_data.g2 = g2
        rc_data.mult1 = mult1
        rc_data.lowerterm = lowerterm
        rc_data.upperterm = upperterm
    elif atom_ion_name == 'o_iii_aeff':
        rc_template={'term':'', 'case1': '',
                     'a2':float(0.0), 'a4':float(0.0), 'a5':float(0.0), 'a6':float(0.0),
                     'b':float(0.0), 'c':float(0.0), 'd':float(0.0)}
        term=data.TERM
        case1 = data.CASE1
        a2 = data.A2
        a4 = data.A4
        a5 = data.A5
        a6 = data.A6
        b = data.B
        c = data.C
        d = data.D
        term = np.char.strip(term)
        case1 = np.char.strip(case1)
        n_line = len(a2)
        rc_data = pd.DataFrame(data=rc_template, index=np.arange(n_line))
        rc_data.term = term
        rc_data.case1 = case1
        rc_data.a2 = a2
        rc_data.a4 = a4
        rc_data.a5 = a5
        rc_data.a6 = a6
        rc_data.b = b
        rc_data.c = c
        rc_data.d = d
    elif atom_ion_name == 'o_iii_br':
        rc_template={'wavelength':float(0.0),
                     'br_a':float(0.0), 'br_b':float(0.0), 'br_c':float(0.0),
                     'g1':int(0), 'g2':int(0),
                     'mult1':'', 'lowerterm':'', 'upperterm':'' }
        wavelength=data.WAVELENGTH
        br_a = data.BR_A
        br_b = data.BR_B
        br_c = data.BR_C
        g1 = data.G1
        g2 = data.G2
        mult1 = data.MULT1
        lowerterm = data.LOWERTERM
        upperterm = data.UPPERTERM
        mult1 = np.char.strip(mult1)
        lowerterm = np.char.strip(lowerterm)
        upperterm = np.char.strip(upperterm)
        n_line = len(br_a)
        rc_data = pd.DataFrame(data=rc_template, index=np.arange(n_line))
        rc_data.wavelength = wavelength
        rc_data.br_a = br_a
        rc_data.br_b = br_b
        rc_data.br_c = br_c
        rc_data.g1 = g1
        rc_data.g2 = g2
        rc_data.mult1 = mult1
        rc_data.lowerterm = lowerterm
        rc_data.upperterm = upperterm
    elif atom_ion_name == 'ne_iii_aeff':
        rc_template={'wavelength': float(0.0), 'a': float(0.0), 'b':float(0.0), 'c': float(0.0), 'd':float(0.0), 'f':float(0.0), 'br':float(0.0)}
        wavelength=data.WAVELENGTH
        a = data.A
        b = data.B
        c = data.C
        d = data.D
        f = data.F
        br = data.BR
        n_line = len(wavelength)
        rc_data = pd.DataFrame(data=rc_template, index=np.arange(n_line))
        rc_data.wavelength = wavelength
        rc_data.a = a
        rc_data.b = b
        rc_data.c = c
        rc_data.d = d
        rc_data.f = f
        rc_data.br = br
    else:
        raise RuntimeError('no match found for expression')
    # endfor
    return rc_data

def read_aeff_collection_list(atom_rc_file):
    """
         This function returns the list of effective recombination coefficients (Aeff) from the 1st binary table extension
         of the FITS data file ('rc_collection.fits')

     :Private:

     :Returns:
        type=an array of data. This function returns the aeff_data_list:
               { Aeff_Data:'',
                 Extension:0.0}

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_collection.fits')

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_template = {'aeff_data': '', 'extension': 0 }
    data, hdr = getdata(atom_rc_file, 1, header=True)
    aeff_data=data.AEFF_DATA
    extension = data.EXTENSION
    element_length = len(aeff_data)

    element_data = pd.DataFrame(data=element_template, index=np.arange(element_length))
    #aij_data=aij_data.strip
    aeff_data=np.char.strip(aeff_data)
    element_data.aeff_data=aeff_data
    element_data.extension = extension

    return element_data

def search_aeff_collection(atom_rc_file, atom, ion, br=None):
    """
         This function searches effective recombination coefficients (Aeff) for given element
         and ionic levels in the FITS data file ('rc_collection.fits'), and returns the data entry.

     :Returns:
        type=array of data. This function returns the Aeff_Data.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_collection.fits')
         atom          : in, required, type=string
                         atom name e.g. 'c'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Keywords:
         br            : in, type=boolean
                         set for the branching ratios (Br), may not necessary

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data-rc')
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_collection.fits')
         >>> atom='c'
         >>> ion='iii' # C III
         >>> list_cii_aeff_data = atomneb.search_aeff_collection(atom_rc_file, atom, ion)
         >>> print(list_cii_aeff_data)
            c_iii_aeff

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    if (br is not None) == 1:
        prefix = '_br'
    else:
        prefix = '_aeff'

    element_data_list = read_aeff_collection_list(atom_rc_file)
    atom_ion_name = atom.lower() + '_' + ion.lower() + prefix
    aeff_list = np.asarray(element_data_list.aeff_data)
    #aeff_res = [x for x in aeff_list if re.search(atom_ion_name, x)]
    aeff_res = [x for x in aeff_list if atom_ion_name in x]
    if len(aeff_res) == 0:
        print('could not find the given element or ion')
        return 0
    select_aeff_data = aeff_res
    return select_aeff_data

def read_aeff_collection_references(atom_rc_file):
    """
         This function returns the reference list of recombination coefficients (Aeff) from the 2nd binary table extension
         of the FITS data file ('rc_collection.fits').

     :Private:

     :Returns:
        type=an array of data. This function returns the aeff_data_reference:
              { Reference:'',
                Citation:''}

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_collection.fits')

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    reference_template = {'reference': '', 'citation': 0}
    data, hdr = getdata(atom_rc_file, 2, header=True)
    reference = data.ATOMICDATA
    citation = data.REFERENCE
    reference_length = len(reference)

    reference_data = pd.DataFrame(data=reference_template, index=np.arange(reference_length))
    # aij_data=aij_data.strip
    reference = np.char.strip(reference)
    citation = np.char.strip(citation)
    reference_data.reference = reference
    reference_data.citation = citation
    return reference_data

def list_aeff_collection_references(atom_rc_file, atom, ion, br=None):
    """
         This function returns a list for all references of recombination coefficients (Aeff)
         for given element and ionic level from the FITS data file ('rc_collection.fits').

     :Returns:
        type=an array of strings. This function returns the references.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_collection.fits')
         atom          : in, required, type=string
                         atom name e.g. 'c'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Keywords:
         br            : in, type=boolean
                         set for the branching ratios (Br), may not necessary

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data-rc')
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_collection.fits')
         >>> atom='c'
         >>> ion='iii' # C III
         >>> list_cii_aeff_reference = atomneb.list_aeff_collection_references(atom_rc_file, atom, ion)
         >>> print(list_cii_aeff_reference)


     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aeff_collection_list(atom_rc_file)
    if (br is not None) == 1:
        prefix = '_br'
    else:
        prefix = '_aeff'
    atom_ion_name = atom.lower() + '_' + ion.lower() + prefix
    aeff_list = np.asarray(element_data_list.aeff_data)
    #aeff_res = [x for x in aeff_list if re.search(atom_ion_name, x)]
    aeff_res = [x for x in aeff_list if atom_ion_name in x]
    ii_length = len(aeff_res)
    if ii_length == 0:
        print('could not find the given element or ion')
        return 0
    nlines1 = ii_length  # temp[0]
    references = ['' for i in range(nlines1)]
    for i in range(0, nlines1):
        select_aeff_data_str1 = aeff_res[i].split('_')
        ref_num = len(select_aeff_data_str1)
        if ref_num > 3:
            references[i] = select_aeff_data_str1[3]
        else:
            references[i] = ''
    return references

def get_aeff_collection_reference_citation(atom_rc_file, atom, ion, br=None, reference=None):
    """
         This function returns the reference citation for a recombination coefficient (Aeff)
         from the 2nd binary table extension of the FITS data file ('rc_collection.fits').

     :Returns:
        type=string. This function returns the Citation.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_collection.fits')
         atom          : in, required, type=string
                         atom name e.g. 'c'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Keywords:
         br            : in, type=boolean
                         set for the branching ratios (Br), may not necessary
         reference     : in, type=string
                         set for the reference, not necessary

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data-rc')
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_collection.fits')
         >>> atom='c'
         >>> ion='iii' # C III
         >>> citation = atomneb.get_aeff_collection_reference_citation(atom_rc_file, atom, ion)
         >>> print(citation)
            Davey, A. R., Storey, P. J. and Kisielius, R., Astron.Astrophys.Suppl., 142, 85, 2000

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_reference = read_aeff_collection_references(atom_rc_file)
    if (br is not None) == 1:
        prefix = '_br'
    else:
        prefix = '_aeff'
    reference_list = np.asarray(element_data_reference.reference)
    if (reference is not None) == 1:
        atom_ion_name = atom.lower() + '_' + ion.lower() + prefix + '_' + reference.upper()
    else:
        atom_ion_name = atom.lower() + '_' + ion.lower() + prefix
    ii = np.where(reference_list == atom_ion_name)[0]  #
    if len(ii) == 0:
        print('could not find the given reference')
        return 0
    citation = np.asarray(element_data_reference.citation[ii])

    return citation

def read_aeff_sh95(atom_rc_file, atom, ion, reference=None, case1=None):
    """
         This function returns the effective recombination coefficients (Aeff) from the table extensions
         of the FITS data file ('rc_SH95.fits').

     :Returns:
        type=an array of data. This function returns the effective recombination coefficients.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_SH95.fits')
         atom          : in, required, type=string
                         atom name e.g. 'h'
         ion           : in, required, type=string
                         ionic level e.g 'ii'

     :Keywords:
         reference     : in, type=string
                         set for the reference,  not necessary
         case1         : in, type=string
                         set for the case 'a' or 'b', defualt 'b'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data-rc')
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_SH95.fits')
         >>> atom='h'
         >>> ion='ii' # H I
         >>> hi_rc_data = atomneb.read_aeff_sh95(atom_rc_file, atom, ion)
         >>> print(hi_rc_data.aeff[0])
            100.00000       500.00000       0.0000000   4.2140000e-27    1.7560000e-27  ...
            ...

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aeff_sh95_list(atom_rc_file)
    if (reference is not None) == 1:
        atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff'
        if (case1 is not None) == 1:
            atom_ion_name = atom_ion_name + '_' + case1
        else:
            atom_ion_name = atom_ion_name + '_b'
        atom_ion_name = atom_ion_name + '_' + reference.upper()
        aeff_list = np.asarray(element_data_list.aeff_data)
        ii = np.where(aeff_list == atom_ion_name)[0]
        if len(ii) == 0:
            print ('could not find the given element or ion')
            return 0
    else:
        atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff'
        if (case1 is not None) == 1:
            atom_ion_name = atom_ion_name + '_' + case1
        else:
            atom_ion_name = atom_ion_name + '_b'
        aeff_list = np.asarray(element_data_list.aeff_data)
        #aeff_res = [x for x in aeff_list if re.search(atom_ion_name, x)]
        aeff_res = [x for x in aeff_list if atom_ion_name in x]
        iii = np.where(aeff_list == aeff_res[0])[0]  #
        if len(iii) == 0:
            print('could not find the given element or ion')
            return 0
        ii = min(iii)
    extension = np.asarray(element_data_list.extension)[ii]
    rc_aeff, header1 = getdata(atom_rc_file, extension, header=True)
    col1 = len(rc_aeff[0])
    row1 = len(rc_aeff)

    aeff_template={'aeff':np.zeros((row1,col1))}
    rc_data = pd.concat([pd.DataFrame(v) for k, v in aeff_template.items()], axis=1, keys=list(aeff_template.keys()))
    rc_data.aeff = rc_aeff

    return rc_data

def read_aeff_sh95_list(atom_rc_file):
    """
         This function returns the list of effective recombination coefficients (Aeff) from the 1st binary table extension
         of the FITS data file ('rc_SH95.fits')

     :Private:

     :Returns:
        type=an array of data. This function returns the aeff_data_list:
              { Aeff_Data:'',
                Extension:0.0}

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_SH95.fits')

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_template = {'aeff_data': '', 'extension': 0 }
    data, hdr = getdata(atom_rc_file, 1, header=True)
    aeff_data=data.AEFF_DATA
    extension = data.EXTENSION
    element_length = len(aeff_data)

    element_data = pd.DataFrame(data=element_template, index=np.arange(element_length))
    #aij_data=aij_data.strip
    aeff_data=np.char.strip(aeff_data)
    element_data.aeff_data=aeff_data
    element_data.extension = extension

    return element_data

def search_aeff_sh95(atom_rc_file, atom, ion):
    """
         This function searches effective recombination coefficients (Aeff) for given element
         and ionic levels in the FITS data file ('rec_SH95.fits'), and returns the data entry.

     :Returns:
        type=array of data. This function returns the Aeff_Data.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_SH95.fits')
         atom          : in, required, type=string
                         atom name e.g. 'h'
         ion           : in, required, type=string
                         ionic level e.g 'ii'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data-rc')
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_SH95.fits')
         >>> atom='h'
         >>> ion='ii' # H I
         >>> list_hi_aeff_data = atomneb.search_aeff_sh95(atom_rc_file, atom, ion)
         >>> print(list_hi_aeff_data)
            h_ii_aeff_a h_ii_aeff_b

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aeff_sh95_list(atom_rc_file)
    atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff'
    aeff_list = np.asarray(element_data_list.aeff_data)
    #aeff_res = [x for x in aeff_list if re.search(atom_ion_name, x)]
    aeff_res = [x for x in aeff_list if atom_ion_name in x]
    if len(aeff_res) == 0:
        print('could not find the given element or ion')
        return 0
    select_aeff_data = aeff_res
    return select_aeff_data

def read_aeff_sh95_references(atom_rc_file):
    """
         This function returns the reference list of recombination coefficients (Aeff) from the 2nd binary table extension
         of the FITS data file ('rc_SH95.fits').

     :Private:

     :Returns:
        type=an array of data. This function returns the aeff_data_reference:
              { Reference:'',
                Citation:''}

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_SH95.fits')

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    reference_template = {'reference': '', 'citation': 0}
    data, hdr = getdata(atom_rc_file, 2, header=True)
    reference = data.ATOMICDATA
    citation = data.REFERENCE
    reference_length = len(reference)

    reference_data = pd.DataFrame(data=reference_template, index=np.arange(reference_length))
    reference = np.char.strip(reference)
    citation = np.char.strip(citation)
    reference_data.reference = reference
    reference_data.citation = citation
    return reference_data

def list_aeff_sh95_references(atom_rc_file, atom, ion):
    """
         This function returns a list for all references of recombination coefficients (Aeff)
         for given element and ionic level from the FITS data file ('rc_SH95.fits').

     :Returns:
        type=an array of strings. This function returns the references.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_SH95.fits')
         atom          : in, required, type=string
                         atom name e.g. 'h'
         ion           : in, required, type=string
                         ionic level e.g 'ii'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = ['atomic-data-rc']
         >>> atom_rc_file = os.path.join(base_dir,data_dir,'rc_SH95.fits')
         >>> atom='h'
         >>> ion='ii' # H I
         >>> list_hi_aeff_reference = atomneb.list_aeff_sh95_references(atom_rc_file, atom, ion)
         >>> print(list_hi_aeff_reference)


     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aeff_sh95_list(atom_rc_file)
    atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff'
    aeff_list = np.asarray(element_data_list.aeff_data)
    #aeff_res = [x for x in aeff_list if re.search(atom_ion_name, x)]
    aeff_res = [x for x in aeff_list if atom_ion_name in x]
    ii_length = len(aeff_res)
    if ii_length == 0:
        print('could not find the given element or ion')
        return 0
    nlines1 = ii_length  # temp[0]
    references = ['' for i in range(nlines1)]
    for i in range(0, nlines1):
        select_aeff_data_str1 = aeff_res[i].split('_')
        ref_num = len(select_aeff_data_str1)
        if ref_num > 4:
            references[i] = select_aeff_data_str1[3]
        else:
            references[i] = ''
    return references

def get_aeff_sh95_reference_citation(atom_rc_file, atom, ion, reference=None, case1=None):
    """
         This function returns the reference citation for a recombination coefficient (Aeff)
         from the 2nd binary table extension of the FITS data file ('rc_SH95.fits').

     :Returns:
        type=string. This function returns the Citation.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_SH95.fits')
         atom          : in, required, type=string
                         atom name e.g. 'h'
         ion           : in, required, type=string
                         ionic level e.g 'ii'

     :Keywords:
         reference     : in, type=string
                         set for the reference,  not necessary
         case1         : in, type=string
                         set for the case 'a' or 'b', defualt 'b'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = ['atomic-data-rc']
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_SH95.fits')
         >>> atom='h'
         >>> ion='ii' # H I
         >>> citation = atomneb.get_aeff_sh95_reference_citation(atom_rc_file, atom, ion)
         >>> print(citation)
            Storey, P. J. and Hummer, D. G., MNRAS, 272, 41S, 1995

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_reference = read_aeff_sh95_references(atom_rc_file)
    atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff'
    reference_list = np.asarray(element_data_reference.reference)
    if (case1 is not None) == 1:
        atom_ion_name = atom_ion_name + '_' + case1
    else:
        atom_ion_name = atom_ion_name + '_b'
    if (reference is not None) == 1:
        atom_ion_name = atom_ion_name + '_' + reference.upper()
    ii = np.where(reference_list == atom_ion_name)[0]  #
    if len(ii) == 0:
        print('could not find the given reference')
        return 0
    citation = np.asarray(element_data_reference.citation[ii])

    return citation

def read_aeff_ppb91(atom_rc_file, atom, ion, reference=None):
    """
         This function returns the effective recombination coefficients (Aeff) from the table extensions
         of the FITS data file ('rc_PPB91.fits').

     :Returns:
        type=an array of data. This function returns the effective recombination coefficients:
              { Ion: ' '
                Case1:''
                Wavelength:0.0,
                a: 0.0,
                b: 0.0,
                c: 0.0,
                d: 0.0,
                br: 0.0,
                y: 0.0}

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_PPB91.fits')
         atom          : in, required, type=string
                         atom name e.g. 'c'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Keywords:
         reference     : in, type=string
                         set for the reference,  not necessary

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data-rc')
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_PPB91.fits')
         >>> atom='c'
         >>> ion='iii' # C II
         >>> cii_rc_data = atomneb.read_aeff_ppb91(atom_rc_file, atom, ion)
         >>> n_line = len(cii_rc_data.wavelength)
         >>> for i in range(0, n_line):
         >>>           print(cii_rc_data.ion[i], cii_rc_data.case1[i], cii_rc_data.wavelength[i],
         >>>                 cii_rc_data.a[i], cii_rc_data.b[i], cii_rc_data.c[i],
         >>>                 cii_rc_data.d[i], cii_rc_data.br[i], cii_rc_data.q[i], cii_rc_data.y[i])
            C2+A       9903.4600      0.69700000     -0.78400000      ...
            C2+A       4267.1500       1.0110000     -0.75400000      ...
            ...

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aeff_ppb91_list(atom_rc_file)
    if (reference is not None) == 1:
        atom_ion_name = atom.lower() + '_' + ion.lower()+ '_aeff_' + reference.upper()
        aeff_list = np.asarray(element_data_list.aeff_data)
        ii = np.where(aeff_list == atom_ion_name)[0]
        if len(ii) == 0:
            print ('could not find the given element or ion')
            return 0
    else:
        atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff'
        aeff_list = np.asarray(element_data_list.aeff_data)
        #aeff_res = [x for x in aeff_list if re.search(atom_ion_name, x)]
        aeff_res = [x for x in aeff_list if atom_ion_name in x]
        iii = np.where(aeff_list == aeff_res[0])[0]  #
        if len(iii) == 0:
            print('could not find the given element or ion')
            return 0
        ii = min(iii)

    extension = np.asarray(element_data_list.extension)[ii]

    rc_template={'ion':'', 'case1': '', 'wavelength': float(0.0),
                 'a': float(0.0), 'b':float(0.0), 'c': float(0.0),
                 'd':float(0.0), 'br':float(0.0), 'q':'', 'y':float(0.0)}
    data, header1 = getdata(atom_rc_file, int(extension), header=True)
    ion = data.ION
    case1 = data.CASE1
    wavelength = data.WAVELENGTH
    a = data.A
    b = data.B
    c = data.C
    d = data.D
    br = data.BR
    q = data.Q
    y = data.Y
    ion = np.char.strip(ion)
    case1 = np.char.strip(case1)
    n_line = len(wavelength)
    rc_data = pd.DataFrame(data=rc_template, index=np.arange(n_line))
    rc_data.ion = ion
    rc_data.case1 = case1
    rc_data.wavelength = wavelength
    rc_data.a = a
    rc_data.b = b
    rc_data.c = c
    rc_data.d = d
    rc_data.br = br
    rc_data.y = y
    return rc_data

def read_aeff_ppb91_list(atom_rc_file):
    """
         This function returns the list of effective recombination coefficients (Aeff) from the 1st binary table extension
         of the FITS data file ('rc_PPB91.fits')

     :Private:

     :Returns:
        type=an array of data. This function returns the aeff_data_list:
              { Aeff_Data:'',
                Extension:0.0}

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_PPB91.fits')

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_template = {'aeff_data': '', 'extension': 0}
    data, hdr = getdata(atom_rc_file, 1, header=True)
    aeff_data=data.AEFF_DATA
    extension = data.EXTENSION
    element_length = len(aeff_data)

    element_data = pd.DataFrame(data=element_template, index=np.arange(element_length))
    aeff_data=np.char.strip(aeff_data)
    element_data.aeff_data=aeff_data
    element_data.extension = extension
    return element_data

def search_aeff_ppb91(atom_rc_file, atom, ion):
    """
         This function searches effective recombination coefficients (Aeff) for given element
         and ionic levels in the FITS data file ('rec_PPB91.fits'), and returns the data entry.

     :Returns:
        type=array of data. This function returns the Aeff_Data.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_PPB91.fits')
         atom          : in, required, type=string
                         atom name e.g. 'c'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data-rc')
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_PPB91.fits')
         >>> atom='c'
         >>> ion='iii'
         >>> list_cii_aeff_data = atomneb.search_aeff_ppb91(atom_rc_file, atom, ion)
         >>> print(list_cii_aeff_data)
            c_iii_aeff

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aeff_ppb91_list(atom_rc_file)
    atom_ion_name = atom.lower() + '_' + ion.lower()  + '_aeff'
    aeff_list = np.asarray(element_data_list.aeff_data)
    #aeff_res = [x for x in aeff_list if re.search(atom_ion_name, x)]
    aeff_res = [x for x in aeff_list if atom_ion_name in x]
    if len(aeff_res) == 0:
        print('could not find the given element or ion')
        return 0
    select_aeff_data = aeff_res
    return select_aeff_data

def read_aeff_ppb91_references(atom_rc_file):
    """
         This function returns the reference list of recombination coefficients (Aeff) from the 2nd binary table extension
         of the FITS data file ('rc_PPB91.fits').

     :Private:

     :Returns:
        type=an array of data. This function returns the aeff_data_reference:
              { Reference:'',
                Citation:''}

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_PPB91.fits')

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    reference_template = {'reference': '', 'citation': 0}
    data, hdr = getdata(atom_rc_file, 2, header=True)
    reference = data.ATOMICDATA
    citation = data.REFERENCE
    reference_length = len(reference)

    reference_data = pd.DataFrame(data=reference_template, index=np.arange(reference_length))
    reference = np.char.strip(reference)
    citation = np.char.strip(citation)
    reference_data.reference = reference
    reference_data.citation = citation
    return reference_data

def list_aeff_ppb91_references(atom_rc_file, atom, ion):
    """
         This function returns a list for all references of recombination coefficients (Aeff)
         for given element and ionic level from the FITS data file ('rc_PPB91.fits').

     :Returns:
        type=an array of strings. This function returns the references.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_PPB91.fits')
         atom          : in, required, type=string
                         atom name e.g. 'c'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = ['atomic-data-rc']
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_PPB91.fits')
         >>> atom='c'
         >>> ion='iii'
         >>> list_cii_aeff_reference = atomneb.list_aeff_ppb91_references(atom_rc_file, atom, ion)
         >>> print(list_cii_aeff_reference)


     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aeff_ppb91_list(atom_rc_file)
    atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff'
    aeff_list = np.asarray(element_data_list.aeff_data)
    #aeff_res = [x for x in aeff_list if re.search(atom_ion_name, x)]
    aeff_res = [x for x in aeff_list if atom_ion_name in x]
    ii_length = len(aeff_res)
    if ii_length == 0:
        print('could not find the given element or ion')
        return 0
    nlines1 = ii_length  # temp[0]
    references = ['' for i in range(nlines1)]
    for i in range(0, nlines1):
        select_aeff_data_str1 = aeff_res[i].split('_')
        ref_num = len(select_aeff_data_str1)
        if ref_num > 3:
            references[i] = select_aeff_data_str1[3]
        else:
            references[i] = ''
    return references

def get_aeff_ppb91_reference_citation(atom_rc_file, atom, ion, reference=None):
    """
        This function returns the reference citation for a recombination coefficient (Aeff)
        from the 2nd binary table extension of the FITS data file ('rc_PPB91.fits').

    :Returns:
       type=string. This function returns the Citation.

    :Params:
        Atom_RC_file  : in, required, type=string
                        the FITS data file name ('rc_PPB91.fits')
        atom          : in, required, type=string
                        atom name e.g. 'c'
        ion           : in, required, type=string
                        ionic level e.g 'iii'

    :Keywords:
        reference     : in, type=string
                        set for the reference,  not necessary

    :Examples:
       For example::

        >>> import atomneb
        >>> import numpy as np
        >>> import os
        >>> base_dir = os.getcwd()
        >>> data_dir = ['atomic-data-rc']
        >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_PPB91.fits')
        >>> atom='c'
        >>> ion='iii'
        >>> citation = atomneb.get_aeff_ppb91_reference_citation(atom_rc_file, atom, ion)
        >>> print(citation)
           Pequignot, D., Petitjean, P. and Boisson, C. Astron.Astrophys., 251, 680, 1991

    :Categories:
      Recombination Lines

    :Dirs:
     ./
         Main routines

    :Author:
      Ashkbiz Danehkar

    :Copyright:
      This library is released under a GNU General Public License.

    :Version:
      0.2.0

    :History:
        15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """

    element_data_reference = read_aeff_ppb91_references(atom_rc_file)
    reference_list = np.asarray(element_data_reference.reference)
    if (reference is not None) == 1:
        atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff_' + '_' + reference.upper()
    else:
        atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff'
    ii = np.where(reference_list == atom_ion_name)[0]  #
    if len(ii) == 0:
        print('could not find the given reference')
        return 0
    citation = np.asarray(element_data_reference.citation[ii])

    return citation

def read_aeff_he_i_pfsd12(atom_rc_file, atom, ion, wavelength=None, reference=None):
    """
         This function returns the effective recombination coefficients (Aeff) from the table extensions
         of the FITS data file ('rc_he_ii_PFSD12.fits').

     :Returns:
        type=an array of data. This function returns the effective recombination coefficients.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_he_ii_PFSD12.fits')
         atom          : in, required, type=string
                         atom name e.g. 'he'
         ion           : in, required, type=string
                         ionic level e.g 'ii'

     :Keywords:
         wavelength    : in, type=boolean
                         set for returning the wavelengths
         reference     : in, type=string
                         set for the reference, not necessary

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data-rc')
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_he_ii_PFSD12.fits')
         >>> atom='he'
         >>> ion='ii' # He I
         >>> hei_rc_data = atomneb.read_aeff_he_i_pfsd12(atom_rc_file, atom, ion)
         >>> hei_rc_data_wave = atomneb.read_aeff_he_i_pfsd12(atom_rc_file, atom, ion, wavelength=True)
         >>> print(hei_rc_data.aeff[0])
            5000.0000       10.000000      -25.379540      -25.058970      -25.948440      ...
         >>> n_line = len(hei_rc_data_wave.wavelength)
         >>> for i in range(0, n_line):
         >>>     print(hei_rc_data_wave.wavelength[i],
         >>>           hei_rc_data_wave.lowerterm[i], hei_rc_data_wave.upperterm[i])
            2945.00005p^{3}P2s^{3}S
            3188.00004p^{3}P2s^{3}S
            3614.00005p^{1}P2s^{1}S
            ...

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aeff_he_i_pfsd12_list(atom_rc_file)
    if (wavelength is not None) == 1:
        prefix = '_wavelength'
    else:
        prefix = '_aeff'
    if (reference is not None) == 1:
        atom_ion_name = atom.lower() + '_' + ion.lower() + prefix + '_' + reference.upper()
        aeff_list = np.asarray(element_data_list.aeff_data)
        ii = np.where(aeff_list == atom_ion_name)[0]  #
        if len(ii) == 0:
            print('could not find the given reference')
            return
    else:
        reference = 'PFSD13'
        if (wavelength is not None) == 1:
            atom_ion_name = atom.lower() + '_' + ion.lower() + prefix
        else:
            atom_ion_name = atom.lower() + '_' + ion.lower() + prefix + '_' + reference.upper()
        aeff_list=np.asarray(element_data_list.aeff_data)
        #aeff_res = [x for x in aeff_list if re.search(atom_ion_name, x)]
        aeff_res = [x for x in aeff_list if atom_ion_name in x]
        iii = np.where(aeff_list == aeff_res[0])[0]  #
        if len(iii) == 0:
            print('could not find the given element or ion')
            return 0
        ii = min(iii)
    extension = np.asarray(element_data_list.extension)[ii]

    atom_ion_name = atom.lower() + '_' + ion.lower()+ prefix
    if atom_ion_name == 'he_ii_aeff':
        rc_aeff, header1 = getdata(atom_rc_file, extension, header=True)
        col1 = len(rc_aeff[0])
        row1 = len(rc_aeff)
        aeff_template = {'aeff': np.zeros((row1, col1))}
        rc_data = pd.concat([pd.DataFrame(v) for k, v in aeff_template.items()], axis=1,
                            keys=list(aeff_template.keys()))
        rc_data.aeff = rc_aeff
    elif atom_ion_name == 'he_ii_wavelength':
        rc_template={'wavelength': float(0.0), 'lowerterm':'', 'upperterm':''}
        data, header1 = getdata(atom_rc_file, extension, header=True)
        wavelength = data.WAVELENGTH
        lowerterm = data.LOWERTERM
        upperterm = data.UPPERTERM
        lowerterm = np.char.strip(lowerterm)
        upperterm = np.char.strip(upperterm)
        n_line = len(wavelength)
        rc_data = pd.DataFrame(data=rc_template, index=np.arange(n_line))
        rc_data.wavelength = wavelength
        rc_data.lowerterm = lowerterm
        rc_data.upperterm = upperterm
    else:
        raise RuntimeError('no match found for expression')
    return rc_data

def read_aeff_he_i_pfsd12_list(atom_rc_file):
    """
         This function returns the list of effective recombination coefficients (Aeff) from the 1st binary table extension
         of the FITS data file ('rc_he_ii_PFSD12.fits')

     :Private:

     :Returns:
        type=an array of data. This function returns the aeff_data_list:
               { Aeff_Data:'',
                 Extension:0.0}

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_he_ii_PFSD12.fits')

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_template = {'aeff_data': '', 'extension': 0 }
    data, hdr = getdata(atom_rc_file, 1, header=True)
    aeff_data=data.AEFF_DATA
    extension = data.EXTENSION
    element_length = len(aeff_data)

    element_data = pd.DataFrame(data=element_template, index=np.arange(element_length))
    aeff_data=np.char.strip(aeff_data)
    element_data.aeff_data=aeff_data
    element_data.extension = extension

    return element_data

def search_aeff_he_i_pfsd12(atom_rc_file, atom, ion):
    """
         This function searches effective recombination coefficients (Aeff) for given element
         and ionic levels in the FITS data file ('rec_he_ii_PFSD12.fits'), and returns the data entry.

     :Returns:
        type=array of data. This function returns the Aeff_Data.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_he_ii_PFSD12.fits')
         atom          : in, required, type=string
                         atom name e.g. 'he'
         ion           : in, required, type=string
                         ionic level e.g 'ii'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data-rc')
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_he_ii_PFSD12.fits')
         >>> atom='he'
         >>> ion='ii' # He I
         >>> list_hei_aeff_data = atomneb.search_aeff_he_i_pfsd12(atom_rc_file, atom, ion)
         >>> print(list_hei_aeff_data)
            he_ii_aeff_PFSD12 he_ii_aeff_PFSD13

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aeff_he_i_pfsd12_list(atom_rc_file)
    atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff'
    aeff_list = np.asarray(element_data_list.aeff_data)
    #aeff_res = [x for x in aeff_list if re.search(atom_ion_name, x)]
    aeff_res = [x for x in aeff_list if atom_ion_name in x]
    if len(aeff_res) == 0:
        print('could not find the given element or ion')
        return 0
    select_aeff_data = aeff_res
    return select_aeff_data

def read_aeff_he_i_pfsd12_references(atom_rc_file):
    """
         This function returns the reference list of recombination coefficients (Aeff) from the 2nd binary table extension
         of the FITS data file ('rc_he_ii_PFSD12.fits').

     :Private:

     :Returns:
        type=an array of data. This function returns the aeff_data_reference:
              { Reference:'',
                Citation:''}

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_he_ii_PFSD12.fits')

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    reference_template = {'reference': '', 'citation': 0}
    data, hdr = getdata(atom_rc_file, 2, header=True)
    reference = data.ATOMICDATA
    citation = data.REFERENCE
    reference_length = len(reference)

    reference_data = pd.DataFrame(data=reference_template, index=np.arange(reference_length))
    reference = np.char.strip(reference)
    citation = np.char.strip(citation)
    reference_data.reference = reference
    reference_data.citation = citation
    return reference_data

def list_aeff_he_i_pfsd12_references(atom_rc_file, atom, ion):
    """
         This function returns a list for all references of recombination coefficients (Aeff)
         for given element and ionic level from the FITS data file ('rc_he_ii_PFSD12.fits').

     :Returns:
        type=an array of strings. This function returns the references.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_he_ii_PFSD12.fits')
         atom          : in, required, type=string
                         atom name e.g. 'he'
         ion           : in, required, type=string
                         ionic level e.g 'ii'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = ['atomic-data-rc']
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_he_ii_PFSD12.fits')
         >>> atom='he'
         >>> ion='ii' # He I
         >>> list_hei_aeff_reference = atomneb.list_aeff_he_i_pfsd12_references(atom_rc_file, atom, ion)
         >>> print(list_hei_aeff_reference)
            PFSD12 PFSD13

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aeff_he_i_pfsd12_list(atom_rc_file)
    atom_ion_name = atom.lower() + '_' + ion.lower()+ '_aeff'
    aeff_list = np.asarray(element_data_list.aeff_data)
    #aeff_res = [x for x in aeff_list if re.search(atom_ion_name, x)]
    aeff_res = [x for x in aeff_list if atom_ion_name in x]
    ii_length = len(aeff_res)
    if ii_length == 0:
        print('could not find the given element or ion')
        return 0
    nlines1 = ii_length  # temp[0]
    references = ['' for i in range(nlines1)]
    for i in range(0, nlines1):
        select_aeff_data_str1 = aeff_res[i].split('_')
        ref_num = len(select_aeff_data_str1)
        if ref_num > 3:
            references[i] = select_aeff_data_str1[3]
        else:
            references[i] = ''
    return references

def get_aeff_he_i_pfsd12_reference_citation(atom_rc_file, atom, ion, reference=None):
    """
         This function returns the reference citation for a recombination coefficient (Aeff)
         from the 2nd binary table extension of the FITS data file ('rc_he_ii_PFSD12.fits').

     :Returns:
        type=string. This function returns the Citation.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_he_ii_PFSD12.fits')
         atom          : in, required, type=string
                         atom name e.g. 'he'
         ion           : in, required, type=string
                         ionic level e.g 'ii'

     :Keywords:
         reference     : in, type=string
                         set for the reference e.g. 'PFSD13',  may not necessary

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = ['atomic-data-rc']
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_he_ii_PFSD12.fits')
         >>> atom='he'
         >>> ion='ii' # He I
         >>> reference='PFSD13'
         >>> citation = atomneb.get_aeff_he_i_pfsd12_reference_citation(atom_rc_file, atom, ion, reference=reference)
         >>> print(citation)
            Porter, R. L., Ferland, G. J., Storey, P. J. and Detisch, M. J., MNRAS, 433L, 89, 2013

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         15/01/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """

    #  reference_template={Reference:'', Citation: ''}
    element_data_reference = read_aeff_he_i_pfsd12_references(atom_rc_file)
    reference_list = np.asarray(element_data_reference.reference)
    if (reference is not None) == 1:
        atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff_' + reference.upper()
    else:
        atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff'
    ii = np.where(reference_list == atom_ion_name)[0]  #
    if len(ii) == 0:
        print('could not find the given reference')
        return 0
    citation = np.asarray(element_data_reference.citation[ii])

    return citation

def read_aeff_n_ii_fsl13(atom_rc_file, atom, ion, wavelength_range, wavelength=None, reference=None):
    """
         This function returns the effective recombination coefficients (Aeff) from the table extensions
         of the FITS data file ('rc_n_iii_FSL13.fits').

     :Returns:
        type=an array of data. This function returns the effective recombination coefficients.

     :Params:
         Atom_RC_file     : in, required, type=string
                            the FITS data file name ('rc_n_iii_FSL13.fits')
         atom             : in, required, type=string
                            atom name e.g. 'n'
         ion              : in, required, type=string
                            ionic level e.g 'iii'
         wavelength_range : in, required, type=array
                            wavelength range e.g. [4400.0, 7100.0]

     :Keywords:
         wavelength    : in, type=boolean
                         set for returning the wavelengths
         reference     : in, type=string
                         set for the reference, not necessary

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data-rc')
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_n_iii_FSL13.fits')
         >>> atom='n'
         >>> ion='iii' # N II
         >>> wavelength_range=[4400.0, 7100.0]
         >>> nii_rc_data = atomneb.read_aeff_n_ii_fsl13(atom_rc_file, atom, ion, wavelength_range)
         >>> nii_rc_data_wave = atomneb.read_aeff_n_ii_fsl13(atom_rc_file, atom, ion, wavelength_range, wavelength=True)
         >>> print(nii_rc_data.aeff[0])
            255.000      79.5000      47.3000      12.5000      ...
         >>> n_line = len(nii_rc_data_wave.wavelength)
         >>> for i in range(0, n_line):
         >>>     print(nii_rc_data_wave.wavelength[i], nii_rc_data_wave.tr[i], nii_rc_data_wave.trans[i])
            6413.236g - 4f2p6g G[9/2]o4 - 2p4f F[7/2]e3
            6556.326g - 4f2p6g G[9/2]o5 - 2p4f G[7/2]e4
            6456.976g - 4f2p6g G[9/2]o5 - 2p4f F[7/2]e4
            ...

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         03/07/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aeff_n_ii_fsl13_list(atom_rc_file)
    wave_min = min(wavelength_range)
    wave_max = max(wavelength_range)
    wave_loc1 = np.where((element_data_list.wavelength >= wave_min)&(element_data_list.wavelength <= wave_max))[0]
    wave_size = len(wave_loc1)
    if wave_size == 0:
        return 0
    rc_element_template={'wavelength': float(0.0), 'aeff':[np.zeros((7,4))]}
    rc_data = pd.DataFrame(data=rc_element_template, index=np.arange(wave_size))

    extension1 = np.asarray(element_data_list.extension)[wave_loc1]
    wavelength1 = np.asarray(element_data_list.wavelength)[wave_loc1]

    if (wavelength is not None) == 1:
        rc_wave_template={'wavelength': float(0.0), 'tr':'',  'trans': '', 't_x': ''}
        rc_wave = pd.DataFrame(data=rc_wave_template, index=np.arange(wave_size))
        rc_wave.wavelength = np.asarray(element_data_list.wavelength)[wave_loc1]
        rc_wave.tr =  np.asarray(element_data_list.tr)[wave_loc1]
        rc_wave.trans =  np.asarray(element_data_list.trans)[wave_loc1]
        rc_wave.t_x =  np.asarray(element_data_list.t_x)[wave_loc1]
        return rc_wave
    for i in range(0, wave_size ):
        rc_aeff, header1 = getdata(atom_rc_file, int(extension1[i]), header=True)
        # temp=size(rc_aeff,/DIMENSIONS)
        # col1=temp[0]
        # row1=temp[1]
        # aeff_template={Aeff:dblarr(col1,row1)}
        # rc_data=replicate(aeff_template, 1)
        #rc_data.wavelength[i] = wavelength1[i]
        rc_data.loc[i, 'wavelength']=wavelength1[i]
        #rc_data.aeff[i] = rc_aeff
        rc_data.loc[i, 'aeff']=[rc_aeff]

    return rc_data

def read_aeff_n_ii_fsl13_list(atom_rc_file):
    """
         This function returns the list of effective recombination coefficients (Aeff) from the 1st binary table extension
         of the FITS data file ('rc_n_iii_FSL13.fits')

     :Private:

     :Returns:
        type=an array of data. This function returns the aeff_data_list:
              {Aeff_Data:'',  Extension:0, $
               IND:long(0), Wavelength: float(0.0), $
               Tr:'',  Trans: '', T_X: ''}

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_n_iii_FSL13.fits')

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         03/07/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_template={'aeff_data':'',  'extension':0,
                      'ind':int(0), 'wavelength': float(0.0),
                      'tr':'',  'trans': '', 't_x': ''}
    data, header1 = getdata(atom_rc_file, 1, header=True)
    aeff_data=data.AEFF_DATA
    extension=data.EXTENSION
    ind = data.IND
    wavelength = data.WAVELENGTH
    tr = data.TR
    trans = data.TRANS

    element_length = len(aeff_data)
    element_data = pd.DataFrame(data=element_template, index=np.arange(element_length))

    aeff_data=np.char.strip(aeff_data)
    tr = np.char.strip(tr)
    trans = np.char.strip(trans)

    element_data.aeff_data = aeff_data
    element_data.extension = extension
    element_data.ind = ind
    element_data.wavelength = wavelength
    element_data.tr = tr
    element_data.trans = trans
    return element_data

def search_aeff_n_ii_fsl13(atom_rc_file, atom, ion, wavelength):
    """
         This function searches effective recombination coefficients (Aeff) for given element
         and ionic levels in the FITS data file ('rc_n_iii_FSL13.fits'), and returns the data entry.

     :Returns:
        type=array of data. This function returns the Aeff_Data.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_n_iii_FSL13.fits')
         atom          : in, required, type=string
                         atom name e.g. 'n'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Keywords:
         wavelength    : in, type=float
                         set the wavelengths
     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data-rc')
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_n_iii_FSL13.fits')
         >>> atom='n'
         >>> ion='iii' # N II
         >>> wavelength=5679.56
         >>> list_nii_aeff_data = atomneb.search_aeff_n_ii_fsl13(atom_rc_file, atom, ion, wavelength)
         >>> print(np.asarray(list_nii_aeff_data.wavelength))
            5679.56
         >>> print(np.asarray(list_nii_aeff_data.aeff))
            7810.00      1780.00      850.000      151.000      74.4000      53.1000      47.4000
            7370.00      1700.00      886.000      206.000      110.000      80.1000      70.8000
            7730.00      1680.00      900.000      239.000      138.000      103.000      92.9000
            8520.00      1710.00      905.000      244.000      142.000      107.000      97.0000

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         03/07/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aeff_n_ii_fsl13_list(atom_rc_file)

    ii = np.where((abs(element_data_list.wavelength - wavelength) < 0.025))[0]
    ii_length = len(ii)
    if len(ii) == 0:
        print('could not find the given wavelength')
        return 0
    if ii_length > 1:
        value_min = min(abs(element_data_list[ii].wavelength - wavelength))
        ii = np.where((abs(element_data_list.wavelength - wavelength) == value_min))[0]
        ii_length = 1
    extension1 = np.asarray(element_data_list.extension)[ii]
    wavelength1 = np.asarray(element_data_list.wavelength)[ii]

    rc_element_template={'wavelength': float(0.0), 'aeff':[np.zeros((7,4))]}
    select_aeff_data = pd.DataFrame(data=rc_element_template, index=np.arange(ii_length))
    for i in range(0, ii_length):
        rc_aeff, header1 = getdata(atom_rc_file, extension1[i], header=True)
        #select_aeff_data.wavelength[i] = wavelength1[i]
        select_aeff_data.loc[i, 'wavelength']=wavelength1[i]
        # select_aeff_data.aeff[i] = rc_aeff
        select_aeff_data.loc[i, 'aeff'] = [rc_aeff]
    return select_aeff_data

def read_aeff_n_ii_fsl13_references(atom_rc_file):
    """
         This function returns the reference list of recombination coefficients (Aeff) from the 2nd binary table extension
         of the FITS data file ('rc_n_iii_FSL13.fits').

     :Private:

     :Returns:
        type=an array of data. This function returns the aeff_data_reference:
              { Reference:'',
                Citation:''}

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_n_iii_FSL13.fits')

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         03/07/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    reference_template = {'reference': '', 'citation': 0}
    data, hdr = getdata(atom_rc_file, 2, header=True)
    reference = data.ATOMICDATA
    citation = data.REFERENCE
    reference_length = len(reference)

    reference_data = pd.DataFrame(data=reference_template, index=np.arange(reference_length))
    reference = np.char.strip(reference)
    citation = np.char.strip(citation)
    reference_data.reference = reference
    reference_data.citation = citation
    return reference_data

def list_aeff_n_ii_fsl13_references(atom_rc_file, atom, ion):
    """
         This function returns a list for all references of recombination coefficients (Aeff)
         for given element and ionic level from the FITS data file ('rc_n_iii_FSL13.fits').

     :Returns:
        type=an array of strings. This function returns the references.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_n_iii_FSL13.fits')
         atom          : in, required, type=string
                         atom name e.g. 'n'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = ['atomic-data-rc']
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_n_iii_FSL13.fits')
         >>> atom='n'
         >>> ion='iii' # N II
         >>> list_nii_aeff_reference = atomneb.list_aeff_n_ii_fsl13_references(atom_rc_file, atom, ion)
         >>> print(list_nii_aeff_reference)

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         03/07/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    rc_reference_data = read_aeff_n_ii_fsl13_references(atom_rc_file)

    atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff'
    rc_list = np.asarray(rc_reference_data.reference)
    #rc_res = [x for x in rc_list if re.search(atom_ion_name, x)]
    rc_res = [x for x in rc_list if atom_ion_name in x]
    ii_length = len(rc_res)
    if ii_length == 0:
        print('could not find the given element or ion')
        return 0
    nlines1 = ii_length  # temp[0]
    references = ['' for i in range(nlines1)]
    for i in range(0, nlines1):
        select_rc_data_str1 = rc_res[i].split('_')
        ref_num = len(select_rc_data_str1)
        if ref_num > 3:
            references[i] = select_rc_data_str1[3]
        else:
            references[i] = ''
    return references

def get_aeff_n_ii_fsl13_reference_citation(atom_rc_file, atom, ion, reference=None):
    """
         This function returns the reference citation for a recombination coefficient (Aeff)
         from the 2nd binary table extension of the FITS data file ('rc_n_iii_FSL13.fits').

     :Returns:
        type=string. This function returns the Citation.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_n_iii_FSL13.fits')
         atom          : in, required, type=string
                         atom name e.g. 'n'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Keywords:
         reference     : in, type=string
                         set for the reference e.g. 'FSL13',  may not necessary

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = ['atomic-data-rc']
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_n_iii_FSL13.fits')
         >>> atom='n'
         >>> ion='iii' # N II
         >>> reference='FSL13'
         >>> citation = atomneb.get_aeff_n_ii_fsl13_reference_citation(atom_rc_file, atom, ion)
         >>> print(citation)
            Fang X., Storey P.J., and Liu X.-W., R. 2011, Astron.Astrophys. 530, A18; 2013, Astron.Astrophys. 550, C2

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         03/07/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_reference = read_aeff_n_ii_fsl13_references(atom_rc_file)
    reference_list = np.asarray(element_data_reference.reference)
    atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff'
    if (reference is not None) == 1:
        atom_ion_name = atom_ion_name + '_' + reference.upper()

    ii = np.where(reference_list == atom_ion_name)[0]  #
    if len(ii) == 0:
        print('could not find the given reference')
        return 0
    citation = np.asarray(element_data_reference.citation[ii])

    return citation

def read_aeff_o_ii_ssb17(atom_rc_file, atom, ion, case1, wavelength_range, wavelength=None, reference=None):
    """
         This function returns the effective recombination coefficients (Aeff) from the table extensions
         of the FITS data file ('rc_o_iii_SSB17.fits').

     :Returns:
        type=an array of data. This function returns the effective recombination coefficients.

     :Params:
         Atom_RC_file     : in, required, type=string
                            the FITS data file name ('rc_o_iii_SSB17.fits')
         atom             : in, required, type=string
                            atom name e.g. 'o'
         ion              : in, required, type=string
                            ionic level e.g 'iii'
         case1            : in, type=string
                            set for the case 'a' or 'b', defualt 'b'
         wavelength_range : in, required, type=array
                            wavelength range e.g. [5320.0, 5330.0]

     :Keywords:
         wavelength    : in, type=boolean
                         set for returning the wavelengths
         reference     : in, type=string
                         set for the reference, not necessary

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data-rc')
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_o_iii_SSB17_orl_case_b.fits')
         >>> atom='o'
         >>> ion='iii' # O II
         >>> case1='B'
         >>> wavelength_range=[5320.0, 5330.0]
         >>> oii_rc_data = atomneb.read_aeff_o_ii_ssb17(atom_rc_file, atom, ion, case1, wavelength_range)
         >>> oii_rc_data_wave = atomneb.read_aeff_o_ii_ssb17(atom_rc_file, atom, ion,
         >>>                                                 case1, wavelength_range, wavelength=True)
         >>> print(oii_rc_data.aeff[0])
            1.64100e-30  1.60000e-30  1.56400e-30  1.54100e-30 ...
         >>> n_line = len(oii_rc_data_wave.wavelength)
         >>> for i in range(0, n_line):
         >>>      print(oii_rc_data_wave.wavelength[i], oii_rc_data_wave.lower_term[i], oii_rc_data_wave.upper_term[i])
            5327.172s22p2(1S)3p 2Po
            5325.422s22p2(1S)3p 2Po
            5327.182s22p2(1D)3d 2Ge
            ...

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         03/07/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aeff_o_ii_ssb17_list(atom_rc_file)
    wave_min = min(wavelength_range)
    wave_max = max(wavelength_range)
    wave_loc1 = np.where((element_data_list.wavelength >= wave_min)&(element_data_list.wavelength <= wave_max))[0]
    wave_size = len(wave_loc1)
    if wave_size == 0:
        return 0
    rc_element_template={'wavelength': float(0.0), 'aeff':[np.zeros((7,4))]}
    rc_data = pd.DataFrame(data=rc_element_template, index=np.arange(wave_size))

    extension1 = np.asarray(element_data_list.extension)[wave_loc1]
    wavelength1 = np.asarray(element_data_list.wavelength)[wave_loc1]

    if (wavelength is not None) == 1:
        rc_wave_template={'wavelength': float(0.0), 'lower_term':'',  'upper_term': ''}
        rc_wave = pd.DataFrame(data=rc_wave_template, index=np.arange(wave_size))
        rc_wave.wavelength = np.asarray(element_data_list.wavelength)[wave_loc1]
        rc_wave.lower_term = np.asarray(element_data_list.lower_term)[wave_loc1]
        rc_wave.upper_term = np.asarray(element_data_list.upper_term)[wave_loc1]
        return rc_wave
    for i in range(0, wave_size ):
        rc_aeff, header1 = getdata(atom_rc_file, extension1[i], header=True)
        rc_data.loc[i, 'wavelength']=wavelength1[i]
        rc_data.loc[i, 'aeff']=[rc_aeff]

    return rc_data

def read_aeff_o_ii_ssb17_list(atom_rc_file):
    """
         This function returns the list of effective recombination coefficients (Aeff) from the 1st binary table extension
         of the FITS data file ('rc_o_iii_SSB17.fits')

     :Private:

     :Returns:
        type=an array of data. This function returns the aeff_data_list:
              {Aeff_Data:'',  Extension:0, $
               IND:long(0), Wavelength: float(0.0), $
               Case1:'',  lower_term: '', upper_term: ''}

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_o_iii_SSB17.fits')

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         03/07/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_template={'aeff_data':'',  'extension':0,
                      'ind':int(0), 'wavelength': float(0.0),
                      'case1': '', 'lower_term':'', 'upper_term':''}
    data, header1 = getdata(atom_rc_file, 1, header=True)
    aeff_data=data.AEFF_DATA
    extension=data.EXTENSION
    ind = data.IND
    wavelength = data.WAVELENGTH
    case1 = data.CASE1
    lower_term = data.LOWER_TERM
    upper_term = data.UPPER_TERM

    element_length = len(aeff_data)
    element_data = pd.DataFrame(data=element_template, index=np.arange(element_length))

    aeff_data=np.char.strip(aeff_data)
    case1 = np.char.strip(case1)
    lower_term = np.char.strip(lower_term)
    upper_term = np.char.strip(upper_term)

    element_data.aeff_data = aeff_data
    element_data.extension = extension
    element_data.ind = ind
    element_data.wavelength = wavelength
    element_data.case1 = case1
    element_data.lower_term = lower_term
    element_data.upper_term = upper_term

    return element_data

def search_aeff_o_ii_ssb17(atom_rc_file, atom, ion, case1, wavelength):
    """
         This function searches effective recombination coefficients (Aeff) for given element
         and ionic levels in the FITS data file ('rc_o_iii_SSB17.fits'), and returns the data entry.

     :Returns:
        type=array of data. This function returns the Aeff_Data.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_o_iii_SSB17.fits')
         atom          : in, required, type=string
                         atom name e.g. 'o'
         ion           : in, required, type=string
                         ionic level e.g 'iii'
         case1            : in, type=string
                            set for the case 'a' or 'b', defualt 'b'
         wavelength    : in, type=float
                         set the wavelengths

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = os.path.join('..','atomic-data-rc')
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_o_iii_SSB17_orl_case_b.fits')
         >>> atom='o'
         >>> ion='iii' # O II
         >>> case1='B'
         >>> wavelength=5325.42
         >>> list_oii_aeff_data = atomneb.search_aeff_o_ii_ssb17(atom_rc_file, atom, ion, case1, wavelength)
         >>> print(np.asarray(list_oii_aeff_data.wavelength))
            5325.42
         >>> print(np.asarray(list_oii_aeff_data.aeff))
            3.41800e-32  3.33300e-32  3.25700e-32  3.20900e-32  3.16800e-32 ...

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         03/07/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_list = read_aeff_o_ii_ssb17_list(atom_rc_file)

    ii = np.where((abs(element_data_list.wavelength - wavelength) < 0.02)&(element_data_list.case1 == case1))[0]

    ii_length = len(ii)
    if len(ii) == 0:
        print('could not find the given wavelength')
        return 0
    if ii_length > 1:
        ii_min = min(abs(element_data_list[ii].wavelength - wavelength))
        ii = np.where((abs(element_data_list.wavelength - wavelength) == ii_min)&(element_data_list.case1 == case1))[0]
        ii_length = 1
    extension1 = np.asarray(element_data_list.extension)[ii]
    wavelength1 = np.asarray(element_data_list.wavelength)[ii]

    rc_element_template={'wavelength': float(0.0), 'aeff':[np.zeros((7,4))]}
    select_aeff_data = pd.DataFrame(data=rc_element_template, index=np.arange(ii_length))
    for i in range(0, ii_length):
        rc_aeff, header1 = getdata(atom_rc_file, extension1[i], header=True)
        #select_aeff_data.wavelength[i] = wavelength1[i]
        select_aeff_data.loc[i, 'wavelength']=wavelength1[i]
        # select_aeff_data.aeff[i] = rc_aeff
        select_aeff_data.loc[i, 'aeff'] = [rc_aeff]
    return select_aeff_data

def read_aeff_o_ii_ssb17_references(atom_rc_file):
    """
         This function returns the reference list of recombination coefficients (Aeff) from the 2nd binary table extension
         of the FITS data file ('rc_o_iii_SSB17.fits').

     :Private:

     :Returns:
        type=an array of data. This function returns the aeff_data_reference:
              { Reference:'',
                Citation:''}

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_o_iii_SSB17.fits')

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         03/07/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    reference_template = {'reference': '', 'citation': 0}
    data, hdr = getdata(atom_rc_file, 2, header=True)
    reference = data.ATOMICDATA
    citation = data.REFERENCE
    reference_length = len(reference)

    reference_data = pd.DataFrame(data=reference_template, index=np.arange(reference_length))
    # aij_data=aij_data.strip
    reference = np.char.strip(reference)
    citation = np.char.strip(citation)
    reference_data.reference = reference
    reference_data.citation = citation
    return reference_data

def list_aeff_o_ii_ssb17_references(atom_rc_file, atom, ion):
    """
         This function returns a list for all references of recombination coefficients (Aeff)
         for given element and ionic level from the FITS data file ('rc_o_iii_SSB17.fits').

     :Returns:
        type=an array of strings. This function returns the references.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_o_iii_SSB17.fits')
         atom          : in, required, type=string
                         atom name e.g. 'o'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = ['atomic-data-rc']
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_o_iii_SSB17_orl_case_b.fits')
         >>> atom='o'
         >>> ion='iii' # O II
         >>> list_oii_aeff_reference = atomneb.list_aeff_o_ii_ssb17_references(atom_rc_file, atom, ion)
         >>> print(list_oii_aeff_reference)


     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         03/07/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    rc_reference_data = read_aeff_n_ii_fsl13_references(atom_rc_file)

    atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff'
    rc_list = np.asarray(rc_reference_data.reference)
    #rc_res = [x for x in rc_list if re.search(atom_ion_name, x)]
    rc_res = [x for x in rc_list if atom_ion_name in x]
    ii_length = len(rc_res)
    if ii_length == 0:
        print('could not find the given element or ion')
        return 0
    nlines1 = ii_length  # temp[0]
    references = ['' for i in range(nlines1)]
    for i in range(0, nlines1):
        select_rc_data_str1 = rc_res[i].split('_')
        ref_num = len(select_rc_data_str1)
        if ref_num > 3:
            references[i] = select_rc_data_str1[3]
        else:
            references[i] = ''
    return references

def get_aeff_o_ii_ssb17_reference_citation(atom_rc_file, atom, ion, reference=None):
    """
         This function returns the reference citation for a recombination coefficient (Aeff)
         from the 2nd binary table extension of the FITS data file ('rc_o_iii_SSB17.fits').

     :Returns:
        type=string. This function returns the Citation.

     :Params:
         Atom_RC_file  : in, required, type=string
                         the FITS data file name ('rc_o_iii_SSB17.fits')
         atom          : in, required, type=string
                         atom name e.g. 'o'
         ion           : in, required, type=string
                         ionic level e.g 'iii'

     :Keywords:
         reference     : in, type=string
                         set for the reference e.g. 'SSB17',  may not necessary

     :Examples:
        For example::

         >>> import atomneb
         >>> import numpy as np
         >>> import os
         >>> base_dir = os.getcwd()
         >>> data_dir = ['atomic-data-rc']
         >>> atom_rc_file = os.path.join(base_dir,data_dir, 'rc_o_iii_SSB17_orl_case_b.fits')
         >>> atom='o'
         >>> ion='iii' # O II
         >>> reference='SSB17'
         >>> citation = atomneb.get_aeff_o_ii_ssb17_reference_citation(atom_rc_file, atom, ion)
         >>> print(citation)
            Storey, P.J., Sochi, T. and Bastin, R. 2017, MNRAS, 470, 379; VizieR On-line Data Catalog: VI/150

     :Categories:
       Recombination Lines

     :Dirs:
      ./
          Main routines

     :Author:
       Ashkbiz Danehkar

     :Copyright:
       This library is released under a GNU General Public License.

     :Version:
       0.2.0

     :History:
         03/07/2017, IDL code by A. Danehkar

         10/09/2020, A. Danehkar, Transferred from IDL to Python
    """
    element_data_reference = read_aeff_n_ii_fsl13_references(atom_rc_file)
    reference_list = np.asarray(element_data_reference.reference)
    atom_ion_name = atom.lower() + '_' + ion.lower() + '_aeff'
    if (reference is not None) == 1:
        atom_ion_name = atom_ion_name + '_' + reference.upper()

    ii = np.where(reference_list == atom_ion_name)[0]  #
    if len(ii) == 0:
        print('could not find the given reference')
        return 0
    citation = np.asarray(element_data_reference.citation[ii])

    return citation

