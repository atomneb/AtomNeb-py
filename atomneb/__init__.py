"""
atomneb - Python Package for Atomic Data of Ionized Nebulae
"""

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

from .atomneb import read_aij,read_aij_list,search_aij,read_aij_references,list_aij_references, \
           get_aij_reference_citation, \
           read_elj,read_elj_list,read_elj_references,get_elj_reference_citation, \
           read_omij,read_omij_list,search_omij,read_omij_references, \
           list_omij_references,get_omij_reference_citation, \
           read_aeff_collection,read_aeff_collection_list, \
           search_aeff_collection,read_aeff_collection_references, \
           list_aeff_collection_references,get_aeff_collection_reference_citation, \
           read_aeff_sh95,read_aeff_sh95_list,search_aeff_sh95,read_aeff_sh95_references, \
           list_aeff_sh95_references,get_aeff_sh95_reference_citation, \
           read_aeff_ppb91,read_aeff_ppb91_list,search_aeff_ppb91,read_aeff_ppb91_references, \
           list_aeff_ppb91_references,get_aeff_ppb91_reference_citation, \
           read_aeff_he_i_pfsd12,read_aeff_he_i_pfsd12_list,search_aeff_he_i_pfsd12, \
           read_aeff_he_i_pfsd12_references,list_aeff_he_i_pfsd12_references, \
           get_aeff_he_i_pfsd12_reference_citation, \
           read_aeff_n_ii_fsl13,read_aeff_n_ii_fsl13_list,search_aeff_n_ii_fsl13, \
           read_aeff_n_ii_fsl13_references,list_aeff_n_ii_fsl13_references, \
           get_aeff_n_ii_fsl13_reference_citation, \
           read_aeff_o_ii_ssb17,read_aeff_o_ii_ssb17_list,search_aeff_o_ii_ssb17, \
           read_aeff_o_ii_ssb17_references,list_aeff_o_ii_ssb17_references, \
           get_aeff_o_ii_ssb17_reference_citation

from .version import __version__

import sys
from numpy.version import version as numpy_version

if sys.version_info[0:2] < (2, 6):
    log_.warn('atomneb requires Python version >= 2.6, but it is version {0}'.format(sys.version_info), calling='atomneb')
try:
    if [int(n) for n in (numpy_version.split('.')[:3])] < [1, 5, 1] :
        log_.warn('atomneb Numpy version >= 1.5.1, but it is version {0}'.format(numpy_version), calling='atomneb')
except:
    log_.warn('Cannot find Numpy version {0}, report the bug'.format(numpy_version), calling='pyemcee')
    


