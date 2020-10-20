---
title: "AtomNeb: Python Package for Atomic Data of Ionized Nebulae"
tags:
  - python
  - astrophysics
  - atomic data
  - gaseous nebulae
  - spectral analysis
authors:
  - name: Ashkbiz Danehkar
    orcid: 0000-0003-4552-5997
    affiliation: 1, 2, 3
affiliations:
 - name: Department of Physics and Astronomy, Macquarie University, Sydney, NSW 2109, Australia
   index: 1
 - name: Harvard-Smithsonian Center for Astrophysics, 60 Garden Street, Cambridge, MA 02138, USA 
   index: 2
 - name: Department of Astronomy, University of Michigan, 1085 S. University Avenue, Ann Arbor, MI 48109, USA 
   index: 3
date: 20 October 2018
bibliography: paper.bib
---

# Introduction

Emission-line spectra emitted from gaseous nebulae can be used to probe the physical and chemical proprieties of the interstellar medium in our Galaxy and other galaxies. The collisional excitation and recombination processes in ionized gases produces collisionally excited lines (CEL) and recombination lines (RL), respectively. The spectral analysis of CELs can be performed by solving statistical equilibrium equations, which requires the energy levels ($E_{j}$), collision strengths ($\Omega_{ij}$), and transition probabilities ($A_{ij}$) of ions, and yields atomic level populations and emissivities. The results of these calculations are employed to determine physical conditions and chemical abundances from CELs [see e.g. @Danehkar:2013; @Danehkar:2014; @Danehkar:2014b]. Moreover, the emissivities of RLs can be determined from effective recombination coefficients ($\alpha_{\rm eff}$) of ions, which are utilized to estimate chemical abundances from RLs.

# Statement of Need

The atomic data of the collisional excitation and recombination processes are essential for spectral analysis tools and photoionization programs. The collisional excitation calculations were implemented in some spectral analysis packages [e.g. @Howarth:1981; @Shaw:1994; @Shaw:1998; @Wesson:2012; @Luridiana:2015; @Howarth:2016; @Danehkar:2018b; @Danehkar:2020] and photoionization programs [e.g. @Ferland:1998; @Kallman:2001; @Ercolano:2003; @Ercolano:2005; @Ercolano:2008], which require to have the collision strengths ($\Omega_{ij}$) and transition probabilities ($A_{ij}$). Moreover, the recombination calculations in spectral analysis and photoionization codes [e.g. @Ferland:1998; @Kallman:2001; @Ercolano:2003; @Ercolano:2005; @Ercolano:2008; @Wesson:2012; @Luridiana:2015; @Danehkar:2018b; @Danehkar:2020] are performed using effective recombination coefficients ($\alpha_{\rm eff}$). The previous library `AtomNeb` [@Danehkar:2019] needs the Interactive Data Language (IDL) compiler, so it was necessary to develop an identical package for the high-level programming language Python that can be used by those astrophysicists, who intend to analyze nebular emission lines by developing codes in Python.

# Description

`AtomNeb` is a Python open-source package containing atomic data for gaseous nebulae stored in the Flexible Image Transport System (FITS) file format [@Wells:1981; @Hanisch:2001; @Pence:2010]. 
These FITS files offer easy access to the atomic data required for emissivity calculations in the collisional excitation and recombination processes usually occurred in ionized gases of planetary nebulae and H II regions. This package has several application programming interface (API) functions developed in Python for retrieving the energy levels, collision strengths, transition probabilities, and recombination coefficients from its FITS files. This package has two units, namely:

- the _collisional excitation unit_ that contains API functions for providing the energy levels ($E_{j}$), collision strengths  ($\Omega_{ij}$), and transition probabilities ($A_{ij}$) of ions that are used to calculate atomic level populations and emissivities in the collisional excitation process. These API functions can be employed to obtain the electron temperature, electron concentration, and chemical abundances from observed fluxes of CELs. The collisional excitation unit includes atomic data from the CHIANTI database version 5.2 [@Landi:2006], version 6.0 [@Dere:2009], version 7.0 [@Landi:2012], and version 9.0 [@Dere:2019], which were compiled based on the atomic data used in the photoionization code `MOCASSIN` [@Ercolano:2003; @Ercolano:2005; @Ercolano:2008] and the nebular empirical analysis tool `NEAT` [@Wesson:2012]. It also includes the atomic data used in the Python emission-line analysis package `pyNeb` [@Luridiana:2015].

- the _recombination unit_ that contains API functions for providing _effective recombination coefficients_ ($\alpha_{\rm eff}$) and _branching ratios_ ($Br$) of ions that are used to compute emissivities in the recombination process. The recombination unit includes effective recombination coefficients for C II [@Davey:2000], N II [@Escalante:1990], O II [@Storey:1994; @Liu:1995], and Ne II [@Kisielius:1998], which were used in the photoionization code `MOCASSIN`. It also contains recombination coefficients of hydrogenic ions for Z=1 to 8 [@Storey:1995], H, He, C, N, O, and Ne ions [@Pequignot:1991], He I [@Porter:2012; @Porter:2013], N II [@Fang:2011; @Fang:2013], and O II [@Storey:2017].

`AtomNeb` uses the FITS handling routines of the Python package `Astropy` [@Astropy:2013; @Astropy:2018] to retrieve the atomic data from its FITS files. It also requires the Python packages `NumPy` [@Walt:2011; @Harris:2020] and `pandas` [@McKinney:2010; @McKinney:2011; @McKinney:2017]. API functions of the IDL library `AtomNeb` [@Danehkar:2019] are used by the IDL library `proEQUIB` [@Danehkar:2018b] to conduct plasma diagnostics and abundance analysis of emission-line spectra from ionized nebulae [see e.g. @Danehkar:2016; @Danehkar:2018a]. Similarly, the `AtomNeb` Python functions can be used by the Python package `pyEQUIB` [@Danehkar:2020] to analyze emission-line spectra. This package is released under the GNU General Public License, and its source code is publicly available on its GitHub repository. Its latest version can be installed directly from its repository on the GitHub, and the stable version from the Python Package Index (PyPi) via ``pip install atomneb`` or alternatively from the Conda Python package manager via ``conda install -c conda-forge atomneb``. The online documentation, tutorials and examples are provided on the GitHub platform (https://github.com/atomneb/AtomNeb-py) and the Read the Docs documentation host (https://atomneb-py.readthedocs.io/).


# Acknowledgements

AD acknowledges the support of Research Excellence Scholarship from Macquarie University.

# References
