---
title: "AtomNeb Python Package, an addendum to AtomNeb: IDL Library for Atomic Data of Ionized Nebulae"
tags:
  - python
  - astrophysics
  - atomic data
  - gaseous nebulae
  - spectral analysis
authors:
  - name: Ashkbiz Danehkar
    orcid: 0000-0003-4552-5997
    affiliation: "1, 2, 3"
affiliations:
 - name: Department of Physics and Astronomy, Macquarie University, Sydney, NSW 2109, Australia
   index: 1
 - name: Harvard-Smithsonian Center for Astrophysics, 60 Garden Street, Cambridge, MA 02138, USA 
   index: 2
 - name: Department of Astronomy, University of Michigan, 1085 S. University Avenue, Ann Arbor, MI 48109, USA 
   index: 3
date: 20 October 2020
bibliography: paper.bib
---

# Addendum

`AtomNeb` is a Python open-source package containing atomic data for gaseous nebulae stored in the Flexible Image Transport System (FITS) file format [@Wells:1981; @Hanisch:2001; @Pence:2010]. These FITS files offer easy access to the atomic data required for emissivity calculations in the collisional excitation and recombination processes usually occurred in ionized gases of planetary nebulae and H II regions. This package has several application programming interface (API) functions developed in Python for retrieving the energy levels, collision strengths, transition probabilities, and recombination coefficients from its FITS files. The previous library `AtomNeb` [@Danehkar:2019] coupled to the library `proEQUIB` [@Danehkar:2018b] needs the Interactive Data Language (IDL) compiler, so this package offers an identical package for the high-level programming language Python that can be used by those astrophysicists, who intend to analyze nebular emission lines by developing codes in Python. The `AtomNeb` Python functions can be used by the Python package `pyEQUIB` [@Danehkar:2020] to analyze emission-line spectra. 

`AtomNeb` uses the FITS handling routines of the Python package `Astropy` [@Astropy:2013; @Astropy:2018] to retrieve the atomic data from its FITS files. It also requires the Python packages `NumPy` [@Walt:2011; @Harris:2020] and `pandas` [@McKinney:2010; @McKinney:2011; @McKinney:2017]. This package is released under the GNU General Public License, and its source code is publicly available on its GitHub repository. Its latest version can be installed directly from its repository on the GitHub, and the stable version from the Python Package Index (PyPi) via ``pip install atomneb`` or alternatively from the Conda Python package manager via ``conda install -c conda-forge atomneb``. The online documentation, tutorials and examples are provided on the GitHub platform (https://github.com/atomneb/AtomNeb-py) and the Read the Docs documentation host (https://atomneb-py.readthedocs.io/).


# Acknowledgements

AD acknowledges the support of Research Excellence Scholarship from Macquarie University.

# References
