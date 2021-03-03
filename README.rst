======================
AtomNeb Python Package
======================

.. image:: https://img.shields.io/pypi/v/atomneb.svg?style=flat
    :target: https://pypi.python.org/pypi/atomneb/
    :alt: PyPI Version
    
.. image:: https://travis-ci.org/atomneb/AtomNeb-py.svg?branch=master
    :target: https://travis-ci.org/atomneb/AtomNeb-py
    :alt: Build Status
    
.. image:: https://ci.appveyor.com/api/projects/status/gi4ok3wy7jjn1ekb?svg=true
    :target: https://ci.appveyor.com/project/danehkar/atomneb-py
    :alt: Build Status
    
.. image:: https://coveralls.io/repos/github/atomneb/AtomNeb-py/badge.svg?branch=master
    :target: https://coveralls.io/github/atomneb/AtomNeb-py?branch=master
    :alt: Coverage Status
    
.. image:: https://img.shields.io/badge/license-GPL-blue.svg
    :target: https://github.com/atomneb/AtomNeb-py/blob/master/LICENSE
    :alt: GitHub license
    
.. image:: https://img.shields.io/conda/vn/conda-forge/atomneb.svg
    :target: https://anaconda.org/conda-forge/atomneb
    :alt: Anaconda Cloud
    
.. image:: https://readthedocs.org/projects/atomneb-py/badge/?version=latest
    :target: https://atomneb-py.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
    
.. image:: https://img.shields.io/badge/python-2.7%2C%203.8-blue.svg
    :alt: Support Python versions 2.7 and 3.8
    
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4287566.svg
    :target: https://doi.org/10.5281/zenodo.4287566
    :alt: Zenodo

.. image:: http://joss.theoj.org/papers/10.21105/joss.02797/status.svg
    :target: https://doi.org/10.21105/joss.02797
    :alt: JOSS

Description
===========

**AtomNeb-py** is a library written in Python for reading atomic data from *AtomNeb*, which is a database containing atomic data stored in the Flexible Image Transport System (FITS) file format for *collisionally excited lines* and *recombination lines* typically observed in spectra of ionized gaseous nebulae. The AtomNeb database were generated for use in `pyEQUIB <https://github.com/equib/pyEQUIB>`_, `proEQUIB <https://github.com/equib/proEQUIB>`_, and other nebular spectral analysis tools. 


Collisionally Excited Lines
---------------------------

*AtomNeb for collisionally excited lines*  contains sets of `atomic datasets <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data>`_, which include energy levels (Ej), collision strengths (Ωij), and transition probabilities (Aij) of the most ions commonly observed in ionized nebulae.

The atomic datasets for collisionally excited lines are as follows:

* `Collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/collection>`_ from the `National Institute of Standards and Technology (NIST) Atomic Spectra Database <https://www.nist.gov/pml/atomic-spectra-database>`_, the `CHIANTI atomic database <http://www.chiantidatabase.org/>`_, and some improved atomic data from `Cloudy v13.04 <https://www.nublado.org/>`_ and pyNeb v1.0. This collection was compiled according to the atomic data used in `pyNeb v1.0 <http://www.iac.es/proyecto/PyNeb/>`_.

* `Chianti52 <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti52>`_ from the `CHIANTI atomic database <http://www.chiantidatabase.org/>`_ version 5.2. This dataset was compiled according to the atomic data used in `MOCASSIN <https://github.com/mocassin/MOCASSIN-2.0>`_.

* `Chianti60 <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti60>`_ from the `CHIANTI atomic database <http://www.chiantidatabase.org/>`_ version 6.0. This dataset was compiled according to the atomic data used in `MOCASSIN <https://github.com/mocassin/MOCASSIN-2.0>`_.

* `Chianti70 <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti70>`_ from the `CHIANTI atomic database <http://www.chiantidatabase.org/>`_ version 7.0. This dataset was compiled according to the atomic data used in `MOCASSIN <https://github.com/mocassin/MOCASSIN-2.0>`_.

* `Chianti90 <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti90>`_ from the `CHIANTI atomic database <http://www.chiantidatabase.org/>`_ version 9.0. This dataset was compiled according to the atomic data used in `NEAT <https://github.com/rwesson/NEAT>`_.

Each dataset contains the following `atomic data FITS files <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti70>`_: ``AtomElj.fits`` for *Energy Levels* (Ej), ``AtomOmij.fits`` for *Collision Strengths* (Ωij), and ``AtomAij.fits`` for *Transition Probabilities* (Aij).


Recombination Lines
-------------------

*AtomNeb for recombination lines* contains sets of `effective recombination coefficients <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_ (αeff) of recombination lines of H I, He I, He II, C I, C II, C III, C VI, N II, N III, N IV, N V, N VI, N VII, O II, O III, O IV, O V, O VI, O VIII, and Ne II ions typically observed in ionized nebulae, as well as Branching ratios (Br) of O II and N II lines.

The atomic datasets for recombination lines are as follows:

* `RC Collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, effective recombination coefficients for C II (`Davey et al. 2000 <http://adsabs.harvard.edu/abs/2000A%26AS..142...85D>`_), N II (`Escalante and Victor 1990 <http://adsabs.harvard.edu/abs/1990ApJS...73..513E>`_), O II (`Storey 1994 <http://adsabs.harvard.edu/abs/1994A%26A...282..999S>`_; `Liu et al. 1995 <http://adsabs.harvard.edu/abs/1995MNRAS.272..369L>`_), and Ne II ions (`Kisielius et al. 1998 <http://adsabs.harvard.edu/abs/1998A%26AS..133..257K>`_), including Branching ratios (Br) for O II and N II ions. This collection was compiled according to the atomic data used in `MOCASSIN <https://github.com/mocassin/MOCASSIN-2.0>`_.

* `SH95 Collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, hydrogenic ions for Z=1 to 8, namely H I, He II, Li III, Be IV, B V, C VI, N VII, and O VIII ions from `Storey and Hummer (1995) <http://adsabs.harvard.edu/abs/1995MNRAS.272...41S>`_.

* `PPB91 Collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, effective recombination coefficients for H, He, C, N, O, Ne ions from `Pequignot, Petitjean and Boisson (1991) <http://adsabs.harvard.edu/abs/1991A%26A...251..680P>`_.

* `PFSD12 He I data <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, effective He I recombination coefficients from `Porter et al (2012) <http://adsabs.harvard.edu/abs/2012MNRAS.425L..28P>`_ and `(2013a) <http://adsabs.harvard.edu/abs/2013MNRAS.433L..89P>`_.

* `FSL13 N II data <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, effective N II recombination coefficients (corrigendum) from `Fang, Storey and Liu (2011) <http://adsabs.harvard.edu/abs/2011A%26A...530A..18F>`_ and `(2013b) <http://adsabs.harvard.edu/abs/2013A%26A...550C...2F>`_.

* `SSB17 O II data <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, effective O II recombination coefficients of 8889 recombination lines for Cases A, B, and C, and 2433 optical (3500-9000Å) recombination lines for Case B from `Storey, Sochi and Bastin (2017) <http://adsabs.harvard.edu/abs/2017MNRAS.470..379S>`_.



Installation
============

Dependent Python Packages
-------------------------

 This package requires the following packages:

    - `NumPy <https://numpy.org/>`_
    - `Astropy <https://www.astropy.org/>`_

The previous version relied on `pandas <https://pandas.pydata.org/>`_, but all the data structures were changed from pandas.DataFrame to those defined by `NumPy <https://numpy.org/>`_ that speed up the computations and reduce the memory usage.
    
* To get this package with all the FITS file, you can simply use ``git`` command as follows::

        git clone https://github.com/atomneb/AtomNeb-py

* If you plan to use the recent O II recombination coefficients (`Storey, Sochi and Bastin 2017 <http://adsabs.harvard.edu/abs/2017MNRAS.470..379S>`_), you need to unpack them::

        cd AtomNeb-py/atomic-data-rc/
        tar -xvf *.fits.tar.gz


To install the last version, all you should need to do is

.. code-block::

    $ python setup.py install

To install the stable version, you can use the preferred installer program (pip):

.. code-block::

    $ pip install atomneb

or you can install it from the cross-platform package manager *conda*:

.. code-block::

    $ conda install -c conda-forge atomneb

How to Use
==========

The Documentation of the functions provides in detail in the *API Documentation* (`atomneb.github.io/AtomNeb-py/doc <https://atomneb.github.io/AtomNeb-py/doc>`_). There are two main categories: *collisionally excited lines (CEL)* and *recombination lines (RC)*.

* The atomic data for **collisionally excited lines (CEL)** contain Energy Levels (Ej), Collision Strengths (Ωij), and Transition Probabilities (Aij). We have four atomic datasets for them: `collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/collection>`_, `chianti52 <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti52>`_, `chianti60 <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti60>`_, and `chianti70 <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti70>`_. 
    
    You need to load the **atomneb** library as follows::
    
        import atomneb
     
    Also::

        import numpy as np
        import os
        
        atom_elj_file = os.path.join(base_dir,data_dir, 'AtomElj.fits')
        atom_omij_file = os.path.join(base_dir,data_dir, 'AtomOmij.fits')
        atom_aij_file = os.path.join(base_dir,data_dir, 'AtomAij.fits')
        elj_data_list = atomneb.read_elj_list(atom_elj_file)
        omij_data_list = atomneb.read_omij_list(atom_omij_file)
        aij_data_list = atomneb.read_aij_list(atom_aij_file)
   
    Now you have access to:
     
    - *Energy Levels* (Ej)::
    
        atom='o'
        ion='iii'
        oiii_elj_data = atomneb.read_elj(atom_elj_file, atom, ion, level_num=6)
        print(oiii_elj_data['j_v'])
        print(oiii_elj_data['ej'])
    
      which gives::
    
        0.00000      1.00000      2.00000      2.00000      0.00000      2.00000
        0.00000      113.200      306.200      20273.30     43185.69     60324.80
    
    - *Collision Strengths* (Ωij)::
    
        atom='o'
        ion='iii'
        oiii_omij_data = atomneb.read_omij(atom_omij_file, atom, ion)
        print(oiii_omij_data['level1'])
        print(oiii_omij_data['level2'])
        print(oiii_omij_data['strength'][0])
    
      which gives::
        
        0       1       1       1       1       ...
        0       2       3       4       5       ...
        100.0      158.50       251.20       398.10       631.0       ...
    
    - *Transition Probabilities* (Aij)::
    
        atom='o'
        ion='iii'
        oiii_aij_data = atomneb.read_aij(atom_aij_file, atom, ion)
        print(oiii_aij_data['aij'][0])
    
      which gives::
        
         0.0000   2.5969e-05       0.0000   2.3220e-06      ...
    
* The atomic data for **recombination lines (RC)** contain effective recombination coefficients (αeff) of emission lines from different collections: `RC Collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, `SH95 Collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, `PPB91 Collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, `PFSD12 He I data <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, `FSL13 N II data <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, and `SSB17 O II data <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_.
    
    You need to load the **atomneb** libary::
    
    
        import atomneb
     
    Also:

        import numpy as np
        import os
    
    Now you have access to effective recombination coefficients (αeff) of the following collections:
     
    - *RC Collection*::

        atom_rc_file = os.path.join(base_dir,data_dir, 'rc_collection.fits')
        atom='c'
        ion='iii'
        cii_rc_data = atomneb.read_aeff_collection(atom_rc_file, atom, ion)
        n_line = len(cii_rc_data['wavelength'])
        for i in range(0, n_line):
             print(cii_rc_data['wavelength'][i], cii_rc_data['a'][i], 
                   cii_rc_data['b'][i], cii_rc_data['c'][i], 
                   cii_rc_data['d'][i], cii_rc_data['f'][i])
        
      which gives::
    
        914.00000      0.69280000     0.021400000    -0.016300000     -0.24310000     -0.88000000
        962.00000       1.0998000   -0.0042000000    -0.027900000     -0.22940000     -0.96560000
        997.00000      0.78210000     -0.36840000   0.00030000000     -0.12170000     -0.78740000
        ...
        
    - *SH95 Collection*::
    
        atom_rc_file = os.path.join(base_dir,data_dir, 'rc_SH95.fits')
        atom='h'
        ion='ii'
        hi_rc_data = atomneb.read_aeff_sh95(atom_rc_file, atom, ion)
        print(hi_rc_data['aeff'][0])
        
      which gives::
    
        100.00000       500.00000       0.0000000   4.2140000e-27   1.7560000e-27   1.0350000e-27
        ...
        
    - *PPB91 Collection*::
    
        atom_rc_file = os.path.join(base_dir,data_dir, 'rc_PPB91.fits')
        atom='c'
        ion='iii'
        cii_rc_data = atomneb.read_aeff_ppb91(atom_rc_file, atom, ion)
        n_line = len(cii_rc_data['wavelength'])
        for i in range(0, n_line):
           print(cii_rc_data['ion'][i], cii_rc_data['case1'][i], cii_rc_data['wavelength'][i],
                 cii_rc_data['a'][i], cii_rc_data['b'][i], cii_rc_data['c'][i],
                 cii_rc_data['d'][i], cii_rc_data['br'][i], cii_rc_data['q'][i], cii_rc_data['y'][i])
           
      which gives::
    
        C2+A       9903.4600      0.69700000     -0.78400000       4.2050000      0.72000000       1.0000000       1.6210000
        C2+A       4267.1500       1.0110000     -0.75400000       2.5870000      0.71900000      0.95000000       2.7950000
        ...
          
    - *PFSD12 He I data*::

        atom_rc_file = os.path.join(base_dir,data_dir, 'rc_he_ii_PFSD12.fits')
        atom='he'
        ion='ii'
        hei_rc_data = atomneb.read_aeff_he_i_pfsd12(atom_rc_file, atom, ion)
        hei_rc_data_wave = atomneb.read_aeff_he_i_pfsd12(atom_rc_file, atom, ion, wavelength=True)
        print(hei_rc_data['aeff'][0])
           
      which gives::
    
        5000.0000       10.000000      -25.379540      -25.058970      -25.948440      -24.651820      -25.637660     
        ...
        
    - *FSL13 N II data*::
    
        atom_rc_file = os.path.join(base_dir,data_dir, 'rc_n_iii_FSL13.fits')
        atom='n'
        ion='iii'
        wavelength_range=[4400.0, 7100.0] 
        nii_rc_data = atomneb.read_aeff_n_ii_fsl13(atom_rc_file, atom, ion, wavelength_range)
        nii_rc_data_wave = atomneb.read_aeff_n_ii_fsl13(atom_rc_file, atom, ion, wavelength_range, wavelength=True)
        print(nii_rc_data.aeff[0])
        n_line = len(hei_rc_data_wave['wavelength'])
        for i in range(0, n_line):
           print(nii_rc_data_wave['wavelength'][i], nii_rc_data_wave['tr'][i], nii_rc_data_wave['trans'][i])
        
      which gives::
    
        255.000      79.5000      47.3000      12.5000      6.20000      4.00000      2.86000
        258.000      54.4000      29.7000      7.92000      4.11000      2.72000      2.00000
        310.000      48.1000      23.7000      5.19000      2.55000      1.65000      1.21000
        434.000      50.3000      23.2000      4.71000      2.26000      1.45000      1.05000
          
        6413.23 6g - 4f2p6g G[9/2]o4 - 2p4f F[7/2]e3
        6556.32 6g - 4f2p6g G[9/2]o5 - 2p4f G[7/2]e4
        6456.97 6g - 4f2p6g G[9/2]o5 - 2p4f F[7/2]e4
        6446.53 6g - 4f2p6g F[7/2]o3 - 2p4f D[5/2]e2
        6445.34 6g - 4f2p6g F[7/2]o4 - 2p4f D[5/2]e3
        ...
        
    - *SSB17 O II data*: You first need to unpack rc_o_iii_SSB17_orl_case_b.fits.tar.gz, e.g.:: 

        tar -xvf rc_o_iii_SSB17_orl_case_b.fits.tar.gz

      If you need to have access to the full dataset (entire wavelengths, case A and B)::

        tar -xvf rc_o_iii_SSB17.fits.tar.gz


      Please note that using the entire atomic data will make your program very slow and you may need to have a higher memory on your system. Without the above comment, as default, the cose uses rc_o_iii_SSB17_orl_case_b.fits::

        aatom_rc_file = os.path.join(base_dir,data_dir, 'rc_o_iii_SSB17_orl_case_b.fits')
        atom='o'
        ion='iii'
        case1='B'
        wavelength_range=[5320.0, 5330.0] 
        oii_rc_data = atomneb.read_aeff_o_ii_ssb17(atom_rc_file, atom, ion, case1, wavelength_range)
        oii_rc_data_wave = atomneb.read_aeff_o_ii_ssb17(atom_rc_file, atom, ion, case1, wavelength_range, wavelength=True)
        print(oii_rc_data['aeff'][0])
        n_line = len(oii_rc_data_wave.wavelength)
        for i in range(0, n_line):
           print(oii_rc_data_wave.wavelength[i], oii_rc_data_wave.lower_term[i], oii_rc_data_wave.upper_term[i])
        
      which gives::
    
        1.64100e-30  1.60000e-30  1.56400e-30  1.54100e-30  1.52100e-30  1.50900e-30
        ...
          
        5327.17 2s22p2(1S)3p 2Po
        5325.42 2s22p2(1S)3p 2Po
        5327.18 2s22p2(1D)3d 2Ge
        5326.84 2s22p2(1D)3d 2Ge
        ...


Documentation
=============

For more information on how to use the API functions from the AtomNeb Python package, please read the `API Documentation  <https://atomneb.github.io/AtomNeb-py/doc>`_ published on `atomneb.github.io/AtomNeb-py <https://atomneb.github.io/AtomNeb-py>`_.


References
==========

* Danehkar, A. (2020). AtomNeb Python Package, an addendum to AtomNeb: IDL Library for Atomic Data of Ionized Nebulae. *J. Open Source Softw.*, **5**, 2797. doi:`10.21105/joss.02797 <https://doi.org/10.21105/joss.02797>`_ ads:`2020JOSS....5.2797D <https://ui.adsabs.harvard.edu/abs/2020JOSS....5.2797D>`_.

* Danehkar, A. (2019). AtomNeb: IDL Library for Atomic Data of Ionized Nebulae. *J. Open Source Softw.*, **4**, 898. doi:`10.21105/joss.00898 <https://doi.org/10.21105/joss.00898>`_ ads:`2019JOSS....4..898D <https://ui.adsabs.harvard.edu/abs/2019JOSS....4..898D>`_.


Citation
========

Using **AtomNeb** in a scholarly publication? Please cite these papers:

.. code-block:: bibtex

   @article{Danehkar2020,
     author = {{Danehkar}, Ashkbiz},
     title = {AtomNeb Python Package, an addendum to AtomNeb: IDL Library for Atomic Data of Ionized Nebulae},
     journal = {Journal of Open Source Software},
     volume = {5},
     number = {55},
     pages = {2797},
     year = {2020},
     doi = {10.21105/joss.02797}
   }

   @article{Danehkar2019,
     author = {{Danehkar}, Ashkbiz},
     title = {AtomNeb: IDL Library for Atomic Data of Ionized Nebulae},
     journal = {Journal of Open Source Software},
     volume = {4},
     number = {35},
     pages = {898},
     year = {2019},
     doi = {10.21105/joss.00898}
   }

Learn More
==========

==================  =============================================
**Documentation**   https://atomneb-py.readthedocs.io/
**Repository**      https://github.com/atomneb/AtomNeb-py
**Issues & Ideas**  https://github.com/atomneb/AtomNeb-py/issues
**Conda-Forge**     https://anaconda.org/conda-forge/atomneb
**PyPI**            https://pypi.org/project/atomneb/
**DOI**             `10.21105/joss.02797 <https://doi.org/10.21105/joss.02797>`_
**Archive**         `10.5281/zenodo.4287566 <https://doi.org/10.5281/zenodo.4287566>`_
==================  =============================================
