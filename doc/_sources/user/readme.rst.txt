Introduction
============

**atomneb** is a library written in Python for reading atomic data from a database containing atomic data stored in the Flexible Image Transport System (FITS) file format for *collisionally excited lines* and *recombination lines* typically observed in spectra of ionized gaseous nebulae. The AtomNeb database were generated for use in `pyEQUIB <https://github.com/equib/pyEQUIB>`_, `proEQUIB <https://github.com/equib/proEQUIB>`_, and other nebular spectral analysis tools. 


Collisional Excitation Unit
---------------------------

*AtomNeb for collisionally excited lines*  contains sets of `atomic datasets <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data>`_, which include energy levels (Ej), collision strengths (Ωij), and transition probabilities (Aij) of the most ions commonly observed in ionized nebulae.

The atomic datasets for collisionally excited lines are as follows:

* `Collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/collection>`_ from the `National Institute of Standards and Technology (NIST) Atomic Spectra Database <https://www.nist.gov/pml/atomic-spectra-database>`_, the `CHIANTI atomic database <http://www.chiantidatabase.org/>`_, and some improved atomic data from `Cloudy v13.04 <https://www.nublado.org/>`_ and pyNeb v1.0. This collection was compiled according to the atomic data used in `pyNeb v1.0 <http://www.iac.es/proyecto/PyNeb/>`_.

* `Chianti52 <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti52>`_ from the `CHIANTI atomic database <http://www.chiantidatabase.org/>`_ version 5.2. This dataset was compiled according to the atomic data used in `MOCASSIN <https://github.com/mocassin/MOCASSIN-2.0>`_.

* `Chianti60 <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti60>`_ from the `CHIANTI atomic database <http://www.chiantidatabase.org/>`_ version 6.0. This dataset was compiled according to the atomic data used in `MOCASSIN <https://github.com/mocassin/MOCASSIN-2.0>`_.

* `Chianti70 <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti70>`_ from the `CHIANTI atomic database <http://www.chiantidatabase.org/>`_ version 7.0. This dataset was compiled according to the atomic data used in `MOCASSIN <https://github.com/mocassin/MOCASSIN-2.0>`_.

* `Chianti90 <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti90>`_ from the `CHIANTI atomic database <http://www.chiantidatabase.org/>`_ version 9.0. This dataset was compiled according to the atomic data used in `NEAT <https://github.com/rwesson/NEAT>`_.

Each dataset contains the following `atomic data FITS files <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti70>`_: ``AtomElj.fits`` for *Energy Levels* (Ej), ``AtomOmij.fits`` for *Collision Strengths* (Ωij), and ``AtomAij.fits`` for *Transition Probabilities* (Aij).


Recombination Unit
------------------

*AtomNeb for recombination lines* contains sets of `effective recombination coefficients <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_ (αeff) of recombination lines of H I, He I, He II, C I, C II, C III, C VI, N II, N III, N IV, N V, N VI, N VII, O II, O III, O IV, O V, O VI, O VIII, and Ne II ions typically observed in ionized nebulae, as well as Branching ratios (Br) of O II and N II lines.

The atomic datasets for recombination lines are as follows:

* `RC Collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, effective recombination coefficients for C II (`Davey et al. 2000 <http://adsabs.harvard.edu/abs/2000A%26AS..142...85D>`_), N II (`Escalante and Victor 1990 <http://adsabs.harvard.edu/abs/1990ApJS...73..513E>`_), O II (`Storey 1994 <http://adsabs.harvard.edu/abs/1994A%26A...282..999S>`_; `Liu et al. 1995 <http://adsabs.harvard.edu/abs/1995MNRAS.272..369L>`_), and Ne II ions (`Kisielius et al. 1998 <http://adsabs.harvard.edu/abs/1998A%26AS..133..257K>`_), including Branching ratios (Br) for O II and N II ions. This collection was compiled according to the atomic data used in `MOCASSIN <https://github.com/mocassin/MOCASSIN-2.0>`_.

* `SH95 Collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, hydrogenic ions for Z=1 to 8, namely H I, He II, Li III, Be IV, B V, C VI, N VII, and O VIII ions from `Storey and Hummer (1995) <http://adsabs.harvard.edu/abs/1995MNRAS.272...41S>`_.

* `PPB91 Collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, effective recombination coefficients for H, He, C, N, O, Ne ions from `Pequignot, Petitjean and Boisson (1991) <http://adsabs.harvard.edu/abs/1991A%26A...251..680P>`_.

* `PFSD12 He I data <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, effective He I recombination coefficients from `Porter et al (2012) <http://adsabs.harvard.edu/abs/2012MNRAS.425L..28P>`_ and `(2013a) <http://adsabs.harvard.edu/abs/2013MNRAS.433L..89P>`_.

* `FSL13 N II data <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, effective N II recombination coefficients (corrigendum) from `Fang, Storey and Liu (2011) <http://adsabs.harvard.edu/abs/2011A%26A...530A..18F>`_ and `(2013b) <http://adsabs.harvard.edu/abs/2013A%26A...550C...2F>`_.

* `SSB17 O II data <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, effective O II recombination coefficients of 8889 recombination lines for Cases A, B, and C, and 2433 optical (3500-9000Å) recombination lines for Case B from `Storey, Sochi and Bastin (2017) <http://adsabs.harvard.edu/abs/2017MNRAS.470..379S>`_.

