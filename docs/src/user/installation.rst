Installation
============

To install the last version, all you should need to do is::

        python setup.py install

To install the stable version, you can use the preferred installer program (pip)::

        pip install atomneb

or you can install it from the cross-platform package manager *conda*::

        conda install -c conda-forge atomneb

To get this package with all the FITS file, you can simply use ``git`` command as follows::

        git clone https://github.com/atomneb/AtomNeb-py

If you plan to use the recent O II recombination coefficients (`Storey, Sochi and Bastin 2017 <http://adsabs.harvard.edu/abs/2017MNRAS.470..379S>`_), you need to unpack them::

        cd AtomNeb-py/atomic-data-rc/
        tar -xvf *.fits.tar.gz
        
This package requires the following packages:

    - `NumPy <https://numpy.org/>`_
    - `pandas <https://pandas.pydata.org/>`_
    - `Astropy <https://www.astropy.org/>`_
