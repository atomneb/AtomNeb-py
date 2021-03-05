Usage
=====

.. role:: python(code)
   :language: python

The Documentation of the functions provides in detail in the *API Documentation* (`atomneb.github.io/AtomNeb-py/doc <https://atomneb.github.io/AtomNeb-py/doc>`_). There are two main categories: *collisionally excited lines (CEL)* and *recombination lines (RC)*.

Collisional Excitation Unit
---------------------------

The atomic data for **collisional excitation unit (CEL)** contain Energy Levels (Ej), Collision Strengths (Ωij), and Transition Probabilities (Aij). We have four atomic datasets for them: `collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/collection>`_, `chianti52 <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti52>`_, `chianti60 <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti60>`_, and `chianti70 <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data/chianti70>`_. 
    
    You need to load the **atomneb** library as follows::
    
        import atomneb
     
    Also::

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
    
    * *Collision Strengths* (Ωij)::
    
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
    
    * *Transition Probabilities* (Aij)::
    
        atom='o'
        ion='iii'
        oiii_aij_data = atomneb.read_aij(atom_aij_file, atom, ion)
        print(oiii_aij_data['aij'][0])
    
      which gives::
        
         0.0000   2.5969e-05       0.0000   2.3220e-06      ...

Recombination Unit
------------------   

The atomic data for **recombination unit (RC)** contain effective recombination coefficients (αeff) of emission lines from different collections: `RC Collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, `SH95 Collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, `PPB91 Collection <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, `PFSD12 He I data <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, `FSL13 N II data <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_, and `SSB17 O II data <https://github.com/atomneb/AtomNeb-py/tree/master/atomic-data-rc>`_.
    
    You need to load the **atomneb** libary::
    
    
        import atomneb
     
    Also::

        import numpy as np
        import os
    
    Now you have access to effective recombination coefficients (αeff) of the following collections:
     
    * *RC Collection*::

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
        
    * *SH95 Collection*::
    
        atom_rc_file = os.path.join(base_dir,data_dir, 'rc_SH95.fits')
        atom='h'
        ion='ii'
        hi_rc_data = atomneb.read_aeff_sh95(atom_rc_file, atom, ion)
        print(hi_rc_data['aeff'][0])
        
      which gives::
    
        100.00000       500.00000       0.0000000   4.2140000e-27   1.7560000e-27   1.0350000e-27
        ...
        
    * *PPB91 Collection*::
    
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
          
    * *PFSD12 He I data*::

        atom_rc_file = os.path.join(base_dir,data_dir, 'rc_he_ii_PFSD12.fits')
        atom='he'
        ion='ii'
        hei_rc_data = atomneb.read_aeff_he_i_pfsd12(atom_rc_file, atom, ion)
        hei_rc_data_wave = atomneb.read_aeff_he_i_pfsd12(atom_rc_file, atom, ion, wavelength=True)
        print(hei_rc_data['aeff'][0])
           
      which gives::
    
        5000.0000       10.000000      -25.379540      -25.058970      -25.948440      -24.651820      -25.637660     
        ...
        
    * *FSL13 N II data*::
    
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
        
    * *SSB17 O II data*: You first need to unpack rc_o_iii_SSB17_orl_case_b.fits.tar.gz, e.g.:: 

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
        n_line = len(oii_rc_data_wave['wavelength'])
        for i in range(0, n_line):
           print(oii_rc_data_wave['wavelength'][i], oii_rc_data_wave['lower_term'][i], oii_rc_data_wave['upper_term'][i])
        
      which gives::
    
        1.64100e-30  1.60000e-30  1.56400e-30  1.54100e-30  1.52100e-30  1.50900e-30
        ...
          
        5327.17 2s22p2(1S)3p 2Po
        5325.42 2s22p2(1S)3p 2Po
        5327.18 2s22p2(1D)3d 2Ge
        5326.84 2s22p2(1D)3d 2Ge
        ...

