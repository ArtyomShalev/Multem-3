'''
this script is used to run and test Multem-3 simulation 
there are different options how to run Multem-3:
! important: multem3 binary file should be in the directory with fort.10 and mutipole_regime_prameters.ini
--------------------------------------------------------
1.  manually build this project: run cmake and make
    copy fort.10 to directory Multem-3/bin
    run multem3 bin file (Multem-3/bin)

2.  manually build this project: run cmake and make
    run this script 

--------------------------------------------------------
'''
import time
import numpy as np
import sys
sys.path.insert(0, '..')
import src.python_simulations.multem3_py_calculating as calc
import shutil


#triangular lattice of dieletric spheres 
input_params = {
        'mts': '0', #multipole_type_selected
        'mos': '0', #multipole_order_selected
        'mps': '0', #m_projection_selected
        'type': '1',
        'order': '1',
        'm': '1',
        'rmax': 8,
        'lmax': 4,
        'lattice_constant': 300,
        'fab': 60,
        'polar': 'S',
        'epssph_re': 15,
        'epssph_im': 0,
        'epsmed_re': 1,
        'epsmed_im': 0,
        'ktype': 1,
        'kscan': 1,
        'npts': 300,
        'r_ratio': 0.4705,
        'mode': '1',
        'multem_version': '3',
        'nlayer': '1'
    }

omega = np.linspace(3.7275, 3.74, 300)/2/np.pi
kx = np.linspace(0, 0.5, 100)
time0 = time.time()
F, T, R, A = calc.calc_spectrum_omega(omega, 2, 0.01, input_params)
print(f'transmittance\n{T}')
print(f'calculation time: {time.time()-time0}s')

