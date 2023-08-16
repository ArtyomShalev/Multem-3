import calculating_lib as cl
import numpy as np


figures_to_calculate = ['fig3']
# ----- calculating data for fig. 2 a) -------------------
if 'fig2a' in figures_to_calculate:
    input_params = {
        'mts': '0', #multipole_type_selected
        'mos': '0', #multipole_order_selected
        'mps': '0', #m_projection_selected
        'type': '1 1',
        'order': '1 1',
        'm': '1 -1',
        'rmax': 7,
        'lmax': 4,
        'lattice_constant': 300,
        'fab': 60,
        'polar': 'S',
        'epssph_re': 15,
        'epssph_im': 0,
        'epsmed_re': 1,
        'epsmed_im': 0,
        'ktype': 2,
        'kscan': 1,
        'zinf': 3.73,
        'zsup': 3.74,
        'npts': 1000,
        'r_ratio': 0.4705,
        'mode': '1', # 2D array of spheres
        'multem_version': '3'
    }
    
    
    omega = np.linspace(input_params['zinf'], input_params['zsup'], input_params['npts'])
    ak1 = np.array((1e-2, 0, 5e-2, 1e-1))/2/np.pi
    ak2 = np.array((0, 5e-2, 0, 0))/2/np.pi
    kpts = len(ak1)
    lmax = 4
    print('calculating data for fig. 2 a)...')
    for i in range(kpts):
        F, T, R, A = cl.calc_spectrum_omega(omega, ak1[i], ak2[i], input_params)        
        cl.save_1D_data(omega, T, dir='data/fig2a', filename=f'akxy={round(ak1[i]*2*np.pi, 4)}_{round(ak2[i]*2*np.pi, 4)}.txt', format='%19.16e')
    print('done')
    

    # ----- calculating data for fig. 2 b) -------------------
if 'fig2b' in figures_to_calculate:
    print('calculating data for fig. 2 b)...')
    input_params = {
        'mts': '0', #multipole_type_selected
        'mos': '0', #multipole_order_selected
        'mps': '0', #m_projection_selected
        'type': '1 1',
        'order': '1 1',
        'm': '1 -1',
        'rmax': 7,
        'lattice_constant': 300,
        'fab': 60,
        'polar': 'S',
        'epssph_re': 15,
        'epssph_im': 0,
        'epsmed_re': 1,
        'epsmed_im': 0,
        'ktype': 2,
        'kscan': 1,
        'zinf': 3.7315,
        'zsup': 3.7356,
        'npts': 2000,
        'r_ratio': 0.4705,
        'mode': '1', # 2D array of spheres
        'multem_version': '3'
    }
    
    omega = np.linspace(input_params['zinf'], input_params['zsup'], input_params['npts'])
    ak1 = np.array((0.01))/2/np.pi
    ak2 = np.array((0.0))/2/np.pi
    LMAX = [4, 7, 10, 13]
    for lmax in LMAX:
        input_params['lmax'] = lmax
        F, T, R, A = cl.calc_spectrum_omega(omega, ak1, ak2, input_params)        
        cl.save_1D_data(omega, T, dir='data/fig2b', filename=f'lmax={lmax}.txt', format='%19.16e')       
    print('done')

  
    # ----- calculating data for fig. 2 c) -------------------
if 'fig2c' in figures_to_calculate:
    print('calculating data for fig. 2 c)...')
    input_params = {
        'mts': '0', #multipole_type_selected
        'mos': '0', #multipole_order_selected
        'mps': '0', #m_projection_selected
        'type': '1 1',
        'order': '1 1',
        'm': '1 -1',
        'rmax': 7,
        'lmax': 10,
        'lattice_constant': 300,
        'fab': 60,
        'polar': 'S',
        'epssph_re': 15,
        'epssph_im': 0,
        'epsmed_re': 1,
        'epsmed_im': 0,
        'ktype': 2,
        'kscan': 1,
        'zinf': 3.731,
        'zsup': 3.735,
        'npts': 2000,
        'r_ratio': 0.4705,
        'mode': '1' # 2D array of spheres
    }
    
    omega = np.linspace(input_params['zinf'], input_params['zsup'], input_params['npts'])
    ak1 = np.array((0.01))/2/np.pi
    ak2 = np.array((0.0))/2/np.pi
    input_params['multem_version'] = '3'
    # T_lmax=10
    F, T, R, A = cl.calc_spectrum_omega(omega, ak1, ak2, input_params)        
    cl.save_1D_data(omega, T, dir='data/fig2c', filename='T.txt', format='%19.16e') 
    # with lapack and faddeeva
    F, T, R, A = cl.calc_spectrum_omega(omega, ak1, ak2, input_params)        
    cl.save_1D_data(omega, (T+R-1), dir='data/fig2c', filename='error_with_lapack_and_faddeeva.txt', format='%19.16e') 
    # wo lapack
    input_params['multem_version'] = 'wo_lapack'
    F, T, R, A = cl.calc_spectrum_omega(omega, ak1, ak2, input_params)        
    cl.save_1D_data(omega, (T+R-1), dir='data/fig2c', filename='error_wo_lapack.txt', format='%19.16e')       
    # with lapack
    input_params['multem_version'] = 'with_lapack'    
    F, T, R, A = cl.calc_spectrum_omega(omega, ak1, ak2, input_params)        
    cl.save_1D_data(omega, (T+R-1), dir='data/fig2c', filename='error_with_lapack.txt', format='%19.16e') 
    print('done')


if 'fig6' in figures_to_calculate:
    import time

    input_params = {
    'mts': '0', #multipole_type_selected
    'mos': '0', #multipole_order_selected
    'mps': '0', #m_projection_selected
    'type': '1',
    'order': '1',
    'm': '1',
    'rmax': 2,
    'lmax': 4,
    'lattice_constant': 300,
    'fab': 60,
    'polar': 'S',
    'epssph_re': 15,
    'epssph_im': 0,
    'epsmed_re': 1,
    'epsmed_im': 0,
    'ktype': 2,
    'kscan': 1,
    'zinf': 0.5,
    'zsup': 1.0,
    'npts': 1000,
    'r_ratio': 0.4705,
    'mode': '1', #multi-layered
    'dist_btw_spheres_and_interface': 0.6
}

    #fig6 a-d
    omega = np.linspace(3.734, 3.737, input_params['npts'])
    input_params['zinf'] = omega[0]
    input_params['zsup'] = omega[-1]
    ak1 = 0.01/2/np.pi
    ak2 = 0
    for input_params['multem_version'] in ['3_cerf', '3']:
        for input_params['lmax'] in [4]:
            time0 = time.time()
            for input_params['rmax'] in [7]:
                F, T, R, A = cl.calc_spectrum_omega(omega, ak1, ak2, input_params)
                cl.save_1D_data(omega, T, dir=f'data/fig4/{input_params["multem_version"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}_ak1={round(ak1, 8)}.txt', format='%19.16e') 
            print(f'{input_params["multem_version"]} : {time.time()-time0}s')

    #fig6 b-e
    omega = np.linspace(3.73555, 3.73562, input_params['npts'])
    input_params['zinf'] = omega[0]
    input_params['zsup'] = omega[-1]
    ak1 = 8e-5/2/np.pi
    ak2 = 0
    for input_params['multem_version'] in ['3_cerf', '3']:
        for input_params['lmax'] in [4]:
            time0 = time.time()
            for input_params['rmax'] in [7]:
                F, T, R, A = cl.calc_spectrum_omega(omega, ak1, ak2, input_params)
                cl.save_1D_data(omega, T, dir=f'data/fig4/{input_params["multem_version"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}_ak1={round(ak1, 8)}.txt', format='%19.16e') 
            print(f'{input_params["multem_version"]} : {time.time()-time0}s')
    
    #fig6 c-f
    omega = np.linspace(3.73558, 3.735585, input_params['npts'])
    input_params['zinf'] = omega[0]
    input_params['zsup'] = omega[-1]
    ak1 = 4e-5/2/np.pi
    ak2 = 0
    for input_params['multem_version'] in ['3_cerf', '3']:
        for input_params['lmax'] in [4]:
            time0 = time.time()
            for input_params['rmax'] in [7]:
                F, T, R, A = cl.calc_spectrum_omega(omega, ak1, ak2, input_params)
                cl.save_1D_data(omega, T, dir=f'data/fig4/{input_params["multem_version"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}_ak1={round(ak1, 8)}.txt', format='%19.16e') 
            print(f'{input_params["multem_version"]} : {time.time()-time0}s')


if 'fig3' in figures_to_calculate:
    import time

    input_params = {
    'mts': '0', #multipole_type_selected
    'mos': '0', #multipole_order_selected
    'mps': '0', #m_projection_selected
    'type': '1',
    'order': '1',
    'm': '1',
    'rmax': 2,
    'lmax': 4,
    'lattice_constant': 300,
    'fab': 60,
    'polar': 'S',
    'epssph_re': 15,
    'epssph_im': 0,
    'epsmed_re': 1,
    'epsmed_im': 0,
    'ktype': 2,
    'kscan': 1,
    'zinf': 0.5,
    'zsup': 1.0,
    'npts': 100,
    'r_ratio': 0.4705,
    'mode': '3', #multi-layered
    'dist_btw_spheres_and_interface': 1.0
}
    zoomed = True
    omega = np.linspace(3.73, 3.74, input_params['npts'])
    input_params['zinf'] = omega[0]
    input_params['zsup'] = omega[-1]
    ak1 = 0.01/2/np.pi
    ak2 = 0
    for input_params['multem_version'] in ['3']:
        for input_params['lmax'] in [13]:
            time0 = time.time()
            for input_params['rmax'] in [8, 22, 26, 27]:
                F, T, R, A = cl.calc_spectrum_omega(omega, ak1, ak2, input_params)
                if zoomed:
                    cl.save_1D_data(omega, T, dir=f'data/fig5/mode={input_params["mode"]}/version={input_params["multem_version"]}/zoomed//d={input_params["dist_btw_spheres_and_interface"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}.txt', format='%19.16e') 
                else:
                    cl.save_1D_data(omega, T, dir=f'data/fig5/mode={input_params["mode"]}/d={input_params["dist_btw_spheres_and_interface"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}.txt', format='%19.16e') 

            print(f'{input_params["multem_version"]} : {time.time()-time0}s')


if 'fig4' in figures_to_calculate:
    import time
    input_parameters = {
        'ktype': 1,
        'lmax': 3,
        'rmax': 20,
        'lattice_constant': 475,
        'radius': 100,
        'epssph_re': -20.1480000,
        'epssph_im': 1.24700000,
        'epsmed_re': 2.1025,
        'epsmed_im': 0,
        'polar': 'S',
        'fab': 60,
        'mts': '0',
        'mos': '0',
        'mps': '0',
        'type': '1',
        'order': '1',
        'm': '1',
        'multem_version': '3'
    }

    input_parameters['r_ratio'] = input_parameters['radius']/input_parameters['lattice_constant']
    fi = 0
    from_sin_theta = 0.0
    to_sin_theta = 0.999
    n_theta = 1000
    sin_theta = np.linspace(from_sin_theta, to_sin_theta, n_theta)
    theta = np.arcsin(sin_theta)*180/np.pi

    for input_parameters['rmax'] in [7, 8, 9]:
        dir = f'data/fig4/'
        factor = 0
        for wl in [650, 750, 900]:
            if wl == 650:
                input_parameters['epssph_re'] = -12.953
                input_parameters['epssph_im'] = 1.1209
            elif wl == 750:
                input_parameters['epssph_re'] = -20.148
                input_parameters['epssph_im'] = 1.2470
            elif wl == 900:
                input_parameters['epssph_re'] = -32.719
                input_parameters['epssph_im'] = 1.9955
            else:
                print('epspsh set not correct!')
            R = cl.calc_spectrum_ak1(ktype=1, ak1=theta, wl=wl, ip=input_parameters)
            cl.save_result(dir, f'{wl}_{input_parameters["rmax"]}', R)
            factor += 1
        cl.save_result(dir, 'sintheta', sin_theta)


if 'diff_orders_test' in figures_to_calculate:
    input_parameters = {
        'ktype': 1,
        'lmax': 3,
        'rmax': 5,
        'lattice_constant': 475,
        'radius': 100,
        'epssph_re': -20.1480000,
        'epssph_im': 1.24700000,
        'epsmed_re': 2.1025,
        'epsmed_im': 0,
        'polar': 'S',
        'fab': 60,
        'mts': '0',
        'mos': '0',
        'mps': '0',
        'type': '1',
        'order': '1',
        'm': '1',
        'multem_version': '3'
    }
    input_parameters['r_ratio'] = input_parameters['radius']/input_parameters['lattice_constant']
    from_sin_theta = 0.0
    to_sin_theta = 0.999
    n_theta = 100
    sin_theta = np.linspace(from_sin_theta, to_sin_theta, n_theta)
    theta = np.arcsin(sin_theta)*180/np.pi
    R = cl.calc_spectrum_ak1(ktype=1, ak1=theta, wl=600, ip=input_parameters)
