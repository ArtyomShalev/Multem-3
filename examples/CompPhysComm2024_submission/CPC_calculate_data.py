import numpy as np
import sys
sys.path.insert(0, '../..')
import src.python_simulations.multem3_py_calculating as calc
import time 
import matplotlib.pyplot as plt


def get_half_width_half_maxima_and_x0(x, y):
        #TODO corner cases
        idx = [np.argmax(y), np.argmin(y)]
        ymax = y[idx[0]]
        ymin = y[idx[1]]
        amp = ymax - ymin
        idx_hwhm = np.abs(y - 0.5*amp).argmin()
        x0 = x[idx[0]]
        hwhm = np.abs(x0 - x[idx_hwhm])

        return hwhm, x0
        
# fig1 - flowchart [draw.io]
# fig2 - typical system design [POV-ray]
# figures_to_calculate = ['fig4', 'fig5', 'fig8', 'fig9', 'fig10', 'fig6']
figures_to_calculate = ['fig8', 'fig9']
 
    # ----- calculating data for fig. 3 -------------------
if 'fig3' in figures_to_calculate:
    print('fig3 data calculation...')
    start = time.time()
    input_params = {
        'mts': '0', #multipole_type_selected
        'mos': '0', #multipole_order_selected
        'mps': '0', #m_projection_selected
        'type': '1 1',
        'order': '1 1',
        'm': '1 -1',
        'rmax': 16,
        'lmax': 7,
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
        'mode': '1', # 2D array of spheres
        'nlayer': '1'
    }
    
    omega = np.linspace(input_params['zinf'], input_params['zsup'], input_params['npts'])
    ak1 = np.array((0.01))/2/np.pi
    ak2 = np.array((0.0))/2/np.pi
    input_params['multem_version'] = '3'
    F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)        
    calc.save_1D_data(omega, T, dir='data/fig3', filename='T.txt', format='%19.16e') 
    # wo lapack
    input_params['multem_version'] = '2'
    F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)        
    calc.save_1D_data(omega, 1-(T+R), dir='data/fig3', filename='error_wo_lapack.txt', format='%19.16e')       
    # with lapack
    input_params['multem_version'] = 'with_lapack'    
    F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)        
    calc.save_1D_data(omega, 1-(T+R), dir='data/fig3', filename='error_with_lapack.txt', format='%19.16e') 
    print('done')
    print(f'CPU time:{time.time()-start}')


if 'fig4' in figures_to_calculate:
    print('fig4 data calculation...')
    start = time.time()
    input_params = {
        'mts': '0', #multipole_type_selected
        'mos': '0', #multipole_order_selected
        'mps': '0', #m_projection_selected
        'type': '1 1',
        'order': '1 1',
        'm': '1 -1',
        'rmax': 16,
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
        'multem_version': '3',
        'nlayer': '1',
        'lmax': 4
    }
    dir_to_save = 'data/fig4'
    omega = np.linspace(input_params['zinf'], input_params['zsup'], input_params['npts'])
    ak1 = np.array((1e-2, 0, 5e-2, 1e-1))/2/np.pi
    ak2 = np.array((0, 5e-2, 0, 0))/2/np.pi
    kpts = len(ak1)
    lmax = 4
    for i in range(kpts):
        F, T, R, A = calc.calc_spectrum_omega(omega, ak1[i], ak2[i], input_params)        
        calc.save_1D_data(omega, T, dir=dir_to_save, filename=f'akxy={round(ak1[i]*2*np.pi, 4)}_{round(ak2[i]*2*np.pi, 4)}.txt', format='%19.16e')
    #fig 4 b
    input_params['zinf'] = 3.7315
    input_params['zsup'] = 3.7356
    omega = np.linspace(input_params['zinf'], input_params['zsup'], input_params['npts'])
    ak1 = np.array((0.01))/2/np.pi
    ak2 = np.array((0.0))/2/np.pi
    LMAX = [4, 7, 10, 13]
    for lmax in LMAX:
        input_params['lmax'] = lmax
        F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)        
        calc.save_1D_data(omega, T, dir=dir_to_save, filename=f'lmax={lmax}.txt', format='%19.16e')       
    print(f'CPU time:{time.time()-start}')


if 'fig5' in figures_to_calculate:
    print('fig5 data calculation...')
    start = time.time()
    input_parameters = {
        'ktype': 1,
        'lmax': 3,
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

    for input_parameters['rmax'] in [7, 12, 14]:
        dir = f'data/fig5/'
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
                print('epspsh is not set correctly!')
            R = calc.calc_spectrum_ak1(ktype=1, ak1=theta, wl=wl, ip=input_parameters)
            calc.save_result(dir, f'{wl}_{input_parameters["rmax"]}', R)
            factor += 1
        calc.save_result(dir, 'sintheta', sin_theta)
        print(f"CPU time {time.time()-start}")


if 'fig6' in figures_to_calculate:
    print('fig6 data calculation...')
    start = time.time()
    input_params = {
    'mts': '0', #multipole_type_selected
    'mos': '0', #multipole_order_selected
    'mps': '0', #m_projection_selectedfig4
    'type': '1',
    'order': '1',
    'm': '1',
    'rmax': 2,
    'lmax': 4,
    'lattice_constant': 500,
    'fab': 60,
    'polar': 'S',
    # 'epssph_re': 2.2231,
    # 'epssph_re': 6,
    'epssph_re': 25,
    'epssph_im': 0.00,
    'epsmed_re': 0.00,
    'epsmed_im': 0.00,
    'ktype': 2,
    'kscan': 1,
    'zinf': 0.5,
    'zsup': 1.0,
    'npts': 300,
    'r_ratio': 0.48,
    # 'r_ratio': 0.5,
    # 'r_ratio': 0.48683,
    'mode': '3_ref', #multi-layered
    'dist_btw_spheres_and_interface': 0.41
}
    # omega = np.linspace(6.2, 6.292, input_params['npts'])
    # omega = np.linspace(6.1, 6.2, input_params['npts'])
    # omega = np.linspace(1, 2.5, input_params['npts'])
    # omega = np.linspace(2.5, 3.5, input_params['npts'])

    # omega = np.linspace(0.8, 1.2, input_params['npts'])
    omega = np.linspace(0.2, 1.8, input_params['npts'])
    # omega = np.linspace(0.01, 10, input_params['npts'])
    

    #just to check
    # dl = [0.22, 0.0, input_params['dist_btw_spheres_and_interface']]
    # dr = [0.22, 0.0, input_params['dist_btw_spheres_and_interface']]

    # print(f'vector length:',examples/CompPhysComm2024_submission/multem_different_versions/multem3_cerf/bin/multem3_cerfnp.sum(np.power(dl + dr, 2)))

    # dl = [0.33, 0.33, input_params['dist_btw_spheres_and_interface']]
    # dr = [0.33, 0.33, input_params['dist_btw_spheres_and_interface']]

    # print(f'vector length:', np.sum(p.power(dl + dr, 2)))

    # print(f'minimal vector length:', input_params['r_ratio'])
    NLAYERS = ['3', '4', '5', '6', '7', '8', '9']
    # NLAYERS = ['7', '8', '9']
    # NLAYERS = ['3', '4', '6']

    for nlayer in NLAYERS:
        nunit = f'_1unit_{nlayer}'
        input_params['nlayer'] = nlayer
        input_params['zinf'] = omega[0]
        input_params['zsup'] = omega[-1]
        ak1 = 0
        ak2 = 0
        for input_params['multem_version'] in ['3']:
            for input_params['lmax'] in [7]:
                time0 = time.time()
                # for input_params['rmax'] in [32]:
                for input_params['rmax'] in [20, 34]:
                    F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)
                    calc.save_1D_data(omega, R, dir=f'data/fig6/{nunit}/mode={input_params["mode"]}/version={input_params["multem_version"]}/d={input_params["dist_btw_spheres_and_interface"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}.txt', format='%19.16e') 
                print(f'{input_params["multem_version"]} : {time.time()-time0}s')
        #10738.055797576904s
        for input_params['multem_version'] in ['2']:
            for input_params['lmax'] in [7]:
                time0 = time.time()
                # for input_params['rmax'] in [35, 50]:
                for input_params['rmax'] in [20, 34]:
                    F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)
                    calc.save_1D_data(omega, R, dir=f'data/fig6/{nunit}/mode={input_params["mode"]}/version={input_params["multem_version"]}/d={input_params["dist_btw_spheres_and_interface"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}.txt', format='%19.16e') 
                print(f'{input_params["multem_version"]} : {time.time()-time0}s')
    # for input_params['multem_version'] in ['2', '3']:
    #     for input_params['lmax'] in [7]:
    #         time0 = time.time()
    #         for input_params['rmax'] in [20, 30]:
    #             F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)
    #             calc.save_1D_data(omega, T, dir=f'data/fig6/mode={input_params["mode"]}/version={input_params["multem_version"]}/d={input_params["dist_btw_spheres_and_interface"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}.txt', format='%19.16e') 
    #         print(f'{input_params["multem_version"]} : {time.time()-time0}s')
    print(f'CPU time:{time.time()-start}')
   

#========================================================================
#fig 7
#calculation is implemented in dir Zaghoul_MIT_faddeeva_comparison
# fig 7 contains data calculated in fig3a, so one needs to calculate data for fig 3a first
# then run bash script cerf_faddeeva_comparison/run_CompPhysComm_data_prep.sh
#========================================================================


if 'fig8' in figures_to_calculate:
    print('fig8 data calculation...')
    input_params = {
    'mts': '0', #multipole_type_selected
    'mos': '0', #multipole_order_selected
    'mps': '0', #m_projection_selected
    'type': '1',
    'order': '1',
    'm': '1',
    'rmax': 16,
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
    'r_ratio': 0.470512,
    'mode': '1', #multi-layered
    'dist_btw_spheres_and_interface': 0.6,
    'nlayer': '1'
}
    start = time.time()
    dir_to_save = 'data/fig8/'
    omega = np.linspace(3.7354, 3.7357, input_params['npts'])
    input_params['zinf'] = omega[0]
    input_params['zsup'] = omega[-1]
    ak1 = 1e-3/2/np.pi
    ak2 = 0
    for input_params['multem_version'] in ['3', '3_cerf']:
        time0 = time.time()
        F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)
        calc.save_1D_data(omega, R, dir=f'{dir_to_save}{input_params["multem_version"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}_ak1={round(ak1, 8)}.txt', format='%19.16e') 
        print(f'{input_params["multem_version"]} : {time.time()-time0}s')

    omega = np.linspace(3.7354, 3.7357, input_params['npts'])
    input_params['zinf'] = omega[0]
    input_params['zsup'] = omega[-1]
    ak1 = 8e-5/2/np.pi
    ak2 = 0
    for input_params['multem_version'] in ['3', '3_cerf']:
        time0 = time.time()
        F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)
        calc.save_1D_data(omega, R, dir=f'{dir_to_save}{input_params["multem_version"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}_ak1={round(ak1, 8)}.txt', format='%19.16e') 
        print(f'{input_params["multem_version"]} : {time.time()-time0}s')
    
    omega = np.linspace(3.73553338, 3.735533395, input_params['npts']) #lmax4
    input_params['zinf'] = omega[0]
    input_params['zsup'] = omega[-1]
    ak1 = 8e-6/2/np.pi
    ak2 = 0
    for input_params['multem_version'] in ['3_cerf', '3']:
        time0 = time.time()
        F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)
        calc.save_1D_data(omega, R, dir=f'{dir_to_save}{input_params["multem_version"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}_ak1={round(ak1, 8)}.txt', format='%19.16e') 
        print(f'{input_params["multem_version"]} : {time.time()-time0}s')
        # ----- fano resonance fitting -----------
        # hwhm, x0 = get_half_width_half_maxima_and_x0(omega, R)
        # print(f'fwhm={hwhm*2}, x0={x0}')
        # plt.plot(omega, R)
        # plt.scatter(x0, 1)
        # plt.show()
    print(f'CPU time: {time.time()-start}')


if 'fig9' in figures_to_calculate:
    print('fig9 data calculation...')
    input_params = {
    'mts': '0', #multipole_type_selected
    'mos': '0', #multipole_order_selected
    'mps': '0', #m_projection_selected
    'type': '1',
    'order': '1',
    'm': '1',
    'rmax': 16,
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
    'zinf': 0.5,
    'zsup': 1.0,
    'npts': 1000,
    # 'r_ratio': 0.470512,
    'r_ratio': 0.470665,
    'mode': '1', #multi-layered
    'dist_btw_spheres_and_interface': 0.6,
    'nlayer': '1'
}
    start = time.time()
    omega = np.linspace(3.7314, 3.732, input_params['npts'])
    input_params['zinf'] = omega[0]
    input_params['zsup'] = omega[-1]
    ak1 = 1e-3/2/np.pi
    ak2 = 0
    for input_params['multem_version'] in ['3_cerf', '3']:
        time0 = time.time()
        F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)
        calc.save_1D_data(omega, R, dir=f'data/fig9/{input_params["multem_version"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}_ak1={round(ak1, 8)}.txt', format='%19.16e') 
        print(f'{input_params["multem_version"]} : {time.time()-time0}s')

    omega = np.linspace(3.7316, 3.7317, input_params['npts'])
    input_params['zinf'] = omega[0]
    input_params['zsup'] = omega[-1]
    ak1 = 8e-5/2/np.pi
    ak2 = 0
    for input_params['multem_version'] in ['3_cerf', '3']:
        time0 = time.time()
        F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)
        calc.save_1D_data(omega, R, dir=f'data/fig9/{input_params["multem_version"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}_ak1={round(ak1, 8)}.txt', format='%19.16e') 
        print(f'{input_params["multem_version"]} : {time.time()-time0}s')

    omega = np.linspace(3.73163055, 3.731630575, input_params['npts'])
    input_params['zinf'] = omega[0]
    input_params['zsup'] = omega[-1]
    ak1 = 8e-6/2/np.pi
    ak2 = 0
    for input_params['multem_version'] in ['3_cerf', '3']:
        time0 = time.time()
        F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)
        calc.save_1D_data(omega, R, dir=f'data/fig9/{input_params["multem_version"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}_ak1={round(ak1, 8)}.txt', format='%19.16e') 
        print(f'{input_params["multem_version"]} : {time.time()-time0}s')
        #--------- BIC spectral position check -------------
        # hwhm, x0 = get_half_width_half_maxima_and_x0(omega, R)
        # print(f'fwhm={hwhm*2}, x0={x0}')
        # plt.plot(omega, R)
        # plt.scatter(x0, 1)
        # plt.show()
    print(f'CPU time: {time.time()-start}')


if 'fig10' in figures_to_calculate:
    print('fig10 data calculation...')
    start = time.time()
    input_parameters = {
        'ktype': 1,
        'lmax': 2,
        'rmax': 10,
        'radius': 100,
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
    dir = f'data/fig10/'
    fi = 0
    from_sin_theta = 0.08
    to_sin_theta = 0.12
    n_theta = 300
    sin_theta = np.linspace(from_sin_theta, to_sin_theta, n_theta)
    theta = np.arcsin(sin_theta)*180/np.pi
    input_parameters['lattice_constant'] = 475.0
    input_parameters['r_ratio'] = input_parameters['radius']/input_parameters['lattice_constant']
    wl = 650
    input_parameters['zinf'] = wl/input_parameters['lattice_constant']
    input_parameters['zsup'] = (wl+0.01)/input_parameters['lattice_constant']
    input_parameters['epssph_re'] = -12.953
    input_parameters['epssph_im'] = 1.1209
    # total
    input_parameters['mts'] = '0'
    input_parameters['mos'] = '0'
    input_parameters['mps'] = '0'
    R = calc.calc_spectrum_ak1(ktype=1, ak1=theta, wl=wl, ip=input_parameters)
    calc.save_1D_data(sin_theta, R, dir, f'{wl}_m_total.txt')
    # q_{sz} multipole
    input_parameters['mts'] = '1'
    input_parameters['mos'] = '1'
    input_parameters['mps'] = '1'
    input_parameters['type'] = '0 0'
    input_parameters['order'] = '2 2'
    input_parameters['m'] = '-1 1'
    R = calc.calc_spectrum_ak1(ktype=1, ak1=theta, wl=wl, ip=input_parameters)  
    calc.save_1D_data(sin_theta, R, dir, f'{wl}_m_qsz.txt')
    # m{k} + q{sz}
    input_parameters['type'] =  '0 0 1 1'
    input_parameters['order'] = '2 2 1 1'
    input_parameters['m'] =     '-1 1 1 -1'
    R = calc.calc_spectrum_ak1(ktype=1, ak1=theta, wl=wl, ip=input_parameters)
    calc.save_1D_data(sin_theta, R, dir, f'{wl}_m_qsz_mk.txt')
    wl = 750
    from_sin_theta = 0.307
    to_sin_theta = 0.313
    sin_theta = np.linspace(from_sin_theta, to_sin_theta, n_theta)
    theta = np.arcsin(sin_theta)*180/np.pi
    input_parameters['zinf'] = wl/input_parameters['lattice_constant']
    input_parameters['zsup'] = (wl+0.01)/input_parameters['lattice_constant']
    input_parameters['epssph_re'] = -20.148
    input_parameters['epssph_im'] = 1.2470
    # total
    input_parameters['mts'] = '0'
    input_parameters['mos'] = '0'
    input_parameters['mps'] = '0'
    R = calc.calc_spectrum_ak1(ktype=1, ak1=theta, wl=wl, ip=input_parameters)
    calc.save_1D_data(sin_theta, R, dir, f'{wl}_m_total.txt')
    # q_{sz} multipole
    input_parameters['mts'] = '1'
    input_parameters['mos'] = '1'
    input_parameters['mps'] = '1'
    input_parameters['type'] = '0 0'
    input_parameters['order'] = '2 2'
    input_parameters['m'] = '-1 1'
    R = calc.calc_spectrum_ak1(ktype=1, ak1=theta, wl=wl, ip=input_parameters)  
    calc.save_1D_data(sin_theta, R, dir, f'{wl}_m_qsz.txt')
    # m{k} + q{sz}
    input_parameters['type'] =  '0 0 1 1'
    input_parameters['order'] = '2 2 1 1'
    input_parameters['m'] =     '-1 1 1 -1'
    R = calc.calc_spectrum_ak1(ktype=1, ak1=theta, wl=wl, ip=input_parameters)
    calc.save_1D_data(sin_theta, R, dir, f'{wl}_m_qsz_mk.txt')

    from_sin_theta = 0.0
    to_sin_theta = 0.4
    sin_theta = np.linspace(from_sin_theta, to_sin_theta, n_theta)
    theta = np.arcsin(sin_theta)*180/np.pi
    # total
    input_parameters['mts'] = '0'
    input_parameters['mos'] = '0'
    input_parameters['mps'] = '0'
    R = calc.calc_spectrum_ak1(ktype=1, ak1=theta, wl=wl, ip=input_parameters)
    calc.save_1D_data(sin_theta, R, dir, f'{wl}_p_total.txt')
    # p_{s} multipole
    input_parameters['mts'] = '1'
    input_parameters['mos'] = '1'
    input_parameters['mps'] = '1'
    input_parameters['type'] = '0 0'
    input_parameters['order'] = '1 1'
    input_parameters['m'] = '-1 1'
    R = calc.calc_spectrum_ak1(ktype=1, ak1=theta, wl=wl, ip=input_parameters)  
    calc.save_1D_data(sin_theta, R, dir, f'{wl}_p_ps.txt')
    # p{s} + m{z} + q{ks}
    input_parameters['type'] =  '0 0 1 0 0'
    input_parameters['order'] = '1 1 1 2 2'
    input_parameters['m'] =     '-1 1 0 -2 2'
    R = calc.calc_spectrum_ak1(ktype=1, ak1=theta, wl=wl, ip=input_parameters)
    calc.save_1D_data(sin_theta, R, dir, f'{wl}_p_qks_mz_ps.txt')

    wl = 900
    from_sin_theta = 0.4
    to_sin_theta = 0.8
    sin_theta = np.linspace(from_sin_theta, to_sin_theta, n_theta)
    theta = np.arcsin(sin_theta)*180/np.pi
    input_parameters['zinf'] = wl/input_parameters['lattice_constant']
    input_parameters['zsup'] = (wl+0.01)/input_parameters['lattice_constant']
    input_parameters['epssph_re'] = -32.719
    input_parameters['epssph_im'] = 1.9955
    # total
    input_parameters['mts'] = '0'
    input_parameters['mos'] = '0'
    input_parameters['mps'] = '0'
    R = calc.calc_spectrum_ak1(ktype=1, ak1=theta, wl=wl, ip=input_parameters)
    calc.save_1D_data(sin_theta, R, dir, f'{wl}_p_total.txt')
    # p_{s} multipole
    input_parameters['mts'] = '1'
    input_parameters['mos'] = '1'
    input_parameters['mps'] = '1'
    input_parameters['type'] = '0 0'
    input_parameters['order'] = '1 1'
    input_parameters['m'] = '-1 1'
    R = calc.calc_spectrum_ak1(ktype=1, ak1=theta, wl=wl, ip=input_parameters)  
    calc.save_1D_data(sin_theta, R, dir, f'{wl}_p_ps.txt')
    # p{s} + m{z} + q{ks}
    input_parameters['type'] =  '0 0 1 0 0'
    input_parameters['order'] = '1 1 1 2 2'
    input_parameters['m'] =     '-1 1 0 -2 2'
    R = calc.calc_spectrum_ak1(ktype=1, ak1=theta, wl=wl, ip=input_parameters)
    calc.save_1D_data(sin_theta, R, dir, f'{wl}_p_qks_mz_ps.txt')
    print(f'CPU time: {time.time()-start}')

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
    R = calc.calc_spectrum_ak1(ktype=1, ak1=theta, wl=600, ip=input_parameters)



