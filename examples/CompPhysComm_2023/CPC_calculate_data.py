import numpy as np
import sys
sys.path.insert(0, '../..')
import src.python_simulations.multem3_py_calculating as calc
import time 

# 'fig3a', 'fig3b', 'fig3c', 'fig5'
figures_to_calculate = ['fig8']

# ----- calculating data for fig. 3 a) -------------------
if 'fig3' in figures_to_calculate:
    start = time.time()
    input_params = {
        'mts': '0', #multipole_type_selected
        'mos': '0', #multipole_order_selected
        'mps': '0', #m_projection_selected
        'type': '1 1',
        'order': '1 1',
        'm': '1 -1',
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
    print('calculating data for fig. 3 a)...')
    for i in range(kpts):
        F, T, R, A = calc.calc_spectrum_omega(omega, ak1[i], ak2[i], input_params)        
        calc.save_1D_data(omega, T, dir='data/fig3', filename=f'akxy={round(ak1[i]*2*np.pi, 4)}_{round(ak2[i]*2*np.pi, 4)}.txt', format='%19.16e')
    print('done')
    
    # ----- calculating data for fig. 3 b) -------------------
    print('calculating data for fig. 3 b)...')
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
        'zinf': 3.7315,
        'zsup': 3.7356,
        'npts': 1000,
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
        F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)        
        calc.save_1D_data(omega, T, dir='data/fig3', filename=f'lmax={lmax}.txt', format='%19.16e')       
    print(f'CPU time:{time.time()-start}')


  
    # ----- calculating data for fig. 4 -------------------
if 'fig4' in figures_to_calculate:
    start = time.time()
    input_params = {
        'mts': '0', #multipole_type_selected
        'mos': '0', #multipole_order_selected
        'mps': '0', #m_projection_selected
        'type': '1 1',
        'order': '1 1',
        'm': '1 -1',
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
    F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)        
    calc.save_1D_data(omega, T, dir='data/fig4', filename='T.txt', format='%19.16e') 
    # with lapack and faddeeva
    F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)        
    calc.save_1D_data(omega, (T+R-1), dir='data/fig4', filename='error_with_lapack_and_faddeeva.txt', format='%19.16e') 
    # wo lapack
    input_params['multem_version'] = 'wo_lapack'
    F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)        
    calc.save_1D_data(omega, (T+R-1), dir='data/fig4', filename='error_wo_lapack.txt', format='%19.16e')       
    # with lapack
    input_params['multem_version'] = 'with_lapack'    
    F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)        
    calc.save_1D_data(omega, (T+R-1), dir='data/fig4', filename='error_with_lapack.txt', format='%19.16e') 
    print('done')
    print(f'CPU time:{time.time()-start}')


#takes approx 43000 seconds to calculate these data
if 'fig6' in figures_to_calculate:
    start = time.time()
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
    'fab': 90,
    'polar': 'S',
    'epssph_re': 5,
    'epssph_im': 0,
    'epsmed_re': 1,
    'epsmed_im': 0,
    'ktype': 2,
    'kscan': 1,
    'zinf': 0.5,
    'zsup': 1.0,
    'npts': 70,
    'r_ratio': 0.35,
    'mode': 'orig_multem', #multi-layered
    'dist_btw_spheres_and_interface': 0.5
}
    omega = np.linspace(1, 2, input_params['npts'])
    input_params['zinf'] = omega[0]
    input_params['zsup'] = omega[-1]
    ak1 = 0
    ak2 = 0
    # for input_params['multem_version'] in ['3']:
    #     for input_params['lmax'] in [7]:
    #         time0 = time.time()
    #         for input_params['rmax'] in [20, 25, 30]:
    #             F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)
    #             calc.save_1D_data(omega, T, dir=f'data/fig6/mode={input_params["mode"]}/version={input_params["multem_version"]}/d={input_params["dist_btw_spheres_and_interface"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}.txt', format='%19.16e') 
    #         print(f'{input_params["multem_version"]} : {time.time()-time0}s')
    #10738.055797576904s
    for input_params['multem_version'] in ['2']:
        for input_params['lmax'] in [7]:
            time0 = time.time()
            for input_params['rmax'] in [20, 25, 30]:
                F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)
                calc.save_1D_data(omega, T, dir=f'data/fig6/mode={input_params["mode"]}/version={input_params["multem_version"]}/d={input_params["dist_btw_spheres_and_interface"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}.txt', format='%19.16e') 
            print(f'{input_params["multem_version"]} : {time.time()-time0}s')
    for input_params['multem_version'] in ['2', '3']:
        for input_params['lmax'] in [7]:
            time0 = time.time()
            for input_params['rmax'] in [20, 30]:
                F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)
                calc.save_1D_data(omega, T, dir=f'data/fig6/mode={input_params["mode"]}/version={input_params["multem_version"]}/d={input_params["dist_btw_spheres_and_interface"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}.txt', format='%19.16e') 
            print(f'{input_params["multem_version"]} : {time.time()-time0}s')
    print(f'CPU time:{time.time()-start}')
   

if 'fig6_5' in figures_to_calculate:
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
        dir = f'data/fig6/'
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

#========================================================================
#fig 7
#calculation is implemented in dir Zaghoul_MIT_faddeeva_comparison
# fig 7 contains data calculated in fig3a, so one needs to calculate data for fig 3a first
# then run bash script cerf_faddeeva_comparison/run_CompPhysComm_data_prep.sh
#========================================================================

if 'fig8' in figures_to_calculate:
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
        
    def q_factor_estimate(hwhm, x0, q_factor_limit):
            is_enough_spectra = 0
            q = x0/(2*hwhm)
            if (q >= q_factor_limit):
                is_enough_spectra = 1


    def show_spectrum(figsize, x, y, sign):
        plt.figure(figsize=figsize)
        plt.plot(x, y, 'r', lw=0.5)
        plt.scatter(x, y)
        plt.title(sign)
        # plt.scatter([3.7355827121317358], [1])
        # 3.7355827103690302
        # 6.121325668573263e-11
        # plt.scatter([3.7355827103690302], [1])

        #only for example
        # plt.xlabel(r'${\omega d / 2\pi c }$')
        plt.ylabel('T')
        plt.show()

    def find_spectrum(ak1, x, th):
        ip = input_params
        from_x = x[0]
        to_x = x[-1]
        while True:
            F , T, R, A = calc.calc_spectrum_omega(x, ak1, ak2, input_params)
            y = T
            x = np.linspace(from_x, to_x, input_params['npts'])
            x_range = to_x - from_x
            show_spectrum((10,10), x, y, 'spectrum')
            hwhm, x0 = get_half_width_half_maxima_and_x0(x, y)
            hwhm_factor = 2
            from_x = x0 - hwhm_factor*hwhm
            to_x = x0 + hwhm_factor*hwhm
            if (hwhm/x_range*100 >= th): break
        F, T, R, A = calc.calc_spectrum_omega(x, ak1, ak2, input_params)   
        y = T
        x = np.linspace(from_x, to_x, input_params['npts'])
        hwhm, x0 = get_half_width_half_maxima_and_x0(x, y)
        is_enough_spectra = q_factor_estimate(hwhm, x0, 1e11)
        # show_spectrum((10,10), x, y, k_value)
        Q = x0/hwhm

        return x, y, is_enough_spectra, Q


    import time

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
    'r_ratio': 0.4705,
    'mode': '1', #multi-layered
    'dist_btw_spheres_and_interface': 0.6
}
    start = time.time()
    LMAX = [4]
    RMAX = [16]
    input_params['multem_version'] = '3'
    #fig8 b-e ?
    omega = np.linspace(3.734, 3.737, input_params['npts'])
    input_params['zinf'] = omega[0]
    input_params['zsup'] = omega[-1]
    ak1 = 2e-3/2/np.pi
    ak2 = 0
    for input_params['multem_version'] in ['3_cerf', '3']:
        for input_params['lmax'] in LMAX:
            time0 = time.time()
            for input_params['rmax'] in RMAX:
                F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)
                calc.save_1D_data(omega, T, dir=f'data/fig8/{input_params["multem_version"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}_ak1={round(ak1, 8)}.txt', format='%19.16e') 
            print(f'{input_params["multem_version"]} : {time.time()-time0}s')


    omega = np.linspace(3.7354, 3.7358, input_params['npts'])
    input_params['zinf'] = omega[0]
    input_params['zsup'] = omega[-1]
    ak1 = 1e-4/2/np.pi
    ak2 = 0
    for input_params['multem_version'] in ['3_cerf', '3']:
        for input_params['lmax'] in LMAX:
            time0 = time.time()
            for input_params['rmax'] in RMAX:
                F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)
                calc.save_1D_data(omega, T, dir=f'data/fig8/{input_params["multem_version"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}_ak1={round(ak1, 8)}.txt', format='%19.16e') 
            print(f'{input_params["multem_version"]} : {time.time()-time0}s')


    omega = np.linspace(3.73558, 3.735586, input_params['npts'])
    input_params['zinf'] = omega[0]
    input_params['zsup'] = omega[-1]
    ak1 = 5e-5/2/np.pi
    ak2 = 0
    for input_params['multem_version'] in ['3_cerf', '3']:
        for input_params['lmax'] in LMAX:
            time0 = time.time()
            for input_params['rmax'] in RMAX:
                F, T, R, A = calc.calc_spectrum_omega(omega, ak1, ak2, input_params)
                calc.save_1D_data(omega, T, dir=f'data/fig8/{input_params["multem_version"]}', filename=f'lmax={input_params["lmax"]}_rmax={input_params["rmax"]}_ak1={round(ak1, 8)}.txt', format='%19.16e') 
            print(f'{input_params["multem_version"]} : {time.time()-time0}s')
    print(f'CPU time: {time.time()-start}')



if 'fig9' in figures_to_calculate:
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
    dir = f'data/fig9/'
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



