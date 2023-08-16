# ------------------------------------------------------------------------------------------------------------------------
#  this script calculated and draws figures from PHYSICAL REVIEW B 95, 195406 (2017)
#     Surface-lattice resonances in two-dimensional arrays of spheres: Multipolar interactions and a mode analysis
#     Sylvia D. Swiecicki and J. E. Sipe
#  using MULTEM and multipole decomposition modification
# ------------------------------------------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
import calculating_lib


# def create_input(ktype, kscan, npts, fab, ak1, ak2, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im, m_type, m_order, m,
#                  mts, mos, mps, mode):
#     '''
#     mode:   1 - 2D array of spheres in homogeneous media
#             2 - interface
#     '''
#
#     if mode == '1':
#         str_fort10 = ('           ********************************************\n'
#                       '           ********INPUT FILE FOR TRANSMISSION*********\n'
#                       '           ********************************************\n'
#                       '   KTYPE ='+'%2i'%(ktype)+'   KSCAN ='+'%2i'%(kscan)+'   KEMB  = 0    LMAX ='+'%2i'%(lmax)+'   NCOMP = 1   NUNIT = 1\n'
#                                                                                                                   ' ALPHA =    1.000000  BETA =    1.000000   FAB =   '+'%9.6f'%(fab)+'  RMAX ='+'%11.6f'%(rmax)+'\n'
#                                                                                                                                                                                                                  '  NP ='+'%4i'%(npts)+'  ZINF ='+
#                       '%19.15f'%(zinf)+'  ZSUP ='+'%19.15f'%(zsup)+'\n'
#                                                                    '  THETA/AK(1) ='+'%19.15f'%(ak1)+'     FI/AK(2) ='+'%19.15f'%(ak2)+'   POLAR ='+polar+'     FEIN =   0.00\n'
#                                                                                                                                                           '\n'
#                                                                                                                                                           'Give information for the "NCOMP" components \n'
#                                                                                                                                                           '\n'
#                                                                                                                                                           '     IT  = 2\n'
#                                                                                                                                                           '     MUMED =   1.00000000   0.00000000     EPSMED=   1.0000000   0.00000000\n'
#                                                                                                                                                           '   NPLAN = 1  NLAYER = 1\n'
#                                                                                                                                                           '       S =   '+'%10.8f'%(r_ratio)+'     MUSPH =   1.00000000   0.00000000     EPSSPH=  '+'%11.6f'%(epssph_re)+'   '+'%11.8f'%(epssph_im)+'\n'
#                                                                                                                                                                                                                                                                                                     'xyzDL 0.0  0.0  0.0\n'
#                                                                                                                                                                                                                                                                                                     'xyzDR 0.0  0.0  1.8\n'
#                                                                                                                                                                                                                                                                                                     '     MUEMBL=   1.00000000   0.00000000    EPSEMBL=   1.00000000   0.00000000\n'
#                                                                                                                                                                                                                                                                                                     '     MUEMBR=   1.00000000   0.00000000    EPSEMBR=   1.00000000   0.00000000\n')
#
#
#
#
#     if mode == '2':
#         str_fort10 = ('           ********************************************\n'
#                       '           ********INPUT FILE FOR TRANSMISSION*********\n'
#                       '           ********************************************\n'
#                       '   KTYPE ='+'%2i'%(ktype)+'   KSCAN ='+'%2i'%(kscan)+'   KEMB  = 0    LMAX ='+'%2i'%(lmax)+'   NCOMP = 2   NUNIT = 1\n'
#                                                                                                                   ' ALPHA =    1.000000  BETA =    1.000000   FAB =   '+'%9.6f'%(fab)+'  RMAX ='+'%11.6f'%(rmax)+'\n'
#                                                                                                                                                                                                                  '  NP ='+'%4i'%(npts)+'  ZINF ='+
#                       '%19.15f'%(zinf)+'  ZSUP ='+'%19.15f'%(zsup)+'\n'
#                                                                    '  THETA/AK(1) ='+'%19.15f'%(ak1)+'     FI/AK(2) ='+'%19.15f'%(ak2)+'   POLAR ='+polar+'     FEIN =   0.00\n'
#                                                                                                                                                           '\n'
#                                                                                                                                                           'Give information for the "NCOMP" components \n'
#                                                                                                                                                           '\n'
#                                                                                                                                                           '     IT  = 2\n'
#                                                                                                                                                           '     MUMED =   1.00000000   0.00000000     EPSMED=   1.00000000   0.00000000\n'
#                                                                                                                                                           '   NPLAN = 1  NLAYER = 1\n'
#                                                                                                                                                           '       S =   '+'%10.8f'%(r_ratio)+'     MUSPH =   1.00000000   0.00000000     EPSSPH=  '+'%11.6f'%(epssph_re)+'   '+'%11.8f'%(epssph_im)+'\n'
#                                                                                                                                                                                                                                                                                                     'xyzDL 0.0  0.0  0.0\n'
#                                                                                                                                                                                                                                                                                                     'xyzDR 0.0  0.0  1.8\n'
#                                                                                                                                                                                                                                                                                                     '     IT  = 1\n'
#                                                                                                                                                                                                                                                                                                     '  DSLAB       =  1.000000000000000\n'
#                                                                                                                                                                                                                                                                                                     '     MU1   =   1.00000000   0.00000000     EPS1  =   1.0000000   0.00000000\n'
#                                                                                                                                                                                                                                                                                                     '     MU2   =   1.00000000   0.00000000     EPS2  =  '+'%11.6f'%(epssph_re)+'   '+'%11.8f'%(epssph_im)+'\n'
#                                                                                                                                                                                                                                                                                                                                                                                                            '     MU3   =   1.00000000   0.00000000     EPS3  =   1.00000000   0.00000000\n'
#                                                                                                                                                                                                                                                                                                                                                                                                            'xyzDL 0.0  0.0  0.0\n'
#                                                                                                                                                                                                                                                                                                                                                                                                            'xyzDR 0.0  0.0  0.0\n')
#
#
#
#
#     with open('fort.10','w') as f:
#         print(str_fort10, file=f)
#
#     str_ini = ('[selectors]\n'
#                'is_multipole_type_selected = '+'%s'%(mts)+'\n'
#                                                           'is_multipole_order_selected = '+'%s'%(mos)+'\n'
#                                                                                                       'is_m_projection_selected = '+'%s'%(mps)+'\n'
#                                                                                                                                                '\n'
#                                                                                                                                                '[regime]\n'
#                                                                                                                                                'multipole_type = '+'%s'%(m_type)+'\n'
#                                                                                                                                                                                  'multipole_order = '+'%s'%(m_order)+'\n'
#                                                                                                                                                                                                                      'm_projection = '+'%s'%(m)+'\n')
#
#     with open('multipole_regime_parameters.ini','w') as f:
#         print(str_ini, file=f)
#
#
#
#
# def eval(): #TODO refactor
#     # create_input(ktype, npts, ap1, ap2, from_y, to_y, polar, lmax, r_ratio, rmax, epssph_re, epssph_im,
#     #              type, order, is_multipole_type_selected, is_multipole_order_selected, is_m_projection_selected)
#     # print("Trying to run multem...")
#     if os.path.isfile('multem3'):
#         my_env = os.environ.copy()
#         my_env["OMP_NUM_THREADS"] = "1"
#         # print("Running multem...")
#         subprocess.run(['./multem3'],
#                        stdout=subprocess.DEVNULL,
#                        env=my_env)
#         # print("Done.")
#     #FREQUENCY   TRANSMITTANCE  Reflectance   Absorbance
#     F, T, R, A = np.loadtxt('fort.8').T
#     return F, T, R, A

if __name__ == "__main__":

    plt.rcParams.update({'font.size': 26, 'font.serif':"Times New Roman"})

    input_parameters = {
        'ktype': 1,
        'lmax': 3,
        'rmax': 20,
        'lattice_constant': 475,
        'radius': 100,
        'epssph_re': -20.1480000,
        'epssph_im': 1.24700000,
        'epsmed_re': 2.1025,
        'polar': 'S',
        'fab': 60,
        'mts': '0',
        'mos': '0',
        'mps': '0',
        'type': '1',
        'order': '1',
        'm': '1'
    }
    input_parameters['r_ratio'] = input_parameters['radius']/input_parameters['lattice_constant']

    #fig 4
    fi = 0
    from_sin_theta = 0.0
    to_sin_theta = 0.999
    n_theta = 1000
    sin_theta = np.linspace(from_sin_theta, to_sin_theta, n_theta)
    theta = np.arcsin(sin_theta)*180/np.pi

    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(1, 1, 1)
    for input_parameters['rmax'] in [2, 5, 7, 10, 16, 20]:
        dir = f'PRB/fig4/{input_parameters["rmax"]}'
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
            R = calculating_lib.calc_spectrum_ak1(ktype=1, ak1=theta, wl=wl, ip=input_parameters)
            calculating_lib.save_result(dir, f'{wl}', R)
            # ax1.plot(sin_theta, R+factor)
            factor += 1
        calculating_lib.save_result(dir, 'sintheta', sin_theta)
        # plt.show()



    # # fig 3 a
    # ax1 = fig.add_subplot(1, 2, 1)
    # ktype = 1
    # lmax = 3
    # rmax = 5
    #
    # #  figure 7 a
    #
    # kpts = n_theta
    # data_arr = np.empty((kpts, 4))
    # a = 475.0
    # s = 100.0
    # r_ratio = s/a
    # lambda_incident = 750
    # epssph_re = -20.1480000
    # epssph_im = 1.24700000
    # zinf = lambda_incident/a
    # zsup = (lambda_incident+0.01)/a
    # npts = 2
    # polar='S' # S or P
    #
    # # ps multipole
    # is_multipole_type_selected = '1'
    # is_multipole_order_selected = '1'
    # is_m_projection_selected = '1'
    # type = '0 0'
    # order = '1 1'
    # m = '-1 1'
    #
    # for i in range(kpts):
    #     eval(i, zinf, zsup, epssph_re, epssph_im, rmax)
    #
    # ax1.plot(x, data_arr[:, 2], 'r--', label='ps', lw = 2.0)
    #
    # # ps + mz + qks
    # is_multipole_type_selected = '1'
    # is_multipole_order_selected = '1'
    # is_m_projection_selected = '1'
    # type =  '0 0 1 0 0'
    # order = '1 1 1 2 2'
    # m =     '-1 1 0 -2 2'
    #
    # for i in range(kpts):
    #     eval(i, zinf, zsup, epssph_re, epssph_im, rmax)
    #
    # ax1.plot(x, data_arr[:, 2], 'r', label='ps+mz+qks', lw = 2.0)
    #
    # # plot config
    # ax1.set_title(f'{lambda_incident} nm')
    # ax1.set_ylabel('R')
    # ax1.set_xlabel(r'sin$\theta$')
    # ax1.set_ylim(-0.01, 1.01)
    # # ax1.legend()
    #
    # #  figure 7b
    # ax2 = fig.add_subplot(1, 2, 2)
    # fi = 0
    # from_sin_theta = 0.4
    # to_sin_theta = 0.8
    # n_theta = 100
    # sin_theta = np.linspace(from_sin_theta, to_sin_theta, n_theta)
    # x = sin_theta
    # theta = np.arcsin(sin_theta) * 180 / np.pi
    # kpts = n_theta
    # data_arr = np.empty((kpts, 4))
    # a = 475.0
    # s = 100.0
    # r_ratio = s / a
    # lambda_incident = 900
    # epssph_re = -20.1480000
    # epssph_im = 1.24700000
    # zinf = lambda_incident / a
    # zsup = (lambda_incident + 0.01) / a
    # npts = 2
    # polar = 'S'  # S or P
    #
    # # ps multipole
    # is_multipole_type_selected = '1'
    # is_multipole_order_selected = '1'
    # is_m_projection_selected = '1'
    # type = '0 0'
    # order = '1 1'
    # m = '-1 1'
    #
    # for i in range(kpts):
    #     eval(i, zinf, zsup, epssph_re, epssph_im, rmax)
    #
    # ax2.plot(x, data_arr[:, 2], 'r--', label='ps', lw=2.0)
    #
    # # ps + mz + qks
    # is_multipole_type_selected = '1'
    # is_multipole_order_selected = '1'
    # is_m_projection_selected = '1'
    # type = '0 0 1 0 0'
    # order = '1 1 1 2 2'
    # m = '-1 1 0 -2 2'
    #
    # for i in range(kpts):
    #     eval(i, zinf, zsup, epssph_re, epssph_im, rmax)
    #
    # ax2.plot(x, data_arr[:, 2], 'r', label='ps+mz+qks', lw=2.0)
    #
    # # plot config
    # ax2.set_title(f'{lambda_incident} nm')
    # # ax2.set_ylabel('R')
    # ax2.set_xlabel(r'sin$\theta$')
    # # ax2.set_ylim(-0.01, 1.01)
    # ax2.set_yticklabels([])
    # ax2.legend()
    #
    # # plt.subplot_tool()
    # plt.subplots_adjust(left=0.068,
    #                     bottom=0.11,
    #                     right=0.978,
    #                     top=0.934,
    #                     wspace=0.064,
    #                     hspace=0.2)
    #
    # # plt.show()
    # plt.savefig('/home/ashalev/Projects/amos-try/multem2mod/multem3article/figures/ready_to_publish/fig4.pdf')
    # plt.clf(); plt.close()