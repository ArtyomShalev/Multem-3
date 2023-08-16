# TRA
# smuthi                    multem
#0.9953473462924776         0.998502
#0.004652653707521094       0.00149796
#1.3418086086680603e-15     8.23994e-18

import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
import visualisation_lib as vl


def create_input(ktype, kscan, npts, fab, ak1, ak2, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im, m_type, m_order, m,
                 is_multipole_type_selected, is_multipole_order_selected, is_m_projection_selected, mode, multem_version):
    '''
    mode:   1 - 2D array of spheres in homogeneous media
            2 - interface
    '''
    if multem_version == '2':
        float_format = '%13.8f'
    if multem_version == '3':
        float_format = '%19.15f'
    if mode == '1':
        str_fort10 = ('           ********************************************\n'
                      '           ********INPUT FILE FOR TRANSMISSION*********\n'
                      '           ********************************************\n'
                      '   KTYPE ='+'%2i'%(ktype)+'   KSCAN ='+'%2i'%(kscan)+'   KEMB  = 0    LMAX ='+'%2i'%(lmax)+'   NCOMP = 1   NUNIT = 1\n'
                      ' ALPHA =    1.000000  BETA =    1.000000   FAB =   '+'%9.6f'%(fab)+'  RMAX ='+'%11.6f'%(rmax)+'\n'
                      '  NP ='+'%4i'%(npts)+'  ZINF ='+
                      float_format%(zinf)+'  ZSUP ='+float_format%(zsup)+'\n'
                      '  THETA/AK(1) ='+float_format%(ak1)+'     FI/AK(2) ='+float_format%(ak2)+'   POLAR ='+polar+'     FEIN =   0.00\n'
                      '\n'
                      'Give information for the "NCOMP" components \n'
                      '\n'
                      '     IT  = 2\n'
                      '     MUMED =   1.00000000   0.00000000     EPSMED=   '+'%11.6f'%(epsmed_re)+'   '+'%11.6f'%(epsmed_im)+'\n'
                      '   NPLAN = 1  NLAYER = 1\n'
                      '       S =   '+'%10.8f'%(r_ratio)+'     MUSPH =   1.00000000   0.00000000     EPSSPH=  '+'%11.6f'%(epssph_re)+'   '+'%11.6f'%(epssph_im)+'\n'
                      'xyzDL 0.0  0.0  0.0\n'
                      'xyzDR 0.0  0.0  1.8\n'
                      '     MUEMBL=   1.00000000   0.00000000    EPSEMBL=   1.00000000   0.00000000\n'
                      '     MUEMBR=   1.00000000   0.00000000    EPSEMBR=   1.00000000   0.00000000\n')


        with open('fort.10','w') as f:
            print(str_fort10, file=f)

    if mode == '2':
        pass


    str_ini = ('[selectors]\n'
                   'is_multipole_type_selected = '+'%s'%(is_multipole_type_selected)+'\n'
                   'is_multipole_order_selected = '+'%s'%(is_multipole_order_selected)+'\n'
                   'is_m_projection_selected = '+'%s'%(is_m_projection_selected)+'\n'
                   '\n'
                   '[regime]\n'
                   'multipole_type = '+'%s'%(m_type)+'\n'
                   'multipole_order = '+'%s'%(m_order)+'\n'
                   'm_projection = '+'%s'%(m)+'\n')

    with open('multipole_regime_parameters.ini','w') as f:
        print(str_ini, file=f)


def eval(multem_version='3'):
    if multem_version == '3':
        if os.path.isfile('multem3'):
            my_env = os.environ.copy()
            my_env["OMP_NUM_THREADS"] = "1"
            # print("Running multem...")
            subprocess.run(['./multem3'],
                           stdout=subprocess.DEVNULL,
                           env=my_env)
    if multem_version == '2':
        if os.path.isfile('multem2'):
            my_env = os.environ.copy()
            my_env["OMP_NUM_THREADS"] = "1"
            subprocess.run(['./multem2'],
                           stdout=subprocess.DEVNULL,
                           env=my_env)
    if multem_version == 'wo_lapack':
        if os.path.isfile('multem2_wo_lapack'):
            print('wo_lapack')
            my_env = os.environ.copy()
            my_env["OMP_NUM_THREADS"] = "1"
            # print("Running multem...")
            subprocess.run(['./multem2_wo_lapack'],
                           stdout=subprocess.DEVNULL,
                           env=my_env)
    if multem_version == 'with_lapack':
        if os.path.isfile('multem2_with_lapack'):
            print('with_lapack')
            my_env = os.environ.copy()
            my_env["OMP_NUM_THREADS"] = "1"
            # print("Running multem...")
            subprocess.run(['./multem2_with_lapack'],
                           stdout=subprocess.DEVNULL,
                           env=my_env)
    if multem_version == 'zprint':
        if os.path.isfile('multem3_zprint'):
            my_env = os.environ.copy()
            my_env["OMP_NUM_THREADS"] = "1"
            # print("Running multem...")
            subprocess.run(['./multem3_zprint'],
                           stdout=subprocess.DEVNULL,
                           env=my_env)
    #FREQUENCY   TRANSMITTANCE  Reflectance   Absorbance
    F, T, R, A = np.loadtxt('fort.8').T
    return F, T, R, A


a = 350
wl1_dimless = a/0.46
wl2_dimless = a/0.445 
zinf = wl1_dimless/a
zsup = wl2_dimless/a
epsmed_re = 1
epsmed_im = 0
npts = 100
AK1 = np.arange(0, 90, 1)
transmittance = np.zeros((npts,len(AK1)))
reflectance = np.zeros((npts,len(AK1)))
absorbance = np.zeros((npts,len(AK1)))
# for idx, ak1 in enumerate(AK1):
#     print(idx)
#     create_input(ktype=1, kscan=2, npts=npts, fab=90, ak1=ak1, ak2=0, zinf=zinf, zsup=zsup, polar='S',
#                 lmax=3, r_ratio=100/a, rmax=7, epssph_re=50, epssph_im=0, m_type='1', m_order='3', m='0', 
#                 is_m_projection_selected=1, is_multipole_type_selected=1, is_multipole_order_selected=1, 
#                 mode='1', multem_version='3')

#     F, T, R, A = eval(multem_version='3')
#     transmittance[:,idx] = T
#     reflectance[:,idx] = R
#     absorbance[:,idx] = A

dir = '/home/ashalev/Projects/Multem-3/examples/smuthi_multem_comparison/magnetic/l=3/m=0'
# np.savetxt(f'{dir}/R_multem.txt', reflectance)

R_mul = np.loadtxt(f'{dir}/R_multem.txt')
# R_smu_0deg_azim = np.loadtxt(f'{dir}/R_smuthi_11row_0deg_azim.txt')
R_103_full_inc = np.loadtxt(f'{dir}/R_103_full_incident.txt')
R_103_full_inc_wide = np.loadtxt(f'{dir}/R_103_full_incident_widerange.txt')

R_smu_0deg_polar = np.loadtxt(f'{dir}/R_smuthi_11row_0deg_polar.txt')
R_smu_45deg_polar = np.loadtxt(f'{dir}/R_smuthi_11row_45deg_polar.txt')
R_smu_30deg_azim = np.loadtxt(f'{dir}/R_smuthi_11row_30deg_azim.txt')
# R_smu_random = np.loadtxt(f'{dir}/R.txt')


WL = np.linspace(wl1_dimless, wl2_dimless, npts)
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(2, 2, 1)
vl.plot_2D_map(AK1, WL, R_mul, logscale=True, vmin=1e-10, vmax=1, is_cb_needed=True, x_label='theta [deg]', y_label=r'$d/\lambda$', z_label='refl', extent=[0, 90, a/wl2_dimless, a/wl1_dimless])
ax.set_title('MULTEM M-103')

# ax = fig.add_subplot(2, 2, 2)
# vl.plot_2D_map(AK1, WL, R_smu_random, logscale=True, vmin=1e-30, vmax=1, x_label='theta [deg]??', y_label=r'$d/\lambda$', z_label='refl', extent=[0, 90, a/wl1_dimless, a/wl2_dimless])
# ax.set_title('SMUTHI M103 + randomTmatrix')

ax = fig.add_subplot(2, 2, 2)
vl.plot_2D_map(AK1, WL, R_103_full_inc, logscale=True, vmin=1e-45, vmax=1e-35, is_cb_needed=True, x_label='theta [deg]', y_label=r'$d/\lambda$', z_label='refl', extent=[0, 90, a/wl2_dimless, a/wl1_dimless])
ax.set_title('SMUTHI M-103 full incident field')


ax = fig.add_subplot(2, 2, 3)
# vl.plot_2D_map(AK1, WL, R_smu_0deg_polar, logscale=True, vmin=1e-30, vmax=1, x_label='theta [deg]??', y_label=r'$d/\lambda$', z_label='refl', extent=[0, 90, a/wl1_dimless, a/wl2_dimless])
vl.plot_2D_map(AK1, WL, R_103_full_inc_wide , logscale=True, vmin=1e-45, vmax=1e-35, is_cb_needed=True, x_label='theta [deg]', y_label=r'$d/\lambda$', z_label='refl', extent=[0, 90, 0.4, 0.6])
ax.set_title('SMUTHI M-103 full incident wide range')


# ax = fig.add_subplot(2, 2, 4)
# vl.plot_2D_map(AK1, WL, R_smu_45deg_polar, logscale=True, vmin=1e-30, vmax=1, x_label='theta [deg]??', y_label=r'$d/\lambda$', z_label='refl', extent=[0, 90, a/wl1_dimless, a/wl2_dimless])
# ax.set_title('SMUTHI M-103 polar tilted multipole 45 deg')

# ax = fig.add_subplot(2, 2, 4)
# vl.plot_2D_map(AK1, WL, R_smu_30deg_azim, logscale=True, vmin=1e-30, vmax=1, x_label='theta [deg]??', y_label=r'$d/\lambda$', z_label='refl', extent=[0, 90, a/wl1_dimless, a/wl2_dimless])
# ax.set_title('SMUTHI M-103 azim tilted multipole 30 deg')


plt.show()


