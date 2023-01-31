import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt

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


def save_1D_data(x, y, dir, filename, format='%e'):
    if not os.path.exists(dir):
        os.makedirs(dir)
    data = np.stack((x, y), axis=-1)
    np.savetxt(dir+'/'+filename, data, delimiter='\t', fmt=format)
    print('data were saved')


# ----- calculating data for fig. 2 a) -------------------
print('calculating data for fig. 2 a)...')
ktype = 2
kscan = 1
fab = 60
rmax = 7
polar = 'S'
epssph_re = 15
epssph_im = 0
epsmed_re = 1
epsmed_im = 0
is_multipole_type_selected = '0'
is_multipole_order_selected = '0'
is_m_projection_selected = '0'
m_type = '0 0'
m_order = '1 1'
m = '-1 1'
npts = 1000
r_ratio = 0.47050000
zinf = 3.73
zsup = 3.74
ak1 = np.array((1e-2, 0, 5e-2, 1e-1))/2/np.pi
ak2 = np.array((0, 5e-2, 0, 0))/2/np.pi
kpts = len(ak1)
lmax = 4
x = np.linspace(zinf, zsup, npts)
# for i in range(kpts):
#     create_input(ktype, kscan, npts, fab, ak1[i], ak2[i], zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im,
#                  m_type, m_order, m, is_multipole_type_selected, is_multipole_order_selected,
#                  is_m_projection_selected, mode='1', multem_version='3')
#     F, T, R, A = eval(multem_version='3')
#     save_1D_data(x, T, dir='data/fig2a', filename=f'akxy={round(ak1[i]*2*np.pi, 4)}_{round(ak2[i]*2*np.pi, 4)}.txt', format='%19.16e')
# print('done')
# ----- calculating data for fig. 2 b) -------------------
# print('calculating data for fig. 2 b)...')
# zinf = 3.7315
# zsup = 3.7356
# npts = 2000
# ak1 = np.array((0.01))/2/np.pi
# ak2 = np.array((0.0))/2/np.pi
# LMAX = [4, 7, 10, 13]
# x = np.linspace(zinf, zsup, npts)
# for lmax in LMAX:
#     create_input(ktype, kscan, npts, fab, ak1, ak2, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im,
#                  m_type, m_order, m, is_multipole_type_selected, is_multipole_order_selected,
#                  is_m_projection_selected, mode='1', multem_version='3')
#     F, T, R, A = eval(multem_version='3')
#     save_1D_data(x, T, dir='data/fig2b', filename=f'lmax={lmax}.txt', format='%19.16e')
# print('done')
# ----- calculating data for fig. 2 c) -------------------
print('calculating data for fig. 2 c)...')
npts = 1000
zinf = 3.731
zsup = 3.735
x = np.linspace(zinf, zsup, npts)
ak1 = np.array((0.01))/2/np.pi
ak2 = np.array((0.0))/2/np.pi
lmax = 10
create_input(ktype, kscan, npts, fab, ak1, ak2, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im,
                 m_type, m_order, m, is_multipole_type_selected, is_multipole_order_selected,
                 is_m_projection_selected, mode='1', multem_version='2')
F, T, R, A = eval(multem_version='wo_lapack')
save_1D_data(x, (T+R-1), dir='data/fig2c', filename='error_wo_lapack.txt', format='%19.16e')
F, T, R, A = eval(multem_version='with_lapack')
save_1D_data(x, (T+R-1), dir='data/fig2c', filename='error_with_lapack.txt', format='%19.16e')
x0 = 3.732167206122217085e+00 #lmax10 #the value got from hand-crafting fano fitting code
x0_idx = np.argmin(np.abs(x - x0))
values_at_resonance_regions = (T + R - 1)[x0_idx-10:x0_idx+10]
print('average error with lapack:', np.mean(values_at_resonance_regions))
print('done')
# spectrum position offset between different lmax
# ktype = 2
# kscan = 1
# fab = 60
# rmax = 5
# polar = 'S'
# epssph_re = 15
# epssph_im = 0
# epsmed_re = 1
# epsmed_im = 0
# is_multipole_type_selected = '0'
# is_multipole_order_selected = '0'
# is_m_projection_selected = '0'
# m_type = '0 0'
# m_order = '1 1'
# m = '-1 1'
# npts = 100
# r_ratio = 0.47050000
# ak1 = np.array((1e-2))/2/np.pi
# ak2 = np.array((0.0))/2/np.pi
# # dict for lmax and approx resonance position: keys - lmax, values - approx resonance position
# X_approx = {4: [3.735, 3.7353],
#             5: [3.7332, 3.7334],
#             7: [3.7321, 3.7323],
#             8: [3.7321, 3.7322],
#             10: [3.7321, 3.7322],
#             11: [3.7321, 3.7322]}
# # dict of spectrum position: keys - lmax, values - resonance position
# X0 = {}
# for lmax, x_approx in X_approx.items():
#     print(f'Calculating lmax={lmax}')
#     from_k = x_approx[0]
#     to_k = x_approx[1]
#     zinf = from_k
#     zsup = to_k
#     create_input(ktype, kscan, npts, fab, ak1, ak2, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im,
#                  m_type, m_order, m, is_multipole_type_selected, is_multipole_order_selected,
#                  is_m_projection_selected, mode='1', multem_version='3')
#     T = eval(multem_version='3')[1]
#     x = np.linspace(from_k, to_k, npts)
#     save_1D_data(x, T, f'figures/Maksimov_Bulgakov/data/lmax_resonance_position/lmaxx={lmax}.txt', format='%19.16e')
#     # plt.plot(x, T)
#     # plt.show()
#     Tmax_idx = np.where(T == np.max(T)) # rude finding for resonance position
#     X0[lmax] = x[Tmax_idx]
#
# # finding x0 FWHM of resonances with fano fitting tool
# if os.path.isfile('fano-fit.py'): # check all dirs in source python file
#     print('Running fano-fit')
#     subprocess.call(['python3', 'fano-fit.py'])
#
# fitting_results_dir = 'figures/Maksimov_Bulgakov/data/lmax_resonance_position/_output'
# x0_after_fit = np.loadtxt(fitting_results_dir+'/'+'f0.txt') # resonance position after fitting
# x0_before_fit = np.array([float(x) for x in X0.values()]) #converting dict values to numpy array
# print('resonance postion finding error [%]', np.abs(x0_before_fit - x0_after_fit)/x0_before_fit * 100)
# Q = np.loadtxt(fitting_results_dir+'/'+'Q0s-fano.txt') # Q factor of the resonance
# fwhm_values = x0_after_fit/Q #calculating FWHM using definition of the Q factor
# # constructing 2 dicts:
# # 1. lmax and fwhm
# # 2. lmax and x0
# FWHM = {}
# x0 = {}
# for i, lmax in enumerate(X_approx.keys()):
#     FWHM[lmax] = fwhm_values[i]
#     x0[lmax] = x0_after_fit[i]
# # constructing dict for resonance position offset between different lmax: keys - lmax i.e#4-5#, values - offset in FWHM
# rp_offset = {'4-5': np.abs(x0[5] - x0[4])/FWHM[4],
#              '7-8': np.abs(x0[8] - x0[7])/FWHM[7],
#              '10-11': np.abs(x0[11] - x0[10])/FWHM[10]}
# print(rp_offset)
# --- End of the Section LMAX ------------------------------------------------------------------------

