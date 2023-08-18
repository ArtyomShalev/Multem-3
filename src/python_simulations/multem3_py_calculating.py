import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import math as m
import cmath
import time
import matplotlib.colors as mcolors
from matplotlib.colors import hsv_to_rgb
from colorsys import hls_to_rgb
# from typing import str


def multi_loadtxt(dir: str, filelist: list):
    output = ()
    for fname in filelist:
        out = np.loadtxt(dir+"/"+fname)
        output += (out,)
    #TODO what are types for filelist and output
    return output


def create_input(ktype, kscan, npts, fab, ak1, ak2, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, 
                 epssph_im, epsmed_re, epsmed_im, m_type, m_order, m, mts, mos, mps, mode, multem_version, dist_btw_spheres_and_interface=1):
    '''
    mode:   1 - 2D array of spheres in homogeneous media
            2 - interface
            3 - multilayered-design
    '''
    if multem_version == '3':
        float_format = '%19.15f'
    else:
        float_format = '%13.8f'

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

    if mode == '2':
        str_fort10 = ('           ********************************************\n'
                      '           ********INPUT FILE FOR TRANSMISSION*********\n'
                      '           ********************************************\n'
                      '   KTYPE ='+'%2i'%(ktype)+'   KSCAN ='+'%2i'%(kscan)+'   KEMB  = 0    LMAX ='+'%2i'%(lmax)+'   NCOMP = 2   NUNIT = 1\n'
                                                                                                                  ' ALPHA =    1.000000  BETA =    1.000000   FAB =   '+'%9.6f'%(fab)+'  RMAX ='+'%11.6f'%(rmax)+'\n'
                                                                                                                                                                                                                 '  NP ='+'%4i'%(npts)+'  ZINF ='+
                      float_format%(zinf)+'  ZSUP ='+float_format%(zsup)+'\n'
                                                                   '  THETA/AK(1) ='+float_format%(ak1)+'     FI/AK(2) ='+float_format%(ak2)+'   POLAR ='+polar+'     FEIN =   0.00\n'
                    '\n'
                    'Give information for the "NCOMP" components \n'
                    '\n'
                    '     IT  = 2\n'
                    '     MUMED =   1.00000000   0.00000000     EPSMED=   1.00000000   0.00000000\n'
                    '   NPLAN = 1  NLAYER = 1\n'
                    '       S =   '+'%10.8f'%(r_ratio)+'     MUSPH =   1.00000000   0.00000000     EPSSPH=  '+'%11.6f'%(epssph_re)+'   '+'%11.8f'%(epssph_im)+'\n'
                    'xyzDL 0.0  0.0  0.0\n'
                    'xyzDR 0.0  0.0  '+'%2.1f'%(dist_btw_spheres_and_interface)+'\n'
                    '     IT  = 1\n'
                    '  DSLAB       =  1.000000000000000\n'
                    '     MU1   =   1.00000000   0.00000000     EPS1  =   1.0000000   0.00000000\n'
                    '     MU2   =   1.00000000   0.00000000     EPS2  =  '+'%11.6f'%(epssph_re)+'   '+'%11.8f'%(epssph_im)+'\n'
                    '     MU3   =   1.00000000   0.00000000     EPS3  =   1.00000000   0.00000000\n'
                    'xyzDL 0.0  0.0  0.0\n'
                    'xyzDR 0.0  0.0  0.0\n')
        
    if mode == '3':
        str_fort10 = ('           ********************************************\n'
                    '           ********INPUT FILE FOR TRANSMISSION*********\n'
                    '           ********************************************\n'
                    '   KTYPE ='+'%2i'%(ktype)+'   KSCAN ='+'%2i'%(kscan)+'   KEMB  = 0    LMAX ='+'%2i'%(lmax)+'   NCOMP = 2   NUNIT = 1\n'
                                                                                                                ' ALPHA =    1.000000  BETA =    1.000000   FAB =   '+'%9.6f'%(fab)+'  RMAX ='+'%11.6f'%(rmax)+'\n'
                                                                                                                                                                                                                '  NP ='+'%4i'%(npts)+'  ZINF ='+
                    float_format%(zinf)+'  ZSUP ='+float_format%(zsup)+'\n'
                                                                '  THETA/AK(1) ='+float_format%(ak1)+'     FI/AK(2) ='+float_format%(ak2)+'   POLAR ='+polar+'     FEIN =   0.00\n'
                '\n'
                'Give information for the "NCOMP" components \n'
                '\n'
                '     IT  = 2\n'
                '     MUMED =   1.00000000   0.00000000     EPSMED=   1.00000000   0.00000000\n'
                '   NPLAN = 1  NLAYER = 1\n'
                '       S =   '+'%10.8f'%(r_ratio)+'     MUSPH =   1.00000000   0.00000000     EPSSPH=  '+'%11.6f'%(epssph_re)+'   '+'%11.8f'%(epssph_im)+'\n'
                'xyzDL 0.0  0.0  0.0\n'
                'xyzDR 0.0  0.0  '+'%2.1f'%(dist_btw_spheres_and_interface)+'\n'
                '     IT  = 2\n'
                '     MUMED =   1.00000000   0.00000000     EPSMED=   1.00000000   0.00000000\n'
                '   NPLAN = 1  NLAYER = 1\n'
                '       S =   '+'%10.8f'%(r_ratio)+'     MUSPH =   1.00000000   0.00000000     EPSSPH=  '+'%11.6f'%(epssph_re)+'   '+'%11.8f'%(epssph_im)+'\n'
                'xyzDL 0.0  0.0  0.0\n'
                'xyzDR 0.0  0.0  0.0\n')

    with open('fort.10','w') as f:
        print(str_fort10, file=f)

    str_ini = ('[selectors]\n'
               'is_multipole_type_selected = '+'%s'%(mts)+'\n'
                'is_multipole_order_selected = '+'%s'%(mos)+'\n'
                'is_m_projection_selected = '+'%s'%(mps)+'\n'
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
    if multem_version == '3_cerf':
        if os.path.isfile('multem3_cerf'):
            my_env = os.environ.copy()
            my_env["OMP_NUM_THREADS"] = "1"
            # print("Running multem...")
            subprocess.run(['./multem3_cerf'],
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



# def calc_2Dmap_from_multem(x, y, zinf, zsup, znpts, is_w0_needed=False):
#     xnpts, ynpts = len(x), len(y)
#     T, R = np.empty((xnpts, ynpts, znpts)), np.empty((xnpts, ynpts, znpts))
#     if is_w0_needed:
#         W, THETA = np.empty((xnpts, ynpts), dtype=np.complex), np.empty((xnpts, ynpts))
#     calc_counter = 0
#     for i in range(xnpts):
#         for j in range(ynpts):
#             create_input(ktype, kscan, znpts, fab, ak1)
#             f, t, r, a = eval(x[i], y[j], zinf, zsup, znpts)
#             if is_w0_needed:
#                 W[i,j], THETA[i,j] = get_W0(zinf, zsup, znpts, polar)
#             T[i, j, :], R[i, j, :] = t, r
#             calc_counter += 1
#             print(f'{calc_counter} from {xnpts*ynpts}')
#     if is_w0_needed:
#         return T, R, W, THETA
#     else:
#         return T, R


def calc_kR_map(k, R, ip):
    #calculating kR map like in PHYSICAL REVIEW B 105, 075404 (2022) fig 1c-d
    #ip - input params dictionary
    klen = len(k)
    Rlen = len(R)
    value = np.empty((klen, Rlen))
    for i in range(Rlen):
        print(f'{i} from {Rlen}')
        create_input(ktype=1, kscan=1, npts=klen, fab=ip['fab'], ak1=0, ak2=0, zinf=k[0], zsup=k[-1],
                     polar=ip['polar'], lmax=ip['lmax'], r_ratio=R[i], rmax=ip['rmax'], epssph_re=ip['epssph_re'],
                     epssph_im=ip['epssph_im'], m_type=ip['type'], m_order=ip['order'],
                     m=ip['m'], mts=ip['mts'], mos=ip['mos'], mps=ip['mps'], mode='1')
        f, t, r, a = eval()
        value[:, i] = t #choose reflectance, transmittance or absorbance
    return value


def calc_kxky_map(kx, ky, ip, is_w0_needed):
    #ip - input parameters
    #w0 - quasimodal coupling strength
    kxlen = len(kx); kylen = len(ky)
    T = np.empty((kxlen, kylen))
    R = np.empty((kxlen, kylen))
    if is_w0_needed:
        W, THETA = np.empty((kxlen, kylen), dtype=np.complex), np.empty((kxlen, kylen))
    counter = 0
    for i in range(kxlen):
        for j in range(kylen):
            print(f'{counter} from {kxlen*kylen}')
            counter += 1
            create_input(ktype=2, kscan=1, npts=ip['npts'], fab=ip['fab'], ak1=kx[i], ak2=ky[j],
                         zinf=ip['zinf'], zsup=ip['zsup'],
                         polar=ip['polar'], lmax=ip['lmax'], r_ratio=ip['r_ratio'],
                         rmax=ip['rmax'], epssph_re=ip['epssph_re'],
                         epssph_im=ip['epssph_im'], m_type=ip['type'], m_order=ip['order'],
                         m=ip['m'], mts=ip['mts'], mos=ip['mos'], mps=ip['mps'], mode='1')
            f, t, r, a = eval()
            T[i, j] = t[0] #choose reflectance, transmittance or absorbance
            R[i, j] = r[0] #choose reflectance, transmittance or absorbance
            if is_w0_needed:
                W[i,j], THETA[i,j] = get_W0(ip['zinf'], ip['zsup'], ip['npts'], ip['polar'])
    if is_w0_needed:
        return T, R, W, THETA
    else:
        return T, R


def calc_ak1_omega_map(ktype, ak1, omega, ip, is_w0_needed):
    #ktype: 1 - ak1 - theta; 2 - ak1 - kx (2pi/alpha)
    ak1len = len(ak1)
    omegalen=len(omega)
    zinf = omega[0]*2*np.pi
    zsup = omega[-1]*2*np.pi
    T = np.empty((ak1len, omegalen))
    R = np.empty((ak1len, omegalen))
    if is_w0_needed:
        W, THETA = np.empty((ak1len, omegalen), dtype=np.complex), np.empty((ak1len, omegalen))
    for i in range(ak1len):
        print(f'{i} from {ak1len}')
        create_input(ktype=ktype, kscan=1, npts=omegalen, fab=ip['fab'], ak1=ak1[i], ak2=0,
                     zinf=zinf, zsup=zsup,
                     polar=ip['polar'], lmax=ip['lmax'], r_ratio=ip['r_ratio'],
                     rmax=ip['rmax'], epssph_re=ip['epssph_re'],
                     epssph_im=ip['epssph_im'], m_type=ip['type'], m_order=ip['order'],
                     m=ip['m'], mts=ip['mts'], mos=ip['mos'], mps=ip['mps'], mode=['mode'])
        f, t, r, a = eval()
        T[i, :] = t #choose reflectance, transmittance or absorbance
        R[i, :] = r #choose reflectance, transmittance or absorbance
        if is_w0_needed:
            W[i, :], THETA[i, :] = get_W0(zinf, zsup, omegalen, ip['polar'])
    if is_w0_needed:
        return T, R, W, THETA
    else:
        return T, R


def calc_spectrum_omega(omega, ak1, ak2, ip):
    zinf = omega[0]
    zsup = omega[-1]
    omegalen=len(omega)
    if ip['mode'] == '1':
        dist_btw_spheres_and_interface = 0
    else:
        dist_btw_spheres_and_interface = ip['dist_btw_spheres_and_interface']
    create_input(ktype=ip['ktype'], kscan=1, npts=omegalen, fab=ip['fab'], ak1=ak1, ak2=ak2,
                 zinf=zinf, zsup=zsup,
                 polar=ip['polar'], lmax=ip['lmax'], r_ratio=ip['r_ratio'],
                 rmax=ip['rmax'], epssph_re=ip['epssph_re'],
                 epssph_im=ip['epssph_im'], epsmed_re=ip['epsmed_re'], epsmed_im=ip['epsmed_im'], 
                 m_type=ip['type'], m_order=ip['order'],
                 m=ip['m'], mts=ip['mts'], mos=ip['mos'], mps=ip['mps'], mode=ip['mode'], multem_version=ip['multem_version'], dist_btw_spheres_and_interface=dist_btw_spheres_and_interface)
    F, T, R, A = eval(ip['multem_version'])
    return F, T, R, A


def calc_spectrum_ak1(ktype, ak1, wl, ip):
    #wl: wavelengths in nm
    zinf = wl/ip['lattice_constant']
    zsup = zinf + 1e-6
    ak1len = len(ak1)
    R = np.empty((ak1len))
    for i in range(ak1len):
        print(f'{i} from {ak1len}')
        create_input(ktype=ktype, kscan=2, npts=2, fab=ip['fab'], ak1=ak1[i], ak2=0,
                     zinf=zinf, zsup=zsup,
                     polar=ip['polar'], lmax=ip['lmax'], r_ratio=ip['r_ratio'],
                     rmax=ip['rmax'], epssph_re=ip['epssph_re'],
                     epssph_im=ip['epssph_im'], epsmed_re=ip['epsmed_re'], epsmed_im=ip['epsmed_im'],
                     m_type=ip['type'], m_order=ip['order'],
                     m=ip['m'], mts=ip['mts'], mos=ip['mos'], mps=ip['mps'], mode='1', multem_version=ip['multem_version'])
        f, t, r, a = eval()
        R[i] = r[0]
    return R


def show_2D_map(x, y, z, zmin, zmax, logscale=False):
    if logscale:
        im = plt.imshow(z.T, extent = (np.min(x), np.max(x), np.min(y), np.max(y)), cmap=cm.rainbow, norm=LogNorm(vmin=zmin, vmax=1), aspect='auto', interpolation = 'none', origin='lower')
    else:
        im = plt.imshow(z.T, extent = (np.min(x), np.max(x), np.min(y), np.max(y)), cmap=cm.rainbow, vmin=zmin, vmax=zmax, aspect='auto', interpolation = 'none', origin='lower')

    cb = plt.colorbar(im)
    plt.show()


def get_W0(from_y, to_y, npts, polar):
    qiii = np.loadtxt('fort.1')
    len_per_freq = int(len(qiii)/npts)
    qiii_dict = {f'{from_y}': qiii[:len_per_freq],
                 f'{to_y}': qiii[len_per_freq:]}
    L = []
    omega = to_y
    for el in qiii_dict[str(omega)]:
        el = float(el[0]) + 1j*float(el[1])
        L.append(el)
    L = np.array(L).reshape(int(np.sqrt(len_per_freq)),int(np.sqrt(len_per_freq)))
    if polar == 'S':
        W0 = L[1][1]
    elif polar == 'P':
        W0 = L[0][0]
    else:
        raise ValueError
    r = abs(W0)
    phi = cmath.phase(W0)+0*np.pi/2

    W0_rect_form = cmath.rect(r, phi)
    theta = cmath.phase(W0)
    theta = cmath.phase(W0_rect_form)
    return W0_rect_form, theta


def save_result(path, filename, data):
    if not (os.path.exists(path)):
        os.makedirs(path)
    np.savetxt(path+'/'+filename+'.txt', data)
    print('results were saved')


def save_1D_data(x, y, dir, filename, format='%e'):
    if not os.path.exists(dir):
        os.makedirs(dir)
    data = np.stack((x, y), axis=-1)
    np.savetxt(dir+'/'+filename, data, delimiter='\t', fmt=format)
    print('data were saved')

def read_1D_data_real(filename):
    X, Y = [], []
    with open(filename) as data:
        for line in data:
            x, y = line.split()
            X.append(float(x))
            Y.append(float(y))
    X = np.array(X); Y = np.array(Y)
    return X, Y


def read_1D_data_complex(filename):
    X_RE, Y_RE = [], []
    X_IM, Y_IM = [], []
    with open(filename) as data:
        for line in data:
            x_re, x_im, y_re, y_im = line.split()
            X_RE.append(float(x_re))
            X_IM.append(float(x_im))
            Y_RE.append(float(y_re))
            Y_IM.append(float(y_im))
    X_RE = np.array(X_RE); X_IM = np.array(X_IM)
    Y_RE = np.array(Y_RE); Y_IM = np.array(Y_IM)
    return X_RE, X_IM, Y_RE, Y_IM


if __name__ == '__main__':
    # is_multipole_type_selected = '1'
    # is_multipole_order_selected = '1'
    # is_m_projection_selected = '1'
    # type =   '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1'
    # order =  '1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4'
    # m =      '-1 0 1 -2 -1 0 1 2 -3 -2 -1 0 1 2 3 -4 -3 -2 -1 0 1 2 3 4 -1 0 1 -2 -1 0 1 2 -3 -2 -1 0 1 2 3 -4 -3 -2 -1 0 1 2 3 4'
    # type =   '1 1'
    # order =  '1 3'
    # m =      '0 0'

    # rmax = 10
    # lmax = 4
    # a = 350
    # s = 100
    # r_ratio = s/a
    # from_y = 0.4
    # to_y = 0.5
    # ny = 300
    # y = np.linspace(from_y, to_y, ny)
    # zinf = from_y*2*np.pi
    # zsup = to_y*2*np.pi
    # nx = 100
    # ktype = 2
    # kx = np.linspace(0.0/2, 0.5/2, nx)
    # ky = np.linspace(-0.3/2, 0.3/2, ny)
    # ky = [0]
    # polar='S' # S or P
    # epssph_re = 50.0
    # epssph_im = 0.0000
    # # T, R, W, THETA = calc_2Dmap_from_multem(kx, ky, zinf, zsup, 2, is_w0_needed=True)
    # # R = R[:, :, 0]
    # T, R = calc_2Dmap_from_multem(kx, ky, zinf, zsup, ny)
    # R = R[:, 0, :]
    # show_2D_map(kx*2, y, R, zmin=1e-10)
    # path = f'/home/ashalev/Projects/bic_visualisation/topological_charge3/lmax4_wo_M101_M103_zoomed_kx0-01/{polar}'
    # save_result(path, 'R', R)
    # save_result(path, 'kx', kx*2)
    # # save_result(path, 'ky', ky*2)
    # # save_result(path, 'w_re', W.real)
    # # save_result(path, 'w_im', W.imag)
    # save_result(path, 'dim-less_omega', y)

    # nx = 30
    # ny = 30
    # kx = np.linspace(-0.3/2, 0.3/2, nx)
    # ky = np.linspace(-0.3/2, 0.3/2, ny)
    # T, R, W, THETA = calc_2Dmap_from_multem(kx, ky, zinf, zsup, 2, is_w0_needed=True)
    # save_result(path, 'kx', kx*2)
    # save_result(path, 'ky', ky*2)
    # save_result(path, 'theta', THETA)





    # BIC ring PRB 2022
    #normal incidence
    # input_params = {
    #     'mts': '0', #multipole_type_selected
    #     'mos': '0', #multipole_order_selected
    #     'mps': '0', #m_projection_selected
    #     'type': '1',
    #     'order': '1',
    #     'm': '1',
    #     'rmax': 5,
    #     'lmax': 8,
    #     'lattice_constant': 300,
    #     'fab': 60,
    #     'polar': 'S',
    #     'epssph_re': 15,
    #     'epssph_im': 0,
    #     'ktype': 1,
    #     'kscan': 1,
    # }
    # from_R = 0.37; to_R = 0.38; n_R = 500
    # R = np.linspace(from_R, to_R, n_R)
    # from_k = 2.95; to_k = 3.1; n_k = 500
    # k = np.linspace(from_k, to_k, n_k)
    # T = calc_kR_map(k, R, input_params)
    # show_2D_map(k, R, T, zmin=0, zmax=1)
    # dir = 'kxky/0.378/P'
    # save_result(dir, 'T', T)
    # save_result(dir, 'R', R)
    # save_result(dir, 'k', k)
    # kx = np.loadtxt(dir+'/'+'kx.txt')
    # ky = np.loadtxt(dir+'/'+'ky.txt')
    # T = np.loadtxt(dir+'/'+'T.txt')
    # R = np.loadtxt(dir+'/'+'R.txt')

    # plt.plot(kx, 0.3764, 'o', color='white')
    # show_2D_map(kx, ky, R, zmin=1e-5, zmax=1, logscale=True)


    # kxky maps
    # input_params = {
    #     'mts': '0', #multipole_type_selected
    #     'mos': '0', #multipole_order_selected
    #     'mps': '0', #m_projection_selected
    #     'type': '1',
    #     'order': '1',
    #     'm': '1',
    #     'rmax': 2,
    #     'lmax': 8,
    #     'lattice_constant': 300,
    #     'fab': 60,
    #     'polar': 'S',
    #     'epssph_re': 15,
    #     'epssph_im': 0,
    #     'ktype': 1,
    #     'kscan': 1,
    #     'zinf': 2.9829,
    #     'zsup': 2.983,
    #     'npts': 2,
    #     'r_ratio': 0.3764
    # }
    # span = 0.4
    # n = 500
    # kx = np.linspace(-span, span, n)
    # ky = np.linspace(-span, span, n)
    # R_RATIO = [0.374]
    # for input_params['r_ratio'] in R_RATIO:
    #     T, R, W, THETA = calc_kxky_map(kx, ky, input_params, is_w0_needed=True)
    #     show_2D_map(kx, ky, R, 1e-5, 1, logscale=True)
    #     dir = f'kxky/{input_params["r_ratio"]}/{input_params["polar"]}'
    #     save_result(dir, 'kx', kx)
    #     save_result(dir, 'ky', ky)
    #     save_result(dir, 'T', T)
    #     save_result(dir, 'R', R)
    #     save_result(dir, 'w_re', W.real)
    #     save_result(dir, 'w_im', W.imag)
    #     save_result(dir, 'theta', THETA)

    #kxmap
    # input_params = {
    #     'mts': '0', #multipole_type_selected
    #     'mos': '0', #multipole_order_selected
    #     'mps': '0', #m_projection_selected
    #     'type': '1',
    #     'order': '1',
    #     'm': '1',
    #     'rmax': 2,
    #     'lmax': 8,
    #     'lattice_constant': 300,
    #     'fab': 60,
    #     'polar': 'S',
    #     'epssph_re': 15,
    #     'epssph_im': 0,
    #     'ktype': 1,
    #     'kscan': 1,
    #     'zinf': 2.985,
    #     'zsup': 2.9851,
    #     'npts': 2,
    #     'r_ratio': 0.374
    # }
    # span = 0.3
    # n = 50
    # # kx = np.linspace(0.1, 0.2, n)
    # theta = np.linspace(50, 51, n)
    # # omega = np.linspace(2.983453, 2.9834532, 1000)/2/np.pi #0.121
    # # omega = np.linspace(2.9827, 2.9835, 2000)/2/np.pi #2.9829
    # omega = np.linspace(2.95, 3.3, 1000)/2/np.pi
    #
    # # T, R = calc_ak1_omega_map(ktype=1, ak1=theta, omega=omega, ip=input_params, is_w0_needed=False)
    # # show_2D_map(theta, omega, R, 1e-5, 1, logscale=True)
    # for ak1 in [50.4]:
    #     for ak2 in [0, 30, 60, 90]:
    #         dir = f'R_0.374_theta={ak1}_phi_{ak2}'
    #         F, T, R, A = calc_spectrum_omega(omega, ktype=1, ak1=ak1, ak2=ak2, ip=input_params)
    #         plt.plot(omega*2*np.pi, R, label='R')
    #         plt.plot(omega*2*np.pi, T, label='T')
    #         plt.plot(omega*2*np.pi, A, label='A')
    #         plt.yscale('log')
    #         plt.legend()
    #         plt.show()
            # save_result(dir, 'T', T)
            # save_result(dir, 'R', R)
            # save_result(dir, 'A', A)
            # save_result(dir, 'omega', omega*2*np.pi)


    # RMAX issue with triangular lattice
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
        'ktype': 1,
        'kscan': 1,
        'zinf': 0.5,
        'zsup': 1.0,
        'npts': 300,
        'r_ratio': 0.4705
    }

    omega = np.linspace(3.7275, 3.74, 300)/2/np.pi
    kx = np.linspace(0, 0.5, 100)

    for input_params['lmax'] in [4, 13]:
        time0 = time.time()
        for input_params['rmax'] in [7, 29]:
            # T, R = calc_ak1_omega_map(ktype=2, ak1=kx, omega=omega, ip=input_params, is_w0_needed=False)
            # show_2D_map(kx, omega, T, 0, 1)
            F, T, R, A = calc_spectrum_omega(omega, 2, 0.01, 0, input_params)
            # plt.plot(omega*2*np.pi, T, label=f'T_rmax{input_params["rmax"]}')
            # plt.legend()
            dir = f'rmax_issue/2layers/ak1_0.01/lmax={input_params["lmax"]}_rmax={input_params["rmax"]}'
            save_result(dir, 'T', T)
            # save_result(dir, 'R', R)
            save_result(dir, 'omega', omega)
            # save_result(dir, 'kx', kx*2)
        # plt.show()
        print(f'{time.time()-time0}s')