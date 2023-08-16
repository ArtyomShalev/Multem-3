import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
import numpy as np
import yaml


def read_refractive_index_from_yaml(filename, vacuum_wavelength, units="mkm", kind=1):
    """Read optical constants in format provided by refractiveindex.info website.

    Args:
            filename (str): path and file name for yaml data
                            downloaded from refractiveindex.info
            vacuum_wavelength (float or np.array): wavelengths where refractive
                                                   index data is needed
            units (str): units for wavelength. currently, microns ('mkm' or 'um')
                         and nanometers ('nm') can be selected
            kind (int): order of interpolation

    Returns:
            A pair (or np.array of pairs) of wavelength and
            corresponding refractive index (complex)
    """
    if units == "nm":
        factor = 1000
    elif units in ("mkm", "um"):
        factor = 1
    else:
        raise NotImplementedError("Converting wavelength into '"+units
                                  +"' units for refractive index data"
                                  +" was not implemented.")

    with open(filename) as f:
        the_file = yaml.load(f, yaml.SafeLoader)['DATA'][0]
    data_type = the_file['type']
    if data_type != 'tabulated nk':
        raise NotImplementedError("Input data type '"+data_type
                                  +"' available in file "+filename
                                  +" was not implemented.")
    data = the_file['data'].splitlines()
    data_split = []
    for wl in data:
        data_split.append(wl.split())
    data_num = []
    for wl in data_split:
        record = []
        for val in wl:
            record.append(float(val))
        data_num.append(record)
    data_np = np.array(data_num)
    data_wl = data_np[:,0]*factor
    index_re = data_np[:,1]
    index_im = data_np[:,2]
    f_re = interp1d(data_wl, index_re, kind=kind)
    f_im = interp1d(data_wl, index_im, kind=kind)
    data_out = np.transpose(np.vstack((vacuum_wavelength, f_re(vacuum_wavelength)+f_im(vacuum_wavelength)*1j)))
    if len(data_out) == 1:
        return data_out[0]
    return data_out


def create_input(ktype, kscan, npts, fab, ak1, ak2, zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im, m_type, m_order, m,
                 is_multipole_type_selected, is_multipole_order_selected, is_m_projection_selected, mode):
    '''
    mode:   1 - 2D array of spheres in homogeneous media
            2 - interface
    '''

    if mode == '1':
        str_fort10 = ('           ********************************************\n'
                      '           ********INPUT FILE FOR TRANSMISSION*********\n'
                      '           ********************************************\n'
                      '   KTYPE ='+'%2i'%(ktype)+'   KSCAN ='+'%2i'%(kscan)+'   KEMB  = 0    LMAX ='+'%2i'%(lmax)+'   NCOMP = 1   NUNIT = 1\n'
                      ' ALPHA =    1.000000  BETA =    1.000000   FAB =   '+'%9.6f'%(fab)+'  RMAX ='+'%11.6f'%(rmax)+'\n'
                      '  NP ='+'%4i'%(npts)+'  ZINF ='+
                      '%19.15f'%(zinf)+'  ZSUP ='+'%19.15f'%(zsup)+'\n'
                      '  THETA/AK(1) ='+'%19.15f'%(ak1)+'     FI/AK(2) ='+'%19.15f'%(ak2)+'   POLAR ='+polar+'     FEIN =   0.00\n'
                      '\n'
                      'Give information for the "NCOMP" components \n'
                      '\n'
                      '     IT  = 2\n'
                      '     MUMED =   1.00000000   0.00000000     EPSMED=   2.2500000   0.00000000\n'
                      '   NPLAN = 1  NLAYER = 1\n'
                      '       S =   '+'%10.8f'%(r_ratio)+'     MUSPH =   1.00000000   0.00000000     EPSSPH=  '+'%11.6f'%(epssph_re)+'   '+'%11.8f'%(epssph_im)+'\n'
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



def eval(i=0, mode='', multem_version='3'):
    #mode: WIF - With Input Function;
    if mode == 'WIF':
        create_input(ktype, kscan, npts, fab, ak1[i], ak2[i], zinf, zsup, polar, lmax, r_ratio, rmax, epssph_re, epssph_im, m_type, m_order, m,
                     is_multipole_type_selected, is_multipole_order_selected, is_m_projection_selected, mode='1')
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
            # print("Running multem...")
            subprocess.run(['./multem2'],
                           stdout=subprocess.DEVNULL,
                           env=my_env)
        # print("Done.")
    #FREQUENCY   TRANSMITTANCE  Reflectance   Absorbance
    d = np.loadtxt('fort.8').T
    return d


def copy_input(src, target):
    #erase target file
    if os.stat(target).st_size != 0:
        open(target, 'w').close()

    with open(src,'r') as firstfile, open(target,'a') as secondfile:
        # read content from first file
        for line in firstfile:
            # append content to second file
            secondfile.write(line)


def save_1D_data(x, y, filename, format='%e'):
    data = np.stack((x, y), axis=-1)
    np.savetxt(filename, data, delimiter='\t', fmt=format)
    print('data were saved')


if __name__ == '__main__':
    # --------- Colloquium: Light scattering by particle and hole arrays by F. J. Garc ́ıa de Abajo ------------
    npts = 300
    lam_start = 715 #715 #nm
    lam_end = 770 #770 #nm
    x = np.linspace(lam_start, lam_end, npts)
    eps = []
    for wl in x:
        nk = read_refractive_index_from_yaml('/home/ashalev/Projects/amos-try/nk-data/Ag-Johnson.yml', wl, units='nm')[1]
        eps.append(nk**2)
    eps = np.array(eps)
    epssph_re = eps.real
    epssph_im = eps.imag
    a = 500 #nm
    A = []
    for i in range(epssph_re.shape[0]):
        print(f'Progress {i} from {epssph_re.shape[0]}')
        create_input(2, 2, 2, 90.0, 0, 0, x[i]/a, x[i]/a+0.001*x[i]/a, 'S', 7, 0.2, 20.0, epssph_re[i], epssph_im[i],
                     '0', '0', '0', '0', '0', '0', mode='1')
        A.append(eval(multem_version='3')[3][0])
    plt.plot(x, A, 'b')
    plt.yscale('log')
    plt.show()



    # fig19(i)
    src = r'/home/ashalev/Projects/amos-try/multem2mod/multem3article/input_files/multem3/deAbajo/fig19(i)'
    target = r'/home/ashalev/Projects/amos-try/multem2mod/multem3article/fort.10'
    copy_input(src, target)
    wl_i, T_i, R_i, A_i = eval()
    # fig19(ii)
    src = r'/home/ashalev/Projects/amos-try/multem2mod/multem3article/input_files/multem3/deAbajo/fig19(ii)'
    target = r'/home/ashalev/Projects/amos-try/multem2mod/multem3article/fort.10'
    copy_input(src, target)
    wl_ii, T_ii, R_ii, A_ii = eval()

    npts = 300
    lam_start = 715 #715 #nm
    lam_end = 770 #770 #nm
    x = np.linspace(lam_start, lam_end, npts)
    plt.plot(x, A_i, 'b')
    plt.plot(x, A_ii, 'g--')
    plt.yscale('log')
    plt.ylim([0.002, 0.14])
    plt.show()

    # --------- Non-local optical response of two-dimensional arrays of metallic nanoparticles by Vassilios Yannopapas ------------

    # don't know yet how to calculate



    # --------- Bound States in the Continuum and Fano Resonances in the Dirac Cone Spectrum by Bulgakov, Maksimov ------------
    # ---- multem3 -----
    plt.figure('Maksimov_Bulgakov_fig3',figsize=(11,13*6))
    ktype = 2
    kscan = 1
    fab = 60
    rmax = 5
    polar = 'S'
    epssph_re = 15
    epssph_im = 0
    is_multipole_type_selected = '0'
    is_multipole_order_selected = '0'
    is_m_projection_selected = '0'
    m_type = '0 0'
    m_order = '1 1'
    m = '-1 1'
    npts = 300
    r_ratio = 0.47050000
    from_k = 3.73
    to_k = 3.7405
    zinf = from_k
    zsup = to_k
    ak1 = np.array((0.01, 0.05, 0.1, 0))/2/np.pi
    ak2 = np.array((0, 0, 0, 0.05))/2/np.pi
    kpts = len(ak1)
    for lmax in range(1, 13, 1):
        T = []
        for i in range(kpts):
            print(i+1, 'of', kpts)
            T.append(eval(i, mode='WIF', multem_version='3')[1])
        x = np.linspace(from_k, to_k, npts)
        step = -1.7
        # step = 0
        plt.annotate(f'lmax={lmax}', (np.min(x), np.min(T[0])+step*((lmax-1)+0.2)), fontsize=12)
        plt.plot(x, T[0]+step*(lmax-1), 'r')
        plt.plot(x, T[3]+step*(lmax-1), color='silver', lw=2)
        plt.plot(x, T[1]+step*(lmax-1), 'b--', dashes=(15, 15))
        plt.plot(x, T[2]+step*(lmax-1), color='magenta', linestyle= (0, (3, 5, 1, 5)))
        plt.ylabel('T')
        plt.xlabel('ak0')
    plt.savefig(f'multem3_lmax{lmax}_MB_fig3_rmax{rmax}.pdf')
    plt.clf(); plt.close()


    # ----- narrow red spectrum in details -------
    plt.figure('Maksimov_Bulgakov_fig3',figsize=(11,6*10))
    from_k = 3.730
    to_k = 3.7304
    zinf = from_k
    zsup = to_k
    ak1 = [0.01]
    ak2 = [0.0]
    x = np.linspace(from_k, to_k, npts)
    lmax_start = 7
    lmax_end = 18
    for lmax in range(lmax_start, lmax_end, 1):
        T = (eval(mode='WIF', multem_version='3')[1])
        step = -1.7
        save_1D_data(x, T, f'spectra/Maksimov_Bulgakov/lmaxx={lmax}.txt', format='%19.16e')
        plt.annotate(f'lmax={lmax}', (np.min(x), (lmax-lmax_start)*step), fontsize=12)
        plt.plot(x, T+step*(lmax-lmax_start), 'r')
        plt.ylabel('T')
        plt.xlabel('ak0')
    plt.savefig(f'multem3_spectrum_lmax{lmax}_MB_fig3_rmax{rmax}.pdf')
    plt.clf(); plt.close()

    # narrow spectra ak0 vs lmax visualisation
    plt.figure(figsize=(11,6))
    x = np.loadtxt('spectra/Maksimov_Bulgakov_output/RLs-fano.txt')
    y = np.loadtxt('spectra/Maksimov_Bulgakov_output/f0.txt')
    plt.scatter(x, y)
    plt.ylabel('resonance position ak0')
    plt.xlabel('lmax')
    plt.savefig(f'ak0_vs_lmax.pdf')
    plt.clf(); plt.close()



    #---- multem2 ------
    # rmax = 20
    # plt.figure('Maksimov_Bulgakov_fig3',figsize=(11,6))
    # target = r'/home/ashalev/Projects/amos-try/multem2mod/multem3article/fort.10'
    # T = []
    # # ---- kx = 0.01 ky = 0. -----
    # src = r'/home/ashalev/Projects/amos-try/multem2mod/multem3article/input_files/multem2/Maksimov_Bulgakov_fig3/kx=0,01.txt'
    # copy_input(src, target)
    # T.append(eval(multem_version='2')[1])
    # # ---- kx = 0.05 ky = 0. -----
    # src = r'/home/ashalev/Projects/amos-try/multem2mod/multem3article/input_files/multem2/Maksimov_Bulgakov_fig3/kx=0,05.txt'
    # copy_input(src, target)
    # T.append(eval(multem_version='2')[1])
    # # ---- kx = 0.1 ky = 0. -----
    # src = r'/home/ashalev/Projects/amos-try/multem2mod/multem3article/input_files/multem2/Maksimov_Bulgakov_fig3/kx=0,1.txt'
    # copy_input(src, target)
    # T.append(eval(multem_version='2')[1])
    # # ---- kx = 0. ky = 0.05 -----
    # src = r'/home/ashalev/Projects/amos-try/multem2mod/multem3article/input_files/multem2/Maksimov_Bulgakov_fig3/ky=0,05.txt'
    # copy_input(src, target)
    # T.append(eval(multem_version='2')[1])
    #
    # from_k = 3.73
    # to_k = 3.7405
    # npts = 300
    # lmax = 4
    # x = np.linspace(from_k, to_k, npts)
    # # step = -1.7
    # step = 0
    # plt.annotate(f'lmax={lmax}', (np.min(x), np.min(T[0])+step*((lmax-1)+0.2)), fontsize=12)
    # plt.plot(x, T[0]+step*(lmax-1), 'r')
    # plt.plot(x, T[3]+step*(lmax-1), color='silver', lw=2)
    # plt.plot(x, T[1]+step*(lmax-1), 'b--', dashes=(15, 15))
    # plt.plot(x, T[2]+step*(lmax-1), color='magenta', linestyle= (0, (3, 5, 1, 5)))
    # plt.ylabel('T')
    # plt.xlabel('ak0')
    # plt.savefig(f'multem2_lmax{lmax}_MB_fig3_rmax{rmax}.pdf')
    # plt.clf(); plt.close()


    # --------- Multipolar Lattice Resonances in Plasmonic Finite-Size Metasurfaces by Kostyukov ------------
    # ---- multem3 -----
    # plt.figure('Maksimov_Bulgakov_fig3',figsize=(11,6))
    # ktype = 2
    # kscan = 2
    # rmax = 20
    # fab = 90
    # polar = 'S'
    # epssph_re = -7.0124
    # epssph_im = 0.21187
    # is_multipole_type_selected = '0'
    # is_multipole_order_selected = '0'
    # is_m_projection_selected = '0'
    # m_type = '0 0'
    # m_order = '1 1'
    # m = '-1 1'
    # npts = 100
    # r = 60
    # a = 280
    # r_ratio = r/a
    # from_wl = 400
    # to_wl = 500
    # zinf = from_wl/a
    # zsup = to_wl/a
    # lmax = 7
    # ak1 = [0]
    # ak2 = [0]
    # wl, T, R, A = eval(mode='WIF', multem_version='3')
    # x = np.linspace(from_wl, to_wl, npts)
    # plt.plot(x, -np.log(T))
    # # plt.plot(x, T)
    # # plt.plot(x, A)
    # plt.show()



