import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import numpy as np
# from colorsys import hls_to_rgb


def multi_loadtxt(dir, filelist):
    output = ()
    for fname in filelist:
        out = np.loadtxt(dir+"/"+fname)
        output += (out,)
    return output


def plot_2D_map(x, y, z, x_label='', y_label='', z_label='', 
                vmin=0, vmax=1, ax=None, logscale=False, 
                cmap=cm.viridis, is_cb_needed=False, extent=None):
    ax = ax or plt.gca()
    if logscale == True:
        norm = LogNorm(vmin=vmin, vmax=vmax)
        im = ax.imshow(z, interpolation='none', norm=norm, aspect='auto', cmap=cmap, extent=extent, origin='lower')
    else:
        im = ax.imshow(z, interpolation='none', vmin=vmin, vmax=vmax, aspect='auto', cmap=cmap, extent=extent, origin='lower')

    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    if is_cb_needed:
        cb = plt.colorbar(im)
        cb.set_label(z_label)
    return im




if __name__ == '__main__':
    # dir = '/home/ashalev/Projects/Multem-3/examples'
   
    # fig = plt.figure(figsize=(10,10))
    # x = np.linspace(500, 800, 10)
    # mz_smu_Tmatrix, mz_smu_original = multi_loadtxt(dir, ('M0_smu_Tmatrix.txt', 'M0_smu_original.txt'))
    # radius = 170.0 / 2.0
    # ax = fig.add_subplot(2, 2, 1)
    # ax.plot(x, mz_smu_original/3.14/radius**2, lw=4, color='lightblue', label='mz(origin_smuthi)')
    # ax.plot(x, mz_smu_Tmatrix/3.14/radius**2, ls='--', lw=2, color='black', label='mz(Tmatrix_mod)')
    # # ax.set_yscale('log')
    # ax.legend()
    # ax.set_title('GOLD SUBSTRATE TE (example from tutorial)')
    # ax.set_xlabel('wl [nm]')
    # ax.set_ylabel('ext cross section')

    # ax = fig.add_subplot(2, 2, 2)
    # ax.plot(x, mz_smu_original/3.14/radius**2, lw=4, color='lightblue')
    # ax.plot(x, mz_smu_Tmatrix/3.14/radius**2, ls='--', lw=2, color='black')
    # ax.set_yscale('log')
    # ax.set_xlabel('wl [nm]')
    # ax.set_ylabel('ext cross section')


    # pz_smu_Tmatrix, pz_smu_original = multi_loadtxt(dir, ('pz_smu_Tmatrix.txt', 'pz_smu_original.txt'))
    # radius = 170.0 / 2.0
    # ax = fig.add_subplot(2, 2, 3)
    # ax.plot(x, pz_smu_original/3.14/radius**2, lw=4, color='red', label='pz(origin_smuthi)')
    # ax.plot(x, pz_smu_Tmatrix/3.14/radius**2, ls='--', lw=2, color='black', label='pz(Tmatrix_mod)')
    # # ax.set_yscale('log')
    # ax.legend()
    # ax.set_title('GOLD SUBSTRATE TM (example from tutorial)')
    # ax.set_xlabel('wl [nm]')
    # ax.set_ylabel('ext cross section')

    # ax = fig.add_subplot(2, 2, 4)
    # ax.plot(x, pz_smu_original/3.14/radius**2, lw=4, color='red')
    # ax.plot(x, pz_smu_Tmatrix/3.14/radius**2, ls='--', lw=2, color='black')
    # ax.set_yscale('log')
    # ax.set_xlabel('wl [nm]')
    # ax.set_ylabel('ext cross section')



    dir = '/home/ashalev/Projects/Multem-3/examples/smuthi_multem_comparison/all_multipoles'
    T_all = np.loadtxt(f'{dir}/T_all.txt')[0]
    T_t1_l1 = np.loadtxt(f'{dir}/T_t1_l1.txt')[0]
    T_t1_l2 = np.loadtxt(f'{dir}/T_t1_l2.txt')[0]
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1, 1, 1)
    x = np.arange(0, 90, 10)
    ax.plot(x, T_all, lw=4, color='red')
    ax.plot(x, T_t1_l2, lw=2, ls='--', color='red')



    plt.show()



