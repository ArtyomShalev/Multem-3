import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.patches import Rectangle
from matplotlib.offsetbox import AnchoredText
from matplotlib.patches import Arc
import seaborn as sns
from matplotlib.gridspec import GridSpec
from matplotlib.patches import ConnectionPatch
from matplotlib.transforms import Bbox, TransformedBbox, \
    blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector,\
    BboxConnectorPatch
from matplotlib.colors import LogNorm


def multi_loadtxt(dir, filelist):
    output = ()
    for fname in filelist:
        out = np.loadtxt(dir+"/"+fname)
        output += (out,)
    return output


def connect_bbox(bbox1, bbox2,
                 loc1a, loc2a, loc1b, loc2b,
                 prop_lines, prop_patches=None):
    if prop_patches is None:
        prop_patches = prop_lines.copy()
        prop_patches["alpha"] = prop_patches.get("alpha", 1)*0.2

    c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
    c1.set_clip_on(False)
    c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
    c2.set_clip_on(False)

    bbox_patch1 = BboxPatch(bbox1, **prop_patches)
    bbox_patch2 = BboxPatch(bbox2, **prop_patches)

    p = BboxConnectorPatch(bbox1, bbox2,
                           # loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                           loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                           **prop_patches)
    p.set_clip_on(False)

    return c1, c2, bbox_patch1, bbox_patch2, p


def zoom_effect01(ax1, ax2, xmin, xmax, **kwargs):
    """
    ax1 : the main axes
    ax1 : the zoomed axes
    (xmin,xmax) : the limits of the colored area in both plot axes.

    connect ax1 & ax2. The x-range of (xmin, xmax) in both axes will
    be marked.  The keywords parameters will be used ti create
    patches.

    """

    trans1 = blended_transform_factory(ax1.transData, ax1.transAxes)
    trans2 = blended_transform_factory(ax2.transData, ax2.transAxes)

    bbox = Bbox.from_extents(xmin, 0, xmax, 1)

    mybbox1 = TransformedBbox(bbox, trans1)
    mybbox2 = TransformedBbox(bbox, trans2)

    prop_patches = kwargs.copy()
    prop_patches["ec"] = "none"
    prop_patches["alpha"] = 0.1
    kwargs["alpha"] = 0.3
    kwargs["ls"] = '--'
    c1, c2, bbox_patch1, bbox_patch2, p = \
        connect_bbox(mybbox1, mybbox2,
                     loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                     prop_lines=kwargs, prop_patches=prop_patches)

    ax1.add_patch(bbox_patch1)
    ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p


def prep_xy_data(data):
    x, y = [], []
    for line in data:
        xi, yi = line[0], line[1]
        x.append(xi)
        y.append(yi)
    return x, y


def load_errors(filename, err_lim):
    #filnames - list of str txt files names
    content = np.loadtxt(filename)
    content[content<err_lim] = err_lim
    return content


def convert_to_square_2d(arr):
    square_shape = int(np.sqrt(len(arr))), int(np.sqrt(len(arr)))
    arr = arr.reshape(square_shape)
    return arr


def plot_error(ax, x,y,z, title='', vextr=None, n_orders_div=1):
    if vextr is not None:
        z = np.clip(a=z, a_min=vextr[0], a_max=vextr[1])
        n_orders = int(np.ceil(np.abs(np.log10(vextr[0]) -
                                      np.log10(vextr[1])))) // n_orders_div #// 2 + 1
        cmap_idv = matplotlib.colors.ListedColormap(sns.color_palette("RdBu_r", n_orders))
        im  = ax.imshow(
            z.astype(np.float64),
            origin="lower",extent=[exp_min,exp_max,exp_min,exp_max],aspect='auto',
            norm=LogNorm(vmin=vextr[0], vmax=vextr[1]),
            cmap=cmap_idv)
    else:
        im = ax.imshow(np.vectorize(mpmath.log10)(z).astype(np.float64),
                       origin="bottom",extent=[exp_min,exp_max,exp_min,exp_max],aspect='auto',
                       cmap=cmap_std)

    ax.set_title(title)
    return im


def create_arc_and_arrow(ax, x1, y1, len_x, len_y, fc, alpha):
    ax.arrow(x1,
             y1,
             len_x,
             len_y,
             head_width=0.03,
             head_length=0.0001,
             color=fc,
             alpha=alpha,
             linewidth=2,
             ec=fc,
             fc=fc)
 

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


#--------- customizing font --------
plt.rcParams.update({'font.size': 28, 'font.serif':"Times New Roman"})


figures_to_plot = ['fig3']

# ------------- Fig. 4 ---------------------------------------------------------------------------------------------
if 'fig4' in figures_to_plot:
    fig = plt.figure(figsize=(20, 20))
    gs = GridSpec(2, 2, figure=fig)
    # # fig 3 a
    ax1 = fig.add_subplot(gs[0, :])
    core_dir = 'data/fig3'
    # red
    data = np.loadtxt(core_dir+'/akxy=0.01_0.0.txt')
    x, y = prep_xy_data(data)
    ax1.plot(x, y, color='red', label=r'$0.01$'+'\t'+ r'$0$', lw=6)
    #silver
    data = np.loadtxt(core_dir+'/akxy=0.0_0.05.txt')
    x, y = prep_xy_data(data)
    ax1.plot(x, y, color='silver', label=r'$0$'+'\t'+ '  0.05', lw=6)
    #blue dashed
    data = np.loadtxt(core_dir+'/akxy=0.05_0.0.txt')
    x, y = prep_xy_data(data)
    ax1.plot(x, y, 'b--', dashes=(15, 15), label=r'$0.05$'+'\t'+ r'$0$', lw=6)
    #magenta dot-dashed
    data = np.loadtxt(core_dir+'/akxy=0.1_0.0.txt')
    x, y = prep_xy_data(data)
    ax1.plot(x, y, color='magenta', linestyle= (0, (3, 5, 1, 5)), label=r'$0.1$'+'\t'+' 0', lw=6)
    #------- customizing
    ax1.set_yticks([0.0, 1.0], minor=False)
    ax1.set_xlabel(r'$ak_0$')
    ax1.set_ylabel(r'$\mathbf{T}$', rotation=0, fontsize=30)
    ax1.set_xlim(3.7285, 3.74)
    ax1.set_ylim(-0.1, 1.1)
    ax1.legend(frameon=False, loc=2, bbox_to_anchor=(0.0,0.86) ,title='\t'+r'$ak_x$'+'\t'+r'$ak_y$')
    # fig.3 b
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1], sharey=ax2)  
    for lmax in [4, 7, 10, 13]:
        data = np.loadtxt(core_dir+f'/lmax={lmax}.txt')
        x, y = prep_xy_data(data)
        if lmax == 4:
            ax2.plot(x, y, label=str(lmax), color='red', lw=6)
            ax3.plot(x, y, label=str(lmax), color='red', lw=6)
        else:
            ax2.plot(x, y, label=str(lmax), lw=6)
            ax3.plot(x, y, label=str(lmax), lw=6)
    ax2.set_ylim(-0.1, 1.1)
    ax2.set_xlim(3.7317, 3.7325)
    ax3.set_xlim(3.7348, 3.7356)
    ax2.legend(frameon=False, loc=2, title='\t'+r'$l_{max}$')
    # hiding the spines between ax2 and ax3
    ax2.spines['right'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.yaxis.tick_right()
    ax3.tick_params(labelright=False)  # don't put tick labels at the top
    # drawing diagonal lines
    d = .015  # how big to make the diagonal lines in axes coordinates
    kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
    ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    kwargs.update(transform=ax3.transAxes)
    ax3.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax3.plot((-d, +d), (-d, +d), **kwargs)
    #------- customizing
    ax2.set_ylabel(r'$\mathbf{T}$', rotation=0, fontsize=30)
    ax2.set_yticks([0.0, 1.0], minor=False)
    ax2.set_xticks([3.7318, 3.7321, 3.7324], minor=False)
    ax3.set_xticks([3.735, 3.7353, 3.7356], minor=False)
    zoom_effect01(ax1, ax3, 3.7349, 3.7355)
    
    # customizing ticks 
    for ax in [ax1, ax2, ax3]:
        ax.tick_params('both', length=10, width=2, which='major')
    
    fig.text(0.05, 0.9, 'a', ha='center')
    fig.text(0.05, 0.47, 'b', ha='center')
    fig.text(0.51, 0.06, r'$ak_0$', ha='center')



    plt.savefig('fig3.pdf')
    plt.clf(); plt.close()


if 'fig3' in figures_to_plot:
    fig = plt.figure(figsize=(16*1.25, 9*1.25))
    ax1 = fig.add_subplot(1, 1, 1)
    data_dir = 'data/fig3/'
    data = np.loadtxt(data_dir+'T.txt')
    x, T = prep_xy_data(data)
    x, T = x, T
    error_wo_lapack = prep_xy_data(np.loadtxt(data_dir+'error_wo_lapack.txt'))[1]
    # error_wo_lapack_new = prep_xy_data(np.loadtxt(data_dir+'error_wo_lapack_my.txt'))[1]

    error_with_lapack = prep_xy_data(np.loadtxt(data_dir+'error_with_lapack.txt'))[1]
    # error_with_lapack_and_fad = prep_xy_data(np.loadtxt(data_dir+'error_with_lapack.txt'))[1]
    ax11 = ax1.twinx()
    l1, = ax11.plot(x[::2], np.array(error_with_lapack[::2])*1e4, color='C2', lw=3)
    l2, = ax11.plot(x, np.array(error_wo_lapack)*1e1, color='C0', lw=6)
    ax1.set_ylabel('error')
    l3, = ax1.plot(x, T, 'red', lw=18, alpha=0.3)
    # create_arc_and_arrow(ax11, 3.733282, 1.002)
    ax1.set_ylabel(r'$\mathbf{T}$', rotation=0, fontsize=30)
    plt.legend([l1, l2], ['LAPACK', 'W/O LAPACK'], frameon=False)
    # plt.annotate(r"x$10^4$",
    #             xy=(3.731546, 0.6), xycoords='data',    # 0.6 1 option; 0.34 2 option
    #             xytext=(-19, 30), textcoords='offset points',
    #             arrowprops=dict(arrowstyle="->", color='black')
                # )
    # customizing
    ax1.set_yticks([0.0, 1.0], minor=False)
    ax1.set_ylim([-0.25, 1.1])
    ax1.set_xlim([3.731, 3.735])
    ax1.set_xlabel(r'$ak_0$')
    ax1.set_xticks([3.731, 3.732, 3.733, 3.734, 3.735], minor=False)
    ax11.set_yticks([0, 0.04], minor=False)  # -0.0235 1 option; -0.005 2 option
    ax11.set_yticklabels(['0', '0.04'])
    ax11.set_ylabel('error', labelpad=-30)
    arrow_len = 0.00012
    create_arc_and_arrow(ax1, 3.7321, 0.9, -arrow_len, 0.0, 'pink', 1.0)
    create_arc_and_arrow(ax1, 3.73465, 0.17, arrow_len, 0.0, 'C2', 1.0)
    create_arc_and_arrow(ax1, 3.73249, -0.18, arrow_len, 0.0, 'C0', 1.0)
    # annotations on the whole figure
    # fig.text(0.49, 0.02, r'$ak_0$', ha='center')
    # fig.text(0.17, 0.28, r'$l_{max}=10$', ha='center')
    # fig.text(0.17, 0.93, r'$l_{max}=4$', ha='center')
    # fig.text(0.02, 0.97, 'a', ha='center')
    # fig.text(0.02, 0.48, 'b', ha='center')
    fig.text(0.81, 0.35, r'$x10^4$', ha='center', color='C2', fontsize=30)
    fig.text(0.36, 0.13, r'$x10$', ha='center', color='C0', fontsize=30)

    # ---------- zooming
    # zoom_effect01(ax2, ax4, 3.7321, 3.7322)
    # plt.subplot_tool()
    plt.subplots_adjust(left=0.07,
                        bottom=0.1,
                        right=0.915,
                        top=0.963,
                        wspace=0.033,
                        hspace=0.297)

      # customizing ticks 
    for ax in [ax1, ax11]:
        ax.tick_params('both', length=10, width=2, which='major')
    plt.savefig('fig3.pdf')
    plt.clf(); plt.close()


if 'fig6_7' in figures_to_plot:
    plt.rcParams.update({'font.size': 50, 'font.serif':"Times New Roman"})
    fig = plt.figure(figsize=(40, 20))
    # --- reading data ------
    data_dir = 'data/fig6/'
    x, y = np.loadtxt(data_dir+'x.txt'), np.loadtxt(data_dir+'y.txt')
    x2d = convert_to_square_2d(x)
    y2d = convert_to_square_2d(y)
    err_lim = 1e-16
    multem3_re_err = convert_to_square_2d(load_errors(data_dir+'MIT_re_err.txt', err_lim))
    multem3_im_err = convert_to_square_2d(load_errors(data_dir+'MIT_im_err.txt', err_lim))
    multem2_re_err = convert_to_square_2d(load_errors(data_dir+'multem2_re_err.txt', err_lim))
    multem2_im_err = convert_to_square_2d(load_errors(data_dir+'multem2_im_err.txt', err_lim))
    cmap_std = matplotlib.colors.ListedColormap(sns.color_palette("RdBu_r", 16))
    exp_min = -8
    exp_max = 8
    implementations = [multem2_re_err, multem2_im_err, multem3_re_err, multem3_im_err]
    impls_str = ['multem2 re err', 'multem2 im err', 'multem3 re err', 'multem3 im err']

    # --- reading args for Faddeeva func from MB Dirac cone
    with open(data_dir + 'fort.1') as f:
        lines = f.readlines()
    args_re, args_im = [], []
    # ---- postprocess - deleting brackets and left only unique positive values
    for line in lines:
        # line = line.replace("(","")
        # line = line.replace(")","")
        re, im = line.split()
        args_re.append(re)
        args_im.append(im)
    args_re = np.array(args_re, dtype=np.float64)
    args_im = np.array(args_im, dtype=np.float64)
    z = args_re + 1j*args_im
    Z = np.unique(z)
    pos_args = []
    for z in Z:
        if z.real > 0 and z.imag > 0:
            pos_args.append(z)
    pos_args = np.array(pos_args)

    # ax1
    ax1 = fig.add_subplot(2,2,1)
    im = plot_error(ax1, x2d, y2d, multem2_re_err, vextr=[err_lim, 1e-6])
    ax1.get_xaxis().set_visible(False)
    ax1.scatter(np.log10(pos_args.real)[::60], np.log10(pos_args.imag)[::60], color='black', s=15)
    ax1.set_yticks([-8, 0, 8])
    ax1.set_yticklabels([])
    # plotting argument regions
    regimes = np.loadtxt(data_dir+'multem2_regimes.txt')
    regimes2d = convert_to_square_2d(regimes)
    ax1.contour(regimes2d, levels=2, linewidths=7, colors='yellow', linestyles='dashed', extent=[exp_min,exp_max,exp_min,exp_max])

    # ax2
    ax2 = fig.add_subplot(2,2,2, sharey=ax1)
    im = plot_error(ax2, x2d, y2d, multem2_im_err, vextr=[err_lim, 1e-6])
    ax2.get_yaxis().set_visible(False)
    ax2.get_xaxis().set_visible(False)
    # plotting argument regions
    ax2.contour(regimes2d, levels=2, linewidths=7, colors='yellow', linestyles='dashed', extent=[exp_min,exp_max,exp_min,exp_max])


    # ax3
    ax3 = fig.add_subplot(2,2,3)
    im = plot_error(ax3, x2d, y2d, multem3_re_err, vextr=[err_lim, 1e-6])
    ax3.set_yticks([-8, 0, 8])
    ax3.set_yticklabels([])
    ax3.set_xticks([-8, 0, 8])
    ax3.set_xticklabels([])

    # ax4
    ax4 = fig.add_subplot(2,2,4, sharey=ax3)
    im = plot_error(ax4, x2d, y2d, multem3_im_err, vextr=[err_lim, 1e-6])
    ax4.get_yaxis().set_visible(False)
    ax4.set_xticks([-8, 0, 8])
    ax4.set_xticklabels([])


    # ax5 (arg values from MB Dirac cone)
    ax5 = fig.add_axes([0.34, 0.6, 0.15, 0.35])
    # plotting argument regions
    ax5.contour(regimes2d, zorder=1, levels=2, linewidths=7, colors='yellow', linestyles='dashed', extent=[exp_min,exp_max,exp_min,exp_max])
    im = plot_error(ax5, x2d, y2d, multem2_re_err, vextr=[err_lim, 1e-6])
    ax5.scatter(np.log10(pos_args.real)[::50], np.log10(pos_args.imag)[::50], color='black', s=30)
    ax5.set_xticks([-0.0092, -0.0074])
    ax5.set_xticklabels([round(np.min(pos_args.real),3), round(np.max(pos_args.real),3)])
    ax5.set_yticks([np.min(np.log10(pos_args.imag)), np.max(np.log10(pos_args.imag))])
    ax5.set_yticks([np.log10(2), np.log10(16)])    
   

    ax5.set_yticklabels([2, 16])
    ax5.set_xlim([np.min(np.log10(pos_args.real)) - (np.max(np.log10(pos_args.real)) - np.min(np.log10(pos_args.real)))/2, np.max(np.log10(pos_args.real)) + (np.max(np.log10(pos_args.real)) - np.min(np.log10(pos_args.real)))/2])
    ax5.set_ylim([np.min(np.log10(pos_args.imag)) - (np.max(np.log10(pos_args.imag)) - np.min(np.log10(pos_args.imag)))/2, np.max(np.log10(pos_args.imag)) + (np.max(np.log10(pos_args.imag)) - np.min(np.log10(pos_args.imag)))/2])

    for ax_side in ['bottom', 'top', 'left', 'right']:
        ax5.spines[ax_side].set_color('white')
        ax5.spines[ax_side].set_linewidth(2)


    ax1.annotate('', xy=(0.48, 0.55),  xycoords='axes fraction',
                xytext=(0.8, 0.56), textcoords='axes fraction',
                arrowprops=dict(width=0.1, headwidth=3, edgecolor='white', facecolor='white', shrink=0.1, lw=2),
                horizontalalignment='right', verticalalignment='top',
                )

    ax1.add_patch(Rectangle((-0.09, 0.1), 0.16, 1.3, color='white', fill=None, alpha=1, lw=2))

    # colorbar
    cb_ax = fig.add_axes([0.92, 0.07, 0.01, 0.895])
    cbar = fig.colorbar(im, cax=cb_ax)
    cbar.set_label('relative error', labelpad=10, rotation=90)

    # annotation and custom axis ticks labeling
    fig.text(0.495, 0.005, r'$\Re(z)$', ha='center')

    left_x_pos = 0.01
    fig.text(left_x_pos, 0.5, r'$\Im(z)$', ha='center', rotation=90)
    fig.text(0.28, 0.97, r'$\Re(w(z))$', ha='center')
    fig.text(0.73, 0.97, r'$\Im(w(z))$', ha='center')

    delta_x_pos = 0.04
    fig.text(left_x_pos + delta_x_pos, 0.94, r'$10^{8}$', ha='center')
    fig.text(left_x_pos + delta_x_pos, 0.735, r'$10^{0}$', ha='center')
    fig.text(left_x_pos + delta_x_pos, 0.54, r'$10^{-8}$', ha='center')
    fig.text(left_x_pos + delta_x_pos, 0.49, r'$10^{8}$', ha='center')
    fig.text(left_x_pos + delta_x_pos, 0.3, r'$10^{0}$', ha='center')
    fig.text(left_x_pos + delta_x_pos, 0.09, r'$10^{-8}$', ha='center')

    bottom_y_pos = 0.045
    fig.text(0.08, bottom_y_pos, r'$10^{-8}$', ha='center')
    fig.text(0.29, bottom_y_pos, r'$10^{0}$', ha='center')
    fig.text(0.48, bottom_y_pos, r'$10^{8}$', ha='center')
    fig.text(0.515, bottom_y_pos, r'$10^{-8}$', ha='center')
    fig.text(0.705, bottom_y_pos, r'$10^{0}$', ha='center')
    fig.text(0.9, bottom_y_pos, r'$10^{8}$', ha='center')


    #annotating argument regions
    fontsize = 55

    # (i)
    fig.text(0.1, 0.67, r'$(i)$', ha='center', fontsize=fontsize)
    fig.text(0.52, 0.67, r'$(i)$', ha='center', fontsize=fontsize)
    fig.text(0.475, 0.64, r'$(i)$', ha='center', fontsize=fontsize)
    # (ii)
    ax1.annotate(r'$(ii)$', xy=(0.05, 0.54),  xycoords='axes fraction',
                xytext=(0.2, 0.7), textcoords='axes fraction', fontsize=fontsize,
                arrowprops=dict(width=0.1, headwidth=6, edgecolor='white', facecolor='white', shrink=0.1, lw=6),
                horizontalalignment='right', verticalalignment='top',
                )
    ax2.annotate(r'$(ii)$', xy=(0.05, 0.54),  xycoords='axes fraction',
                xytext=(0.2, 0.7), textcoords='axes fraction', fontsize=fontsize,
                arrowprops=dict(width=0.1, headwidth=6, edgecolor='white', facecolor='white', shrink=0.1, lw=6),
                horizontalalignment='right', verticalalignment='top',
                )
    fig.text(0.475, 0.75, r'$(ii)$', ha='center', fontsize=fontsize)

    # (iii)
    fig.text(0.1, 0.82, r'$(iii)$', ha='center', fontsize=fontsize)
    fig.text(0.52, 0.82, r'$(iii)$', ha='center', fontsize=fontsize)
    fig.text(0.475, 0.87, r'$(iii)$', ha='center', fontsize=fontsize)

    # figure letters
    for letter, ax in zip(['a', 'b', 'c', 'd', 'e'], [ax1, ax2, ax3, ax4, ax5]):
        at = AnchoredText(
            letter, prop=dict(size=fontsize, fontweight="bold"), frameon=False, loc='upper left')
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.3")
        ax.add_artist(at)
    plt.subplots_adjust(left=0.08,
                        bottom=0.099,
                        right=0.91,    
                        top=0.96,
                        wspace=0.01,
                        hspace=0.02)


    #make ticks bigger
    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.tick_params('both', length=15, width=4, which='major')

    plt.savefig('fig6.pdf')
    plt.clf(); plt.close()


import calculating_lib as cl

if 'fig5' in figures_to_plot:

    dir = f'data/fig5'
    fig = plt.figure(figsize=(20, 20))
    i = 1
    curve_styles = {
        '7': (8, 1, 'red', 1, 'dotted'),
        '12': (18, 0.4, 'green', 1, '-'),
        '14': (6, 1, 'black', 1, '-'),
    }
    RMAX = ['12', '14', '7'] 
    lines = {}

    for rmax in RMAX:
        R1, R2, R3, x = multi_loadtxt(dir, (f'650_{rmax}.txt', f'750_{rmax}.txt', f'900_{rmax}.txt', 'sintheta.txt'))
        plt.plot(x, R1, color=curve_styles[rmax][2], lw=curve_styles[rmax][0], alpha=curve_styles[rmax][1], ls=curve_styles[rmax][4])
        plt.plot([0.63, 0.63], [0, 0.15], 'red', lw=5)
        plt.plot(x, R2+1, color=curve_styles[rmax][2], lw=curve_styles[rmax][0], alpha=curve_styles[rmax][1], ls=curve_styles[rmax][4])
        plt.plot([0.726, 0.726], [1.04, 1.19], 'red', lw=5)
        lines[rmax], = plt.plot(x, R3+2, color=curve_styles[rmax][2], lw=curve_styles[rmax][0], alpha=curve_styles[rmax][1], ls=curve_styles[rmax][4])
        plt.plot([0.872, 0.872], [2.1, 2.25], 'red', lw=5)
        plt.xlabel(r'sin $\theta$')
        plt.ylabel(r'$\mathbf{R_{spec}}$', rotation=0)
        plt.gca().yaxis.set_label_coords(-0.1, 0.5)
    plt.legend(handles=[lines['7'], lines['12'], lines['14']], labels=['RMAX=7', 'RMAX=12', 'RMAX=14'], frameon=False)
    plt.tick_params('both', length=10, width=2, which='major')

    fig.text(0.5, 0.23, r'$\lambda=650$ nm', ha='center')
    fig.text(0.5, 0.48, r'$\lambda=750$ nm', ha='center')
    fig.text(0.5, 0.8, r'$\lambda=900$ nm', ha='center')
    
    plt.savefig('fig5.pdf')
    plt.clf(); plt.close()




if 'fig9' in figures_to_plot:
    from decimal import Decimal
    LMAX = [10]
    RMAX = [16]
    AK1 = [1e-3, 8e-5, 8e-6]
    X = {}; Y = {}
    step = 1
    #reading data to plot
    for version in ['3_cerf', '3']:
        main_dir = dir=f'data/fig9/{version}/'
        for lmax in LMAX:
            for rmax in RMAX:
                for ak1 in AK1:
                    ak1 = round(ak1/2/np.pi, 8)
                    path = main_dir+f'lmax={lmax}_rmax={rmax}_ak1={ak1}.txt'
                    x, y = cl.read_1D_data_real(path)
                    X[f'{version}_{rmax}_{ak1}'] = x; Y[f'{version}_{rmax}_{ak1}'] = y
    
    fig = plt.figure(figsize=(22, 22))
    error_step = 5
    axs = []
    for i, ak1 in enumerate(AK1):
        ak1 = round(ak1/2/np.pi, 8) 
        ax = fig.add_subplot(3, 2, 2*(i+1)-1)
        axs.append(ax)
        key = f'{RMAX[0]}_{ak1}'
        plt.gca().plot(X['3_cerf_'+key], Y['3_cerf_'+key], 'green', lw=18, alpha=0.4, label='CERF') 
        plt.gca().plot(X['3_'+key], Y['3_'+key], 'black', lw=6, label='SciPy')
        plt.gca().set_title(r'$k_x$='+f'{Decimal(AK1[i]):.0e}')
        plt.gca().set_xlabel(r'$ak_0$')
        plt.gca().legend(frameon=False)
        plt.gca().tick_params('both', length=15, width=4, which='major')
        plt.gca().minorticks_off()
        plt.gca().set_ylabel(r'$\mathbf{R}$', rotation=0)
        plt.gca().yaxis.set_label_coords(-0.18, 0.5)
        ax = fig.add_subplot(3, 2, 2*(i+1))
        axs.append(ax)
        error = np.abs(Y['3_'+key]-Y['3_cerf_'+key])
        plt.gca().scatter(X['3_'+key][::error_step], error[::error_step], color='black', s=20)
        plt.gca().set_ylabel('error')
        plt.gca().set_xlabel(r'$ak_0$')
        plt.gca().set_yscale('log')
        plt.gca().tick_params('both', length=15, width=4, which='major')
        plt.gca().minorticks_off()
        plt.gca().yaxis.set_label_coords(-0.18, 0.5)



    fig.text(0.04, 0.96, r'a', ha='center', fontsize=30)
    fig.text(0.54, 0.96, r'b', ha='center', fontsize=30)
    fig.text(0.04, 0.64, r'c', ha='center', fontsize=30)
    fig.text(0.54, 0.64, r'd', ha='center', fontsize=30)
    fig.text(0.04, 0.33, r'e', ha='center', fontsize=30)
    fig.text(0.54, 0.33, r'f', ha='center', fontsize=30)



    # #setting ticks and labels
    axs[0].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    axs[0].set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
    axs[0].set_xticks([3.7314, 3.7316, 3.7318, 3.7320])
    axs[0].set_xticklabels([3.7314, 3.7316, 3.7318, 3.7320])
    axs[1].set_xticks([3.7314, 3.7316, 3.7318, 3.7320])
    axs[1].set_xticklabels([3.7314, 3.7316, 3.7318, 3.7320])
    axs[2].set_xticks([3.7316, 3.73163, 3.73166, 3.73169])
    axs[2].set_xticklabels([3.7316, 3.73163, 3.73166, 3.73169])
    axs[3].set_xticks([3.7316, 3.73163, 3.73166, 3.73169])
    axs[3].set_xticklabels([3.7316, 3.73163, 3.73166, 3.73169])
    axs[4].set_xticks([3.73163055, 3.731630575])
    axs[4].set_xticklabels([3.73163055, 3.731630575])
    axs[5].set_xticks([3.73163055, 3.731630575])
    axs[5].set_xticklabels([3.73163055, 3.731630575])


    # axs[2].set_xticks([3.7354, 3.7356, 3.7358])
    # axs[2].set_xticklabels([3.7354, 3.7356, 3.7358])
    # axs[2].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    # axs[2].set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
    # axs[3].set_xticks([3.7354, 3.7356, 3.7358])
    # axs[3].set_xticklabels([3.7354, 3.7356, 3.7358])
    # axs[4].set_xticks([3.73558, 3.735583, 3.735586])
    # axs[4].set_xticklabels([3.73558, 3.735583, 3.735586])
    # axs[4].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    # axs[4].set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
    # axs[5].set_xticks([3.73558, 3.735583, 3.735586])
    # axs[5].set_xticklabels([3.73558, 3.735583, 3.735586])
    
  
    plt.subplots_adjust(left=0.12,
                        bottom=0.099,
                        right=0.98,
                        top=0.95,
                        wspace=0.3,
                        hspace=0.4)
    plt.savefig(f'fig9.pdf')
    plt.clf(); plt.close()


if 'fig8' in figures_to_plot:
    from decimal import Decimal
    LMAX = [4]
    RMAX = [16]
    AK1 = [1e-3, 8e-5, 8e-6]
    X = {}; Y = {}
    step = 1
    #reading data to plot
    for version in ['3_cerf', '3']:
        main_dir = dir=f'data/fig8/{version}/'
        for lmax in LMAX:
            for rmax in RMAX:
                for ak1 in AK1:
                    ak1 = round(ak1/2/np.pi, 8)
                    path = main_dir+f'lmax={lmax}_rmax={rmax}_ak1={ak1}.txt'
                    x, y = cl.read_1D_data_real(path)
                    X[f'{version}_{rmax}_{ak1}'] = x; Y[f'{version}_{rmax}_{ak1}'] = y
    
    fig = plt.figure(figsize=(22, 22))
    error_step = 5
    axs = []
    for i, ak1 in enumerate(AK1):
        ak1 = round(ak1/2/np.pi, 8) 
        ax = fig.add_subplot(3, 2, 2*(i+1)-1)
        axs.append(ax)
        key = f'{RMAX[0]}_{ak1}'
        plt.gca().plot(X['3_cerf_'+key], Y['3_cerf_'+key], 'green', lw=18, alpha=0.4, label='CERF') 
        plt.gca().plot(X['3_'+key], Y['3_'+key], 'black', lw=6, label='SciPy')
        plt.gca().set_title(r'$k_x$='+f'{Decimal(AK1[i]):.0e}')
        plt.gca().set_xlabel(r'$ak_0$')
        plt.gca().legend(frameon=False)
        plt.gca().tick_params('both', length=15, width=4, which='major')
        plt.gca().minorticks_off()
        plt.gca().set_ylabel(r'$\mathbf{R}$', rotation=0)
        plt.gca().yaxis.set_label_coords(-0.18, 0.5)
        ax = fig.add_subplot(3, 2, 2*(i+1))
        axs.append(ax)
        error = np.abs(Y['3_'+key]-Y['3_cerf_'+key])
        plt.gca().scatter(X['3_'+key][::error_step], error[::error_step], color='black', s=20)
        plt.gca().set_ylabel('error')
        plt.gca().set_xlabel(r'$ak_0$')
        plt.gca().set_yscale('log')
        plt.gca().tick_params('both', length=15, width=4, which='major')
        plt.gca().minorticks_off()
        plt.gca().yaxis.set_label_coords(-0.18, 0.5)



    fig.text(0.04, 0.96, r'a', ha='center', fontsize=30)
    fig.text(0.54, 0.96, r'b', ha='center', fontsize=30)
    fig.text(0.04, 0.64, r'c', ha='center', fontsize=30)
    fig.text(0.54, 0.64, r'd', ha='center', fontsize=30)
    fig.text(0.04, 0.33, r'e', ha='center', fontsize=30)
    fig.text(0.54, 0.33, r'f', ha='center', fontsize=30)



    # #setting ticks and labels
    # axs[0].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    # axs[0].set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
    # axs[0].set_xticks([3.7314, 3.7316, 3.7318, 3.7320])
    # axs[0].set_xticklabels([3.7314, 3.7316, 3.7318, 3.7320])
    # axs[1].set_xticks([3.7314, 3.7316, 3.7318, 3.7320])
    # axs[1].set_xticklabels([3.7314, 3.7316, 3.7318, 3.7320])
    # axs[2].set_xticks([3.7316, 3.73163, 3.73166, 3.73169])
    # axs[2].set_xticklabels([3.7316, 3.73163, 3.73166, 3.73169])
    # axs[3].set_xticks([3.7316, 3.73163, 3.73166, 3.73169])
    # axs[3].set_xticklabels([3.7316, 3.73163, 3.73166, 3.73169])
    # axs[4].set_xticks([3.73163055, 3.731630575])
    # axs[4].set_xticklabels([3.73163055, 3.731630575])
    # axs[5].set_xticks([3.73163055, 3.731630575])
    # axs[5].set_xticklabels([3.73163055, 3.731630575])


    axs[0].set_xticks([3.7354, 3.7355, 3.7356, 3.7357])
    axs[0].set_xticklabels([3.7354, 3.7355, 3.7356, 3.7357])
    axs[0].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    axs[0].set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
    axs[1].set_xticks([3.7354, 3.7355, 3.7356, 3.7357])
    axs[1].set_xticklabels([3.7354, 3.7355, 3.7356, 3.7357])
    axs[2].set_xticks([3.7354, 3.7355, 3.7356, 3.7357])
    axs[2].set_xticklabels([3.7354, 3.7355, 3.7356, 3.7357])
    axs[2].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    axs[2].set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
    axs[3].set_xticks([3.7354, 3.7355, 3.7356, 3.7357])
    axs[3].set_xticklabels([3.7354, 3.7355, 3.7356, 3.7357])
    axs[4].set_xticks([3.73553338, 3.735533395])
    axs[4].set_xticklabels([3.73553338, 3.735533395])
    axs[4].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    axs[4].set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
    axs[5].set_xticks([3.73553338, 3.735533395])
    axs[5].set_xticklabels([3.73553338, 3.735533395])

    
  
    plt.subplots_adjust(left=0.12,
                        bottom=0.099,
                        right=0.98,
                        top=0.95,
                        wspace=0.3,
                        hspace=0.4)
    plt.savefig(f'fig8.pdf')
    plt.clf(); plt.close()



if 'fig6_old' in figures_to_plot:
    import os
    from matplotlib.offsetbox import OffsetImage, AnnotationBbox
    lmax = 7
    fig = plt.figure(figsize=(20, 20))
    X = {}; Y = {}
    d = '0.41'
    nunit = '_1unit_6'
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='red', lw=18, alpha=0.4, ls='solid'),
            Line2D([0], [0], color='red', lw=6, alpha=1, ls='dotted')]


    ax1 = fig.add_subplot(2, 1, 1)
   # #key:version_rmax
    # #value: (lw, alpha, color, zorder, styles)
    curve_styles = {
        '2_20': (9, 1, 'red', 1, 'dotted'),
        '2_30': (9, 1, 'black', 1, 'dashed'),
        '2_26': (18, 0.7, 'green', 1, 'solid'),
        '2_40': (6, 1, 'blue', 1, 'dotted'),
        # '2_45': (6, 1, 'black', 1, '-'),
    }
    for version in ['2']:
        main_dir = f'data/fig6_new{nunit}/mode=3_old/version={version}/d={d}/'
        for rmax in ['20']:
            path = main_dir+f'lmax={lmax}_rmax={rmax}.txt'
            if os.path.exists(path): 
                x, y = cl.read_1D_data_real(path)
            else:
                print('File not Found')
                continue
            X[f'{version}_{rmax}'] = x; Y[f'{version}_{rmax}'] = y
            ax1.plot(X[f'{version}_{rmax}'], Y[f'{version}_{rmax}'], label=f'RMAX={rmax}',
                    lw=curve_styles[f'{version}_{rmax}'][0], 
                    alpha=curve_styles[f'{version}_{rmax}'][1],
                    color=curve_styles[f'{version}_{rmax}'][2],
                    zorder=curve_styles[f'{version}_{rmax}'][3],
                    ls=curve_styles[f'{version}_{rmax}'][4])

    ax1.plot((x[0], x[-1]),(1.0, 1.0), '-', color='blue', lw=3)
    ax1.legend(frameon=False)
    # ax1.set_ylim([-0.05, 0.5])
    fig.text(0.5, 0.6, r'MULTEM2', color='black', ha='center', fontsize=30)
    ax1.set_ylabel(r'$\mathbf{R}$', rotation=0, fontsize=30)
    # ax2.yaxis.set_label_coords(-0.08, 0.5)        


    ax2 = fig.add_subplot(2, 1, 2)
    # X = {}; Y = {}
    # #key:version_rmax
    # #value: (lw, alpha, color, zorder)
    curve_styles = {
        '3_20': (9, 1, 'red', 1, 'dotted'),
        '3_30': (18, 0.7, 'green', 1, 'solid'),
        # '3_26': (9, 0.7, 'green', 1, 'dashed'),
        '3_40': (6, 1, 'red', 1, 'dotted'),
        '3_22': (18, 0.4, 'green', 1, '-'),
        # '2_45': (6, 1, 'black', 1, '-'),
    }
    for version in ['3']:
        main_dir = f'data/fig6_new{nunit}/mode=3_old/version={version}/d={d}/'
        for rmax in ['20', '30', '40']:
        # for rmax in ['47', '45']:
            path = main_dir+f'lmax={lmax}_rmax={rmax}.txt'
            if os.path.exists(path): 
                x, y = cl.read_1D_data_real(path)
            else:
                print('File not Found')
                continue
            X[f'{version}_{rmax}'] = x; Y[f'{version}_{rmax}'] = y
            ax2.plot(X[f'{version}_{rmax}'], Y[f'{version}_{rmax}'], label=f'RMAX={rmax}',
                    lw=curve_styles[f'{version}_{rmax}'][0], 
                    alpha=curve_styles[f'{version}_{rmax}'][1],
                    color=curve_styles[f'{version}_{rmax}'][2],
                    zorder=curve_styles[f'{version}_{rmax}'][3],
                    ls=curve_styles[f'{version}_{rmax}'][4])


    # ax3 = fig.add_subplot(2, 1, 2)
    # error_step = 1
    # error = np.abs(Y['3_20']-Y['2_20'])
    # ax3.scatter(X['3_20'][::error_step], error[::error_step], color='black', s=20)
    # print(Y[f'2_50']-Y[f'3_50'])
    # print(Y[f'2_20']-Y[f'2_30'])
    # print(Y[f'2_16']-Y[f'3_16'])

    # print(Y[f'2_16']-Y[f'2_35'])

    print('----------------------')
    # print(np.max(Y[f'2_16']))
    # print(np.max(Y[f'2_35']))
    # print(np.max(Y[f'2_50']))
    # print(np.max(Y[f'3_16']))
    # print(np.max(Y[f'3_35']))



    # ax3.set_xticks([6.29, 6.291, 6.291])
    # ax3.set_xticklabels([6.29, 6.291, 6.291])
    ax2.plot((x[0], x[-1]),(1.0, 1.0), '-', color='blue', lw=3)
    
    # ax3.set_yticks([0.0, 0.5, 1.0, 1.5])
    # ax3.set_yticklabels([0.0, 0.5, 1.0, 1.5])
    # ax2.tick_params('both', length=7, width=2, which='both')
    # ax2.legend(frameon=False)
    fig.text(0.5, 0.1, r'MULTEM3', color='black', ha='center', fontsize=30)

    # ax2.set_ylim([-0.1, 0.8])
    ax2.set_xlabel(r'$ak_0$')
    ax2.set_ylabel(r'$\mathbf{R}$', rotation=0, fontsize=30)
    # ax3.yaxis.set_label_coords(-0.08, 0.5)
    # arr_img = plt.imread("unit_slice.png")
    # im = OffsetImage(arr_img, zoom=0.55)
    # ab = AnnotationBbox(im, (0.86, 0.74), xycoords='axes fraction', frameon=False)
    # ax3.add_artist(ab)


    # fig.text(0.03, 0.96, 'a', color='black', ha='center', fontsize=30)
    # fig.text(0.03, 0.65, 'b', color='black', ha='center', fontsize=30)
    # fig.text(0.03, 0.33, 'c', color='black', ha='center', fontsize=30)
    # fig.text(0.7, 0.33, 'd', color='black', ha='center', fontsize=30)

    plt.subplots_adjust(left=0.1,
                        bottom=0.05,
                        right=0.95,
                        top=0.98,
                        wspace=0.01,
                        hspace=0.1)


    # for ax in [ax1, ax2, ax3]:
    #     ax.tick_params('both', length=10, width=2, which='major')

    plt.savefig(f'fig6.pdf')
    plt.clf(); plt.close()


if 'fig6_test' in figures_to_plot:
    import os
    from matplotlib.offsetbox import OffsetImage, AnnotationBbox
    lmax = 7
    fig = plt.figure(figsize=(20, 20))
    X = {}; Y = {}
    d = '0.41'
    # nunit = '_1unit_9nlayer_2nplan'
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='red', lw=18, alpha=0.4, ls='solid'),
            Line2D([0], [0], color='red', lw=6, alpha=1, ls='dotted')]

    # nunit = '_1unit_6'
    ax1 = fig.add_subplot(3, 1, 1)
   # #key:version_rmax
    # #value: (lw, alpha, color, zorder, styles)
    curve_styles = {
        '2_32_5': (3, 1, 'green', 1, 'solid'),
        '2_16_6': (9, 1, 'red', 1, 'dashed'),
        '2_30_3': (3, 1, 'black', 1, 'solid'),
        '2_30_9': (6, 1, 'blue', 1, 'dotted'),
        # '2_45': (6, 1, 'black', 1, '-'),
    }
    NLAYERS = ['6']
    for version in ['2']:
        for nlayer in NLAYERS:
                nunit = f'_1unit_{nlayer}'
                for rmax in ['16']:
                    main_dir = f'data/fig6_new{nunit}/mode=3_ref/version={version}/d={d}/'
                    path = main_dir+f'lmax={lmax}_rmax={rmax}.txt'
                    if os.path.exists(path): 
                        x, y = cl.read_1D_data_real(path)
                    else:
                        print('File not Found')
                        continue
                    X[f'{version}_{rmax}_{nlayer}'] = x; Y[f'{version}_{rmax}_{nlayer}'] = y
                    ax1.plot(X[f'{version}_{rmax}_{nlayer}'], Y[f'{version}_{rmax}_{nlayer}'], label=f'nlayer={nlayer}',
                            lw=curve_styles[f'{version}_{rmax}_{nlayer}'][0], 
                            alpha=curve_styles[f'{version}_{rmax}_{nlayer}'][1],
                            color=curve_styles[f'{version}_{rmax}_{nlayer}'][2],
                            zorder=curve_styles[f'{version}_{rmax}_{nlayer}'][3],
                            ls=curve_styles[f'{version}_{rmax}_{nlayer}'][4])

    ax1.plot((x[0], x[-1]),(1.0, 1.0), '-', color='blue', lw=3)
    ax1.legend(frameon=False)
    ax1.set_ylim([-0.05, 1.5])
    fig.text(0.5, 0.85, r'MULTEM2', color='black', ha='center', fontsize=30)
    ax1.set_ylabel(r'$\mathbf{R}$', rotation=0, fontsize=30)
    # ax2.yaxis.set_label_coords(-0.08, 0.5)        


    ax2 = fig.add_subplot(3, 1, 2)
    # X = {}; Y = {}
    # #key:version_rmax
    # #value: (lw, alpha, color, zorder)
    curve_styles = {
        '3_32_5': (3, 1, 'green', 1, 'solid'),
        '3_30_6': (9, 1, 'red', 1, 'dashed'),
        '3_30_3': (3, 1, 'black', 1, 'solid'),
        '3_30_9': (6, 1, 'blue', 1, 'dotted'),
        # '2_45': (6, 1, 'black', 1, '-'),
    }
    for version in ['3']:
        for nlayer in NLAYERS:
            nunit = f'_1unit_{nlayer}'
            for rmax in ['30']:
            # for rmax in ['47', '45']:
                main_dir = f'data/fig6_new{nunit}/mode=3_ref/version={version}/d={d}/'
                path = main_dir+f'lmax={lmax}_rmax={rmax}.txt'
                if os.path.exists(path): 
                    x, y = cl.read_1D_data_real(path)
                else:
                    print('File not Found')
                    continue
                X[f'{version}_{rmax}_{nlayer}'] = x; Y[f'{version}_{rmax}_{nlayer}'] = y
                ax2.plot(X[f'{version}_{rmax}_{nlayer}'], Y[f'{version}_{rmax}_{nlayer}'], label=f'nlayer={nlayer}',
                        lw=curve_styles[f'{version}_{rmax}_{nlayer}'][0], 
                        alpha=curve_styles[f'{version}_{rmax}_{nlayer}'][1],
                        color=curve_styles[f'{version}_{rmax}_{nlayer}'][2],
                        zorder=curve_styles[f'{version}_{rmax}_{nlayer}'][3],
                        ls=curve_styles[f'{version}_{rmax}_{nlayer}'][4])

    ax2.set_ylim([-0.05, 1.5])
    ax2.plot((x[0], x[-1]),(1.0, 1.0), '-', color='blue', lw=3)
    ax2.legend(frameon=False)

    # ax3 = fig.add_subplot(3, 1, 3)
    # error_step = 1
    # error = np.abs(Y['3_30']-Y['3_32'])
    # ax3.scatter(X['2_32'][::error_step], error[::error_step], color='black', s=20)
    # ax3.set_yscale('log')
    # ax3.set_ylabel(r'$error$', rotation=0, fontsize=30)
    # print(Y[f'3_20']-Y[f'3_30'])
    # print(Y[f'2_20']-Y[f'2_30'])
    # print(Y[f'2_16']-Y[f'3_16'])

    # print(Y[f'2_16']-Y[f'2_35'])

    print('----------------------')
    # print(np.max(Y[f'2_16']))
    # print(np.max(Y[f'2_35']))
    # print(np.max(Y[f'2_50']))
    # print(np.max(Y[f'3_16']))
    # print(np.max(Y[f'3_35']))


#     nunit = '_1unit_8'
#     ax1 = fig.add_subplot(3, 2, 2)
#    # #key:version_rmax
#     # #value: (lw, alpha, color, zorder, styles)
#     curve_styles = {
#         '2_20': (9, 1, 'red', 1, 'solid'),
#         '2_30': (12, 0.7, 'black', 1, 'solid'),
#         '2_6': (18, 0.7, 'green', 1, 'solid'),
#         '2_40': (6, 1, 'blue', 1, 'dotted'),
#         # '2_45': (6, 1, 'black', 1, '-'),
#     }
#     for version in ['2']:
#         main_dir = f'data/fig6{nunit}/mode=3_ref/version={version}/d={d}/'
#         for rmax in ['20', '30', '40']:
#             path = main_dir+f'lmax={lmax}_rmax={rmax}.txt'
#             if os.path.exists(path): 
#                 x, y = cl.read_1D_data_real(path)
#             else:
#                 print('File not Found')
#                 continue
#             X[f'{version}_{rmax}'] = x; Y[f'{version}_{rmax}'] = y
#             ax1.plot(X[f'{version}_{rmax}'], Y[f'{version}_{rmax}'], label=f'RMAX={rmax}',
#                     lw=curve_styles[f'{version}_{rmax}'][0], 
#                     alpha=curve_styles[f'{version}_{rmax}'][1],
#                     color=curve_styles[f'{version}_{rmax}'][2],
#                     zorder=curve_styles[f'{version}_{rmax}'][3],
#                     ls=curve_styles[f'{version}_{rmax}'][4])

#     ax1.plot((x[0], x[-1]),(1.0, 1.0), '-', color='blue', lw=3)
#     ax1.legend(frameon=False)
#     ax1.set_ylim([-0.05, 1.5])
#     fig.text(0.5, 0.85, r'MULTEM2', color='black', ha='center', fontsize=30)
#     ax1.set_ylabel(r'$\mathbf{R}$', rotation=0, fontsize=30)
#     # ax2.yaxis.set_label_coords(-0.08, 0.5)        


#     ax2 = fig.add_subplot(3, 2, 4)
#     # X = {}; Y = {}
#     # #key:version_rmax
#     # #value: (lw, alpha, color, zorder)
#     curve_styles = {
#         '3_20': (9, 1, 'red', 1, 'solid'),
#         '3_30': (12, 0.7, 'black', 1, 'solid'),
#         # '3_26': (9, 0.7, 'green', 1, 'dashed'),
#         '3_40': (6, 1, 'red', 1, 'dotted'),
#         '3_22': (18, 0.4, 'green', 1, '-'),
#         # '2_45': (6, 1, 'black', 1, '-'),
#     }
#     for version in ['3']:
#         main_dir = f'data/fig6{nunit}/mode=3_ref/version={version}/d={d}/'
#         for rmax in ['20', '30', '40']:
#         # for rmax in ['47', '45']:
#             path = main_dir+f'lmax={lmax}_rmax={rmax}.txt'
#             if os.path.exists(path): 
#                 x, y = cl.read_1D_data_real(path)
#             else:
#                 print('File not Found')
#                 continue
#             X[f'{version}_{rmax}'] = x; Y[f'{version}_{rmax}'] = y
#             ax2.plot(X[f'{version}_{rmax}'], Y[f'{version}_{rmax}'], label=f'RMAX={rmax}',
#                     lw=curve_styles[f'{version}_{rmax}'][0], 
#                     alpha=curve_styles[f'{version}_{rmax}'][1],
#                     color=curve_styles[f'{version}_{rmax}'][2],
#                     zorder=curve_styles[f'{version}_{rmax}'][3],
#                     ls=curve_styles[f'{version}_{rmax}'][4])

#     ax2.set_ylim([-0.05, 1.5])
#     ax2.plot((x[0], x[-1]),(1.0, 1.0), '-', color='blue', lw=3)
#     ax2.legend(frameon=False)

#     ax3 = fig.add_subplot(3, 2, 6)
#     error_step = 1
#     error = np.abs(Y['2_40']-Y['3_40'])
#     ax3.scatter(X['2_40'][::error_step], error[::error_step], color='black', s=20)
#     ax3.set_yscale('log')
#     ax3.set_ylabel(r'$error$', rotation=0, fontsize=30)
#     print(Y[f'2_30']-Y[f'2_40'])
#     # print(Y[f'2_20']-Y[f'2_30'])
#     # print(Y[f'2_16']-Y[f'3_16'])

#     # print(Y[f'2_16']-Y[f'2_35'])

#     print('----------------------')




#     # ax3.set_xticks([6.29, 6.291, 6.291])
#     # ax3.set_xticklabels([6.29, 6.291, 6.291])
#     # ax2.plot((x[0], x[-1]),(1.0, 1.0), '-', color='blue', lw=3)
    
#     # ax3.set_yticks([0.0, 0.5, 1.0, 1.5])
#     # ax3.set_yticklabels([0.0, 0.5, 1.0, 1.5])
#     # ax2.tick_params('both', length=7, width=2, which='both')
#     # ax2.legend(frameon=False)
#     fig.text(0.5, 0.5, r'MULTEM3', color='black', ha='center', fontsize=30)

#     # ax2.set_ylim([-0.1, 0.8])
#     ax2.set_xlabel(r'$ak_0$')
#     ax2.set_ylabel(r'$\mathbf{R}$', rotation=0, fontsize=30)
    # ax3.yaxis.set_label_coords(-0.08, 0.5)
    # arr_img = plt.imread("unit_slice.png")
    # im = OffsetImage(arr_img, zoom=0.55)
    # ab = AnnotationBbox(im, (0.86, 0.74), xycoords='axes fraction', frameon=False)
    # ax3.add_artist(ab)


    # fig.text(0.03, 0.96, 'a', color='black', ha='center', fontsize=30)
    # fig.text(0.03, 0.65, 'b', color='black', ha='center', fontsize=30)
    # fig.text(0.03, 0.33, 'c', color='black', ha='center', fontsize=30)
    # fig.text(0.7, 0.33, 'd', color='black', ha='center', fontsize=30)

    plt.subplots_adjust(left=0.1,
                        bottom=0.05,
                        right=0.95,
                        top=0.98,
                        wspace=0.01,
                        hspace=0.1)


    # for ax in [ax1, ax2, ax3]:
    #     ax.tick_params('both', length=10, width=2, which='major')

    plt.savefig(f'fig6_test.pdf')
    plt.clf(); plt.close()


if 'fig6' in figures_to_plot:
    import os
    from matplotlib.offsetbox import OffsetImage, AnnotationBbox
    lmax = 7
    fig = plt.figure(figsize=(20, 20))
    X = {}; Y = {}
    d = '0.41'
    # nunit = '_1unit_9nlayer_2nplan'
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='red', lw=18, alpha=0.4, ls='solid'),
            Line2D([0], [0], color='red', lw=6, alpha=1, ls='dotted')]

    ax1 = fig.add_subplot(3, 1, 1)
   # #key:version_rmax
    # #value: (lw, alpha, color, zorder, styles)
    curve_styles = {
        '2_32_3': (9, 1, 'green', 1, 'solid'),
        '2_32_4': (9, 1, 'red', 1, 'dashed'),
        '2_32_6': (3, 1, 'black', 1, 'solid'),
        '2_40': (6, 1, 'blue', 1, 'dotted'),
        # '2_45': (6, 1, 'black', 1, '-'),
    }
    NLAYERS = ['3', '4', '5', '6', '7', '8', '9'][::-1]
    for version in ['2']:
        for rmax in ['20', '32']:
            for nlayer in NLAYERS:
                nunit = f'_1unit_{nlayer}'
                main_dir = f'data/fig6_new{nunit}/mode=3_ref/version={version}/d={d}/'
                path = main_dir+f'lmax={lmax}_rmax={rmax}.txt'
                if os.path.exists(path): 
                    x, y = cl.read_1D_data_real(path)
                else:
                    print('File not Found')
                    continue
                X[f'{version}_{rmax}_{nlayer}'] = x; Y[f'{version}_{rmax}_{nlayer}'] = y
                if rmax == '32' and (nlayer == '3' or nlayer == '4' or nlayer == '6'):
                    ax1.plot(X[f'{version}_{rmax}_{nlayer}'], Y[f'{version}_{rmax}_{nlayer}'], label=f'{2**(int(nlayer)-1)}',
                            lw=curve_styles[f'{version}_{rmax}_{nlayer}'][0], 
                            alpha=curve_styles[f'{version}_{rmax}_{nlayer}'][1],
                            color=curve_styles[f'{version}_{rmax}_{nlayer}'][2],
                            zorder=curve_styles[f'{version}_{rmax}_{nlayer}'][3],
                            ls=curve_styles[f'{version}_{rmax}_{nlayer}'][4])

    ax1.plot((x[0], x[-1]),(1.0, 1.0), '-', color='blue', lw=5, alpha=0.6)
    fig.text(0.2, 0.95, 'planes', color='black', ha='center', fontsize=30)
    ax1.legend(frameon=False, bbox_to_anchor=(0.06,0.5))
    ax1.set_ylim([-0.05, 2.5])
    fig.text(0.5, 0.95, r'MULTEM2', color='black', ha='center', fontsize=30)
    ax1.set_xlabel(r'$ak_0$')
    ax1.set_ylabel(r'$\mathbf{R}$', rotation=0, fontsize=30)
    ax1.yaxis.set_label_coords(-0.08, 0.5)        

    # ax1.set_xticks([0.8, 0.9, 1.0, 1.1, 1.2])
    # ax1.set_xticklabels([0.8, 0.9, 1.0, 1.1, 1.2])


    ax2 = fig.add_subplot(3, 1, 2)
    # X = {}; Y = {}
    # #key:version_rmax
    # #value: (lw, alpha, color, zorder)
    curve_styles = {
        '3_32_3': (9, 1, 'green', 1, 'solid'),
        '3_32_4': (9, 1, 'red', 1, 'dashed'),
        '3_32_6': (3, 1, 'black', 1, 'solid'),
        '3_40': (6, 1, 'red', 1, 'dotted'),
        '3_22': (18, 0.4, 'green', 1, '-'),
        # '2_45': (6, 1, 'black', 1, '-'),
    }
    for version in ['3']:
        for rmax in ['20', '32']:
             for nlayer in NLAYERS:
                nunit = f'_1unit_{nlayer}'
                main_dir = f'data/fig6_new{nunit}/mode=3_ref/version={version}/d={d}/'
                path = main_dir+f'lmax={lmax}_rmax={rmax}.txt'
                if os.path.exists(path): 
                    x, y = cl.read_1D_data_real(path)
                else:
                    print('File not Found')
                    continue
                X[f'{version}_{rmax}_{nlayer}'] = x; Y[f'{version}_{rmax}_{nlayer}'] = y
                if rmax == '32' and (nlayer == '3' or nlayer == '4' or nlayer == '6'):
                    ax2.plot(X[f'{version}_{rmax}_{nlayer}'], Y[f'{version}_{rmax}_{nlayer}'], label=f'{2**(int(nlayer)-1)}',
                            lw=curve_styles[f'{version}_{rmax}_{nlayer}'][0], 
                            alpha=curve_styles[f'{version}_{rmax}_{nlayer}'][1],
                            color=curve_styles[f'{version}_{rmax}_{nlayer}'][2],
                            zorder=curve_styles[f'{version}_{rmax}_{nlayer}'][3],
                            ls=curve_styles[f'{version}_{rmax}_{nlayer}'][4])

    ax2.set_ylim([-0.05, 2.5])
    ax2.plot((x[0], x[-1]),(1.0, 1.0), '-', color='blue',  lw=5, alpha=0.6)
    fig.text(0.2, 0.62, 'planes', color='black', ha='center', fontsize=30)
    ax2.legend(frameon=False, bbox_to_anchor=(0.06,0.5))
    # ax2.set_xticks([0.8, 0.9, 1.0, 1.1, 1.2])
    # ax2.set_xticklabels([0.8, 0.9, 1.0, 1.1, 1.2])

    #TODO automatic error generation
    ax3 = fig.add_subplot(3, 1, 3)
    # error_step = 1 
    # error_rmax20 = [np.max(np.abs(Y['2_20_3']-Y['3_20_3'])), np.max(np.abs(Y['2_20_4']-Y['3_20_4'])), np.max(np.abs(Y['2_20_5']-Y['3_20_5'])), 
    #                 np.max(np.abs(Y['2_20_6']-Y['3_20_6'])), np.max(np.abs(Y['2_20_7']-Y['3_20_7'])), np.max(np.abs(Y['2_20_8']-Y['3_20_8'])), np.max(np.abs(Y['2_20_9']-Y['3_20_9']))]
    # ax3.plot([2**(3-1), 2**(4-1), 2**(5-1),  2**(6-1),  2**(7-1),  2**(8-1),  2**(9-1)], error_rmax20, color='blue', marker='s', ms=20, ls='dashed', label="20", lw=5)
    # # error_rmax30 = [np.max(np.abs(Y['2_30_3']-Y['3_30_3'])), np.max(np.abs(Y['2_30_4']-Y['3_30_4'])), np.max(np.abs(Y['2_30_5']-Y['3_30_5'])), 
    # #                 np.max(np.abs(Y['2_30_6']-Y['3_30_6'])), np.max(np.abs(Y['2_30_7']-Y['3_30_7'])), np.max(np.abs(Y['2_30_8']-Y['3_30_8'])), np.max(np.abs(Y['2_30_9']-Y['3_30_9']))]
    # # ax3.plot([2**(3-1), 2**(4-1), 2**(5-1),  2**(6-1),  2**(7-1),  2**(8-1),  2**(9-1)], error_rmax30, color='red', marker='s', ms=20, ls='dashed', label="30", lw=5)
    # error_rmax32 = [np.max(np.abs(Y['2_32_3']-Y['3_32_3'])), np.max(np.abs(Y['2_32_4']-Y['3_32_4'])), np.max(np.abs(Y['2_32_5']-Y['3_32_5'])), 
    #                 np.max(np.abs(Y['2_32_6']-Y['3_32_6'])), np.max(np.abs(Y['2_32_7']-Y['3_32_7'])), np.max(np.abs(Y['2_32_8']-Y['3_32_8'])), np.max(np.abs(Y['2_32_9']-Y['3_32_9']))]
    # ax3.plot([2**(3-1), 2**(4-1), 2**(5-1),  2**(6-1),  2**(7-1),  2**(8-1),  2**(9-1)], error_rmax32, color='red', marker='s', ms=20, ls='dashed', label="32", lw=5)
    # # error_rmax40 = [np.max(np.abs(Y['2_40_3']-Y['3_40_3'])), np.max(np.abs(Y['2_40_4']-Y['3_40_4'])), np.max(np.abs(Y['2_40_5']-Y['3_40_5'])), 
    #                 np.max(np.abs(Y['2_40_6']-Y['3_40_6'])), np.max(np.abs(Y['2_40_7']-Y['3_40_7'])), np.max(np.abs(Y['2_40_8']-Y['3_40_8'])), np.max(np.abs(Y['2_40_9']-Y['3_40_9']))]
    # ax3.plot([2**(3-1), 2**(4-1), 2**(5-1),  2**(6-1),  2**(7-1),  2**(8-1),  2**(9-1)], error_rmax40, color='green', marker='s', ms=20, ls='dashed', label="40", lw=5)

    ax3.set_yscale('log')
    ax3.set_ylabel('difference', fontsize=30)
    ax3.set_xlabel(r'Number of planes', fontsize=30)
    fig.text(0.81, 0.19, 'RMAX', color='black', ha='center', fontsize=30)
    ax3.legend(frameon=False, bbox_to_anchor=(0.9,0.5))



    ax3.set_xticks([4, 32, 64, 128, 256])
    ax3.set_xticklabels([4, 32, 64, 128, 256])
    
 

    fig.text(0.5, 0.62, r'MULTEM3', color='black', ha='center', fontsize=30)

    ax2.set_xlabel(r'$ak_0$')
    ax2.set_ylabel(r'$\mathbf{R}$', rotation=0, fontsize=30)
    ax2.yaxis.set_label_coords(-0.08, 0.5)        


    fig.text(0.03, 0.97, 'a', color='black', ha='center', fontsize=30)
    fig.text(0.03, 0.65, 'b', color='black', ha='center', fontsize=30)
    fig.text(0.03, 0.33, 'c', color='black', ha='center', fontsize=30)

    plt.subplots_adjust(left=0.1,
                        bottom=0.05,
                        right=0.95,
                        top=0.98,
                        wspace=0.01,
                        hspace=0.2)


    for ax in [ax1, ax2, ax3]:
        ax.tick_params('both', length=10, width=2, which='major')
        ax.tick_params('both', length=7, width=1.5, which='minor')


    plt.savefig(f'fig6.pdf')
    plt.clf(); plt.close()

   

if 'fig10' in figures_to_plot:
    curve_styles = {
        '650_m_qsz': (6, 1, 'green', 1, 'dashed'),
        '650_m_qsz_mk': (6, 1, 'green', 1, 'solid'),
        '650_m_total': (20, 0.2, 'green', 1, 'solid'),
        '750_m_qsz': (6, 1, 'green', 1, 'dashed'),
        '750_m_qsz_mk': (6, 1, 'green', 1, 'solid'),
        '750_m_total': (20, 0.2, 'green', 1, 'solid'),
        '750_p_ps': (6, 1, 'red', 1, 'dashed'),
        '750_p_qks_mz_ps': (6, 1, 'red', 1, 'solid'),
        '750_p_total': (20, 0.2, 'red', 1, 'solid'),
        '900_p_ps': (6, 1, 'red', 1, 'dashed'),
        '900_p_qks_mz_ps': (6, 1, 'red', 1, 'solid'),
        '900_p_total': (20, 0.2, 'red', 1, 'solid')
    }
    import os
    fig = plt.figure(figsize=(20, 20))
    main_dir = f'data/fig10/'
    fig.add_subplot(2, 2, 1)
    for dataset in ['650_m_qsz', '650_m_qsz_mk', '650_m_total']:
        path = main_dir+f'{dataset}.txt'
        if os.path.exists(path): 
            x, y = cl.read_1D_data_real(path)
        else:
            print('File not Found')
            continue
        plt.gca().plot(x, y, 
                    lw=curve_styles[dataset][0], 
                    alpha=curve_styles[dataset][1],
                    color=curve_styles[dataset][2],
                    zorder=curve_styles[dataset][3],
                    ls=curve_styles[dataset][4])
    plt.gca().set_ylabel(r'$\mathbf{R_{spec}}$', rotation=0, fontsize=30)
    plt.gca().yaxis.set_label_coords(-0.15, 0.5)
    plt.gca().set_xlabel(r'sin $\theta$')
    plt.gca().tick_params('both', length=15, width=4, which='major')
    plt.gca().set_yticks([0.0, 0.1, 0.3, 0.5, 0.7])
    fig.text(0.32, 0.96, '650 nm', color='black', ha='center', fontsize=30)
    fig.text(0.08, 0.96, 'a', color='black', ha='center', fontsize=30)
    fig.text(0.28, 0.83, r'$q_{sz}$', color='green', ha='center', fontsize=30)
    fig.text(0.33, 0.9, r'$m_{k}+q_{sz}$', color='green', ha='center', fontsize=30)


    fig.add_subplot(2, 2, 2)
    for dataset in ['750_m_qsz', '750_m_qsz_mk', '750_m_total']:
        path = main_dir+f'{dataset}.txt'
        if os.path.exists(path): 
            x, y = cl.read_1D_data_real(path)
        else:
            print('File not Found')
            continue
        plt.gca().plot(x, y, 
                    lw=curve_styles[dataset][0], 
                    alpha=curve_styles[dataset][1],
                    color=curve_styles[dataset][2],
                    zorder=curve_styles[dataset][3],
                    ls=curve_styles[dataset][4])
    plt.gca().set_ylabel(r'$\mathbf{R_{spec}}$', rotation=0, fontsize=30)
    plt.gca().yaxis.set_label_coords(-0.15, 0.5)
    plt.gca().set_xlabel(r'sin $\theta$')
    plt.gca().tick_params('both', length=15, width=4, which='major')
    plt.gca().set_yticks([0.0, 0.1, 0.3, 0.5, 0.7])
    fig.text(0.77, 0.96, '750 nm', color='black', ha='center', fontsize=30)
    fig.text(0.55, 0.96, 'b', color='black', ha='center', fontsize=30)
    fig.text(0.73, 0.86, r'$q_{sz}$', color='green', ha='center', fontsize=30)
    fig.text(0.78, 0.93, r'$m_{k}+q_{sz}$', color='green', ha='center', fontsize=30)


    fig.add_subplot(2, 2, 3)
    for dataset in ['750_p_ps', '750_p_qks_mz_ps', '750_p_total']:
        path = main_dir+f'{dataset}.txt'
        if os.path.exists(path): 
            x, y = cl.read_1D_data_real(path)
        else:
            print('File not Found')
            continue
        plt.gca().plot(x, y, 
                    lw=curve_styles[dataset][0], 
                    alpha=curve_styles[dataset][1],
                    color=curve_styles[dataset][2],
                    zorder=curve_styles[dataset][3],
                    ls=curve_styles[dataset][4])
    plt.gca().set_ylabel(r'$\mathbf{R_{spec}}$', rotation=0, fontsize=30)
    plt.gca().yaxis.set_label_coords(-0.15, 0.5)
    plt.gca().set_xlabel(r'sin $\theta$')
    plt.gca().tick_params('both', length=15, width=4, which='major')
    fig.text(0.32, 0.48, '750 nm', color='black', ha='center', fontsize=30)
    fig.text(0.08, 0.48, 'c', color='black', ha='center', fontsize=30)
    fig.text(0.33, 0.3, r'$p_{s}$', color='red', ha='center', fontsize=30)
    fig.text(0.42, 0.42, r'$p_{s}+m_{z}+q_{ks}$', color='red', ha='center', fontsize=30)

    
    fig.add_subplot(2, 2, 4)
    for dataset in ['900_p_ps', '900_p_qks_mz_ps', '900_p_total']:
        path = main_dir+f'{dataset}.txt'
        if os.path.exists(path): 
            x, y = cl.read_1D_data_real(path)
        else:
            print('File not Found')
            continue
        plt.gca().plot(x, y, 
                    lw=curve_styles[dataset][0], 
                    alpha=curve_styles[dataset][1],
                    color=curve_styles[dataset][2],
                    zorder=curve_styles[dataset][3],
                    ls=curve_styles[dataset][4])
    plt.gca().set_ylabel(r'$\mathbf{R_{spec}}$', rotation=0, fontsize=30)
    plt.gca().yaxis.set_label_coords(-0.15, 0.5)
    plt.gca().set_xlabel(r'sin $\theta$')
    plt.gca().tick_params('both', length=15, width=4, which='major')
    fig.text(0.77, 0.48, '900 nm', color='black', ha='center', fontsize=30)
    fig.text(0.55, 0.48, 'd', color='black', ha='center', fontsize=30)
    fig.text(0.72, 0.41, r'$p_{s}$', color='red', ha='center', fontsize=30)
    fig.text(0.7, 0.2, r'$p_{s}+m_{z}+q_{ks}$', color='red', ha='center', fontsize=30)

    plt.subplots_adjust(left=0.15,
                        bottom=0.1,
                        right=0.95,
                        top=0.95,
                        wspace=0.3,
                        hspace=0.3)



    plt.savefig(f'fig10.pdf')
    plt.clf(); plt.close()



if 'extra_data' in figures_to_plot:
    import collections
    fad_x_re, fad_x_im, fad_y_re, fad_y_im = cl.read_1D_data_complex('fort.3')
    wofad_x_re, wofad_x_im, wofad_y_re, wofad_y_im = cl.read_1D_data_complex('fort.2') 
    fad_args = fad_x_re +1j*fad_x_im
    wofad_args = wofad_x_re +1j*wofad_x_im
    print('Amount of args with abs(z) > 10')
    third_region = np.abs(fad_args)>10
    print(collections.Counter(third_region)[True])
    idx_third_region = np.argwhere(third_region == True)
    print('MAX RE abs(Faddeeva - G.P.M. Poppe, C.M.J. Wijers, Algorithm 680)')
    re_abs_err = np.abs(fad_y_re[idx_third_region]-wofad_y_re[idx_third_region])
    print(np.max(re_abs_err))
    print('MAX IM abs(Faddeeva - G.P.M. Poppe, C.M.J. Wijers, Algorithm 680)')
    im_abs_err = np.abs(fad_y_im[idx_third_region]-wofad_y_im[idx_third_region])
    print(np.max(im_abs_err))

       



        # i += 1
        # plt.gca().plot(omega, np.abs(T['29']-T['7']))
        # plt.gca().set_ylabel('np.abs(T_rmax29 - T_rmax7)')
        # # plt.gca().set_ylabel(r'T')
        # plt.gca().set_xlabel(r'$ak_0$')
        # plt.gca().set_yscale('log')


        # plt.gca().plot(omega, np.abs(T['20']-T['7']))
        # plt.gca().set_title('np.abs(T_rmax20 - T_rmax2)')
        # plt.gca().set_ylabel(r'T')
        # plt.gca().set_xlabel(r'$ak_0$')
        # plt.gca().set_yscale('log')

    # plt.subplot_tool()
    # plt.show()

    # plt.subplots_adjust(left=0.063,
    #                     bottom=0.094,
    #                     right=0.97,
    #                     top=0.93,  plt.subplots_adjust(left=0.1,
                   

    #                     wspace=0.2,
    #                     hspace=0.4)
    




#     fig = plt.figure(figsize=(18, 9))
#     # --- reading data ------
#     data_dir = 'data/fig2a'
#     x, y = np.loadtxt(data_dir+'x.txt'), np.loadtxt(data_dir+'y.txt')
#     x2d = convert_to_square_2d(x)
#     y2d = convert_to_square_2d(y)
#     err_lim = 1e-16
#     multem3_re_err = convert_to_square_2d(load_errors(data_dir+'MIT_re_err.txt', err_lim))
#     multem3_im_err = convert_to_square_2d(load_errors(data_dir+'MIT_im_err.txt', err_lim))
#     multem2_re_err = convert_to_square_2d(load_errors(data_dir+'multem2_re_err.txt', err_lim))
#     multem2_im_err = convert_to_square_2d(load_errors(data_dir+'multem2_im_err.txt', err_lim))
#     cmap_std = matplotlib.colors.ListedColormap(sns.color_palette("RdBu_r", 16))
#     exp_min = -8
#     exp_max = 8
#     implementations = [multem2_re_err, multem2_im_err, multem3_re_err, multem3_im_err]
#     impls_str = ['multem2 re err', 'multem2 im err', 'multem3 re err', 'multem3 im err']

#     # --- reading args for Faddeeva func from MB Dirac cone
#     with open('fort.1') as f:
#         lines = f.readlines()  plt.subplots_adjust(left=0.1,
        
#     args_re, args_im = [], []
#     # ---- postprocess - deleting brackets and left only unique positive values
#     for line in lines:
#         line = line.replace("(","")
#         line = line.replace(")","")
#         re, im = line.split(",")
#         args_re.append(re)
#         args_im.append(im)
#     args_re = np.array(args_re, dtype=np.float64)
#     args_im = np.array(args_im, dtype=np.float64)
#     z = args_re + 1j*args_im
#     Z = np.unique(z)
#     pos_args = []
#     for z in Z:
#         if z.real > 0 and z.imag > 0:
#             pos_args.append(z)
#     pos_args = np.array(pos_args)

#     # ax1
#     ax1 = fig.add_subplot(2,2,1)
#     im = plot_error(ax1, x2d, y2d, multem2_re_err, vextr=[err_lim, 1e-6])
#     ax1.get_xaxis().set_visible(False)
#     ax1.scatter(np.log10(pos_args.real)[::30], np.log10(pos_args.imag)[::30], color='black', s=1)
#     ax1.set_yticks([-8, 0, 8])
#     ax1.set_yticklabels([])
#     # ax2
#     ax2 = fig.add_subplot(2,2,2, sharey=ax1)
#     im = plot_error(ax2, x2d, y2d, multem2_im_err, vextr=[err_lim, 1e-6])
#     ax2.get_yaxis().set_visible(False)
#     ax2.get_xaxis().set_visible(False)
#     # ax2.scatter(np.log10(pos_args.real)[::30], np.log10(pos_args.imag)[::30], color='black', s=1)

#     # ax3
#     ax3 = fig.add_subplot(2,2,3)
#     im = plot_error(ax3, x2d, y2d, multem3_re_err, vextr=[err_lim, 1e-6])
#     ax3.set_yticks([-8, 0, 8])
#     ax3.set_yticklabels([])
#     ax3.set_xticks([-8, 0, 8])
#     ax3.set_xticklabels([])
#     ax3.scatter(np.log10(pos_args.real)[::30], np.log10(pos_args.imag)[::30], color='black', s=1)

#     # ax4
#     ax4 = fig.add_subplot(2,2,4, sharey=ax3)
#     im = plot_error(ax4, x2d, y2d, multem3_im_err, vextr=[err_lim, 1e-6])
#     ax4.get_yaxis().set_visible(False)
#     ax4.set_xticks([-8, 0, 8])
#     ax4.set_xticklabels([])
#     ax4.scatter(np.log10(pos_args.real)[::30], np.log10(pos_args.imag)[::30], color='black', s=1)


#     # ax5 (arg values from MB Dirac cone)
#     ax5 = fig.add_axes([0.34, 0.6, 0.15, 0.35])
#     im = plot_error(ax5, x2d, y2d, multem2_re_err, vextr=[err_lim, 1e-6])
#     ax5.scatter(np.log10(pos_args.real)[::5], np.log10(pos_args.imag)[::5], color='black', s=10)
#     ax5.set_xticks([-0.0092, -0.0074])
#     ax5.set_xticklabels([round(np.min(pos_args.real),3), round(np.max(pos_args.real),3)])
#     ax5.set_yticks([np.min(np.log10(pos_args.imag)), np.max(np.log10(pos_args.imag))])
#     ax5.set_yticks([np.log10(2), np.log10(16)])
#     ax5.set_yticklabels([2, 16])
#     ax5.set_xlim([np.min(np.log10(pos_args.real)) - (np.max(np.log10(pos_args.real)) - np.min(np.log10(pos_args.real)))/2, np.max(np.log10(pos_args.real)) + (np.max(np.log10(pos_args.real)) - np.min(np.log10(pos_args.real)))/2])
#     ax5.set_ylim([np.min(np.log10(pos_args.imag)) - (np.max(np.log10(pos_args.imag)) - np.min(np.log10(pos_args.imag)))/2, np.max(np.log10(pos_args.imag)) + (np.max(np.log10(pos_args.imag)) - np.min(np.log10(pos_args.imag)))/2])

#     for ax_side in ['bottom', 'top', 'left', 'right']:
#         ax5.spines[ax_side].set_color('white')
#         ax5.spines[ax_side].set_linewidth(2)

#     ax1.annotate('', xy=(0.48, 0.55),  xycoords='axes fraction',
#                 xytext=(0.8, 0.56), textcoords='axes fraction',
#                 arrowprops=dict(width=0.1, headwidth=3, edgecolor='white', facecolor='white', shrink=0.1, lw=2),
#                 horizontalalignment='right', verticalalignment='top',
#                 )

#     ax1.add_patch(Rectangle((-0.09, 0.1), 0.16, 1.3, color='white', fill=None, alpha=1, lw=2))

#     # colorbar
#     cb_ax = fig.add_axes([0.915, 0.07, 0.01, 0.895])
#     cbar = fig.colorbar(im, cax=cb_ax)
#     cbar.set_label('relative error', labelpad=10, rotation=90)

#     # annotation and custom axis ticks labeling
#     fig.text(0.495, 0.01, r'$\Re(z)$', ha='center')

#     left_x_pos = 0.02
#     fig.text(left_x_pos, 0.5, r'$\Im(z)$', ha='center', rotation=90)
#     fig.text(left_x_pos, 0.25, 'Multem 3', ha='center', rotation=90)
#     fig.text(left_x_pos, 0.7, 'Multem 2', ha='center', rotation=90)
#     fig.text(0.28, 0.97, r'$\Re(w(z))$', ha='center')
#     fig.text(0.73, 0.97, r'$\Im(w(z))$', ha='center')  plt.subplots_adjust(left=0.1,
                        


#     delta_x_pos = 0.04
#     fig.text(left_x_pos + delta_x_pos, 0.93, r'$10^{8}$', ha='center')
#     fig.text(left_x_pos + delta_x_pos, 0.735, r'$10^{0}$', ha='center')
#     fig.text(left_x_pos + delta_x_pos, 0.54, r'$10^{-8}$', ha='center')
#     fig.text(left_x_pos + delta_x_pos, 0.48, r'$10^{8}$', ha='center')
#     fig.text(left_x_pos + delta_x_pos, 0.28, r'$10^{0}$', ha='center')
#     fig.text(left_x_pos + delta_x_pos, 0.09, r'$10^{-8}$', ha='center')

#     bottom_y_pos = 0.05
#     fig.text(0.08, bottom_y_pos, r'$10^{-8}$', ha='center')
#     fig.text(0.29, bottom_y_pos, r'$10^{0}$', ha='center')
#     fig.text(0.48, bottom_y_pos, r'$10^{8}$', ha='center')
#     fig.text(0.52, bottom_y_pos, r'$10^{-8}$', ha='center')
#     fig.text(0.705, bottom_y_pos, r'$10^{0}$', ha='center')
#     fig.text(0.895, bottom_y_pos, r'$10^{8}$', ha='center')

#     # figure letters
#     for letter, ax in zip(['a', 'b', 'c', 'd', 'e'], [ax1, ax2, ax3, ax4, ax5]):
#         at = AnchoredText(
#             letter, prop=dict(size=22), frameon=False, loc='upper left')
#         at.patch.set_boxstyle("round,pad=0.,rounding_size=0.3")
#         ax.add_artist(at)

#     plt.subplots_adjust(left=0.08,
#                         bottom=0.099,
#                         right=0.91,
#                         top=0.954,
#                         wspace=0.01,
#                         hspace=0.02)

#     plt.subplot_tool()

#     plt.show()
# #
# plt.savefig('/home/ashalev/Projects/amos-try/multem2mod/multem3article/figures/ready_to_publish/fig1_new.pdf')
# plt.clf(); plt.close()
# ----------------------------------------------------------------------------------------------------------






