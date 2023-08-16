import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.patches import Rectangle
from matplotlib.offsetbox import AnchoredText
from matplotlib.patches import Arc
# import seaborn as sns
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
    # prop_patches["fc"] = "green"
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





# fig.3 insertion (now it differs from MB)-----
# plt.figure(figsize=(12, 7))
# # teal dashed
# data = np.loadtxt('lmax4_all/akxy=0.001_0.0.txt')
# x, y = prep_xy_data(data)
# plt.plot(x, y, color='cyan', linestyle='--')
# # black solid
# data = np.loadtxt('lmax4_all/akxy=8e-05_0.0.txt')
# x, y = prep_xy_data(data)
# plt.plot(x, y, color='black')
#
# plt.show()


def plot_error(ax, x,y,z, title='', vextr=None, n_orders_div=1):
    # fig, ax = plt.subplots(figsize=(8,6))
    if vextr is not None:
        z = np.clip(a=z, a_min=vextr[0], a_max=vextr[1])
        n_orders = int(np.ceil(np.abs(np.log10(vextr[0]) -
                                      np.log10(vextr[1])))) // n_orders_div #// 2 + 1
        cmap_idv = matplotlib.colors.ListedColormap(sns.color_palette("RdBu_r", n_orders))
        im  = ax.imshow(
            # np.vectorize(mpmath.log10)(z).astype(np.float64),
            z.astype(np.float64),
            origin="lower",extent=[exp_min,exp_max,exp_min,exp_max],aspect='auto',
            norm=LogNorm(vmin=vextr[0], vmax=vextr[1]),
            cmap=cmap_idv)
    else:
        im = ax.imshow(np.vectorize(mpmath.log10)(z).astype(np.float64),
                       origin="bottom",extent=[exp_min,exp_max,exp_min,exp_max],aspect='auto',
                       cmap=cmap_std)

    ax.set_title(title)
    # ax.set_xlabel('real')
    # ax.set_ylabel('imaginary')
    # ax.set_yticklabels([3.731, 3.732, 0, 3.734, 3.735])

    # ax.grid(True, c='black', ls='--')

    return im


def create_arc_and_arrow(ax, x1, y1, len_x, len_y, fc, alpha):
    # # Configure arc
    # # center_x = 2            # x coordinate
    # # center_y = 3.8          # y coordinate
    # radius_1 = 0.0001        # radius 1
    # radius_2 = 500*radius_1           # radius 2 >> for cicle: radius_2 = 2 x radius_1
    # angle = 180             # orientation
    # theta_1 = 0          # arc starts at this angle
    # theta_2 = 300         # arc finishes at this angle
    # arc = Arc([center_x, center_y],
    #           radius_1,
    #           radius_2,
    #           angle = angle,
    #           theta1 = theta_1,
    #           theta2=theta_2,
    #           capstyle = 'round',
    #           linestyle='-',
    #           lw=2,
    #           color = 'black')

    # Add arc
    # ax.add_patch(arc)

    # Add arrow
    # x1 = 3.733335          # x coordinate
    # y1 = 1.0            # y coordinate
    # len_x = 0.00005     # length on the x axis (negative so the arrow points to the left)
    # len_y = 0        # length on the y axis
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
    # ax.annotate(xy=(x1, y1), dxdy=(len_x, len_y),
    #             arrowprops={'arrowstyle': '->', 'lw': 2, 'color': fc},
    #             va='center')

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#--------- customizing font --------
# plt.rcParams['text.usetex'] = True
# print(plt.style.available)
# print(plt.rcParams.keys())
plt.rcParams.update({'font.size': 28, 'font.serif':"Times New Roman"})
# rc('font',**{'family':'serif','serif':['Times']})
# rc('text', usetex=True)


figures_to_plot = ['fig5']


# ------------- Fig. 2 ---------------------------------------------------------------------------------------------
if 'fig2' in figures_to_plot:
    fig = plt.figure(figsize=(20, 7*3))
    gs = GridSpec(3, 2, figure=fig)
    # # fig 2 a
    ax1 = fig.add_subplot(gs[0, :])
    core_dir = 'data/fig2a'
    # red
    data = np.loadtxt(core_dir+'/akxy=0.01_0.0.txt')
    x, y = prep_xy_data(data)
    ax1.plot(x, y, color='red', label=r'$0.01$'+'\t'+ r'$0$', lw=3)
    #silver
    data = np.loadtxt(core_dir+'/akxy=0.0_0.05.txt')
    x, y = prep_xy_data(data)
    ax1.plot(x, y, color='silver', label=r'$0$'+'\t'+ '  0.05', lw=3)
    #blue dashed
    data = np.loadtxt(core_dir+'/akxy=0.05_0.0.txt')
    x, y = prep_xy_data(data)
    ax1.plot(x, y, 'b--', dashes=(15, 15), label=r'$0.05$'+'\t'+ r'$0$', lw=3)
    #magenta dot-dashed
    data = np.loadtxt(core_dir+'/akxy=0.1_0.0.txt')
    x, y = prep_xy_data(data)
    ax1.plot(x, y, color='magenta', linestyle= (0, (3, 5, 1, 5)), label=r'$0.1$'+'\t'+' 0', lw=3)
    #------- customizing
    ax1.set_yticks([0.0, 1.0], minor=False)
    ax1.set_xlabel(r'$ak_0$')
    ax1.set_ylabel(r'$\mathbf{T}$', rotation=0, fontsize=30)
    # ax1.set_xlim(3.7295, 3.74)
    ax1.set_xlim(3.7285, 3.74)
    ax1.set_ylim(-0.1, 1.1)
    ax1.legend(frameon=False, loc=2, bbox_to_anchor=(0.0,0.86) ,title='\t'+r'$ak_x$'+'\t'+r'$ak_y$')
    
    # fig.3 b
    core_dir = 'data/fig2b'
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1], sharey=ax2)  
    for lmax in [4, 7, 10, 13]:
        data = np.loadtxt(core_dir+f'/lmax={lmax}.txt')
        x, y = prep_xy_data(data)
        if lmax == 4:
            ax2.plot(x, y, label=str(lmax), color='red', lw=3)
            ax3.plot(x, y, label=str(lmax), color='red', lw=3)
        else:
            ax2.plot(x, y, label=str(lmax), lw=3)
            ax3.plot(x, y, label=str(lmax), lw=3)
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
    # fig.3 c
    ax4 = fig.add_subplot(gs[2, :])
    data_dir = 'data/fig2c/'
    data = np.loadtxt(data_dir+'T.txt')
    x, T = prep_xy_data(data)
    x, T = x, T
    error_wo_lapack = prep_xy_data(np.loadtxt(data_dir+'error_wo_lapack.txt'))[1]
    error_with_lapack = prep_xy_data(np.loadtxt(data_dir+'error_with_lapack.txt'))[1]
    error_with_lapack_and_fad = prep_xy_data(np.loadtxt(data_dir+'error_with_lapack.txt'))[1]
    ax44 = ax4.twinx()
    l1, = ax44.plot(x[::2], np.array(error_with_lapack[::2])*1e4, color='C2', lw=3)
    l2, = ax44.plot(x, error_wo_lapack, color='C0', lw=3)
    ax4.set_ylabel('error')
    l3, = ax4.plot(x, T, 'red', lw=6, alpha=0.3)
    # create_arc_and_arrow(ax11, 3.733282, 1.002)
    ax4.set_ylabel(r'$\mathbf{T}$', rotation=0, fontsize=30)
    plt.legend([l1, l2], ['LAPACK', 'W/O LAPACK'], frameon=False)
    plt.annotate(r"x$10^4$",
                xy=(3.731546, 0.6), xycoords='data',    # 0.6 1 option; 0.34 2 option
                xytext=(-19, 30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->", color='black')
                )
    # customizing
    ax4.set_yticks([0.0, 1.0], minor=False)
    # ax4.tick_params(axis='y', colors='red')
    ax4.set_ylim([-0.25, 1.1])
    ax4.set_xlim([3.731, 3.735])
    ax4.set_xlabel(r'$ak_0$')
    ax4.set_xticks([3.731, 3.732, 3.733, 3.734, 3.735], minor=False)
    ax44.set_ylim([-0.009, 0.042])
    ax44.set_yticks([0, 0.04], minor=False)  # -0.0235 1 option; -0.005 2 option
    ax44.set_yticklabels(['0', '0.04'])
    ax44.set_ylabel('error', labelpad=-30)
    arrow_len = 0.00012
    create_arc_and_arrow(ax4, 3.7321, 0.9, -arrow_len, 0.0, 'pink', 1.0)
    create_arc_and_arrow(ax4, 3.7346, 0.2, arrow_len, 0.0, 'C2', 1.0)
    create_arc_and_arrow(ax4, 3.7333, 0.24, arrow_len, 0.0, 'C0', 1.0)
    # annotations on the whole figure
    fig.text(0.49, 0.34, r'$ak_0$', ha='center')
    fig.text(0.17, 0.28, r'$l_{max}=10$', ha='center')
    fig.text(0.17, 0.93, r'$l_{max}=4$', ha='center')
    fig.text(0.02, 0.97, 'a)', ha='center')
    fig.text(0.02, 0.65, 'b)', ha='center')
    fig.text(0.02, 0.33, 'c)', ha='center')
    fig.text(0.8, 0.14, r'$x10^4$', ha='center', color='C2', fontsize=26)
    # ---------- zooming
    zoom_effect01(ax1, ax3, 3.7349, 3.7355)
    zoom_effect01(ax2, ax4, 3.7321, 3.7322)
    # plt.subplot_tool()
    plt.subplots_adjust(left=0.07,
                        bottom=0.07,
                        right=0.915,
                        top=0.963,
                        wspace=0.033,
                        hspace=0.297)
    
    
    # delete this
    # ax5 = fig.add_subplot(gs[3, :])
    # ax55 = ax5.twinx()

    # l1, = ax55.plot(x, error_wo_lapack, color='C2')
    # l4, = ax55.plot(x, np.array(error_with_lapack_and_fad)*1e4, color='C3')
    # ax5.set_ylim([-0.25, 1.1])
    # ax5.set_xlim([3.731, 3.735])
    # l3, = ax5.plot(x, T, 'red', lw=4, alpha=0.3)
    # plt.legend([l1, l4], ['W/O LAPACK', 'LAPACK and FADDEEVA'], frameon=False)
    # ax55.set_ylim([-0.009, 0.042])
    # ax55.set_yticks([0, 0.04], minor=False)  # -0.0235 1 option; -0.005 2 option
    # ax55.set_yticklabels(['0', '0.04'])
    # ax55.set_ylabel('error', labelpad=-30)


    # customizing ticks 
    for ax in [ax1, ax2, ax3, ax4, ax44]:
        ax.tick_params('both', length=10, width=2, which='major')


    # plt.show()
    plt.savefig('fig2_new.pdf')
    plt.clf(); plt.close()


if 'fig3' in figures_to_plot:
    plt.rcParams.update({'font.size': 50, 'font.serif':"Times New Roman"})

    fig = plt.figure(figsize=(40, 20))
    # --- reading data ------
    data_dir = 'data/fig3/'
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
        line = line.replace("(","")
        line = line.replace(")","")
        re, im = line.split(",")
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
    ax1.scatter(np.log10(pos_args.real)[::30], np.log10(pos_args.imag)[::30], color='black', s=15)
    ax1.set_yticks([-8, 0, 8])
    ax1.set_yticklabels([])
    # plotting argument regions
    data_dir = 'data/fig3/'
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

    # ax2.scatter(np.log10(pos_args.real)[::30], np.log10(pos_args.imag)[::30], color='black', s=1)

    # ax3
    ax3 = fig.add_subplot(2,2,3)
    im = plot_error(ax3, x2d, y2d, multem3_re_err, vextr=[err_lim, 1e-6])
    ax3.set_yticks([-8, 0, 8])
    ax3.set_yticklabels([])
    ax3.set_xticks([-8, 0, 8])
    ax3.set_xticklabels([])
    # ax3.scatter(np.log10(pos_args.real)[::30], np.log10(pos_args.imag)[::30], color='black', s=1)

    # ax4
    ax4 = fig.add_subplot(2,2,4, sharey=ax3)
    im = plot_error(ax4, x2d, y2d, multem3_im_err, vextr=[err_lim, 1e-6])
    ax4.get_yaxis().set_visible(False)
    ax4.set_xticks([-8, 0, 8])
    ax4.set_xticklabels([])
    # ax4.scatter(np.log10(pos_args.real)[::30], np.log10(pos_args.imag)[::30], color='black', s=1)


    # ax5 (arg values from MB Dirac cone)
    ax5 = fig.add_axes([0.34, 0.6, 0.15, 0.35])
    # plotting argument regions
    ax5.contour(regimes2d, zorder=1, levels=2, linewidths=7, colors='yellow', linestyles='dashed', extent=[exp_min,exp_max,exp_min,exp_max])
    im = plot_error(ax5, x2d, y2d, multem2_re_err, vextr=[err_lim, 1e-6])
    ax5.scatter(np.log10(pos_args.real)[::5], np.log10(pos_args.imag)[::5], color='black', s=30)
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
    # fig.text(left_x_pos, 0.25, 'Multem 3', ha='center', rotation=90)
    # fig.text(left_x_pos, 0.7, 'Multem 2', ha='center', rotation=90)
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
            letter, prop=dict(size=fontsize), frameon=False, loc='upper left')
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

    # plt.subplot_tool()
    #


    # plt.show()
    plt.savefig('fig3_new.pdf')
    plt.clf(); plt.close()


import calculating_lib as cl

if 'fig4' in figures_to_plot:

    plt.rcParams.update({'font.size': 20, 'font.serif':"Times New Roman"})


    dir = f'data/fig4/'
    fig = plt.figure(figsize=(16, 9))
    i = 1
    RMAX = [7, 8, 9]
    color = {7: 'red', 8: 'black', 9: 'green'}
    lw = {7: 4, 8: 2, 9: 8}
    alpha = {7: 1, 8: 1, 9: 0.3}
    lines = {7: None, 8: None, 9: None}

    for rmax in RMAX:
        R1, R2, R3, x = multi_loadtxt(dir, (f'650_{rmax}.txt', f'750_{rmax}.txt', f'900_{rmax}.txt', 'sintheta.txt'))
        plt.plot(x, R1, color=color[rmax], lw=lw[rmax], alpha=alpha[rmax])
        plt.plot(x, R2+1, color=color[rmax], lw=lw[rmax], alpha=alpha[rmax])
        lines[rmax], = plt.plot(x, R3+2, color=color[rmax], lw=lw[rmax], alpha=alpha[rmax])
        plt.xlabel(r'$sin\theta$')
        plt.ylabel(r'$R_{spec}$')
    print(lines)
    plt.legend(handles=[lines[7], lines[8], lines[9]], labels=['RMAX=7', 'RMAX=8', 'RMAX=9'])
    plt.tick_params('both', length=7, width=2, which='major')    

    fig.text(0.5, 0.23, r'$\lambda=650 nm$', ha='center')
    fig.text(0.5, 0.48, r'$\lambda=750 nm$', ha='center')
    fig.text(0.5, 0.8, r'$\lambda=900 nm$', ha='center')
    

    plt.savefig('fig4_PRB_2017.pdf')
    plt.clf(); plt.close()
    # plt.show()




if 'fig6' in figures_to_plot:
    # ---- ERROR BTW RMAX 7 AND RMAX 29 --------
    # --- JUST SPECTRA -----
    LMAX = [4]
    RMAX = [7]
    AK1 = [0.01, 8e-5, 4e-5]
    X = {}; Y = {}
    step = 1
    #reading data to plot
    for version in ['3_cerf', '3']:
        main_dir = dir=f'data/fig4/{version}/'
        for lmax in LMAX:
            for rmax in RMAX:
                for ak1 in AK1:
                    ak1 = round(ak1/2/np.pi, 8)
                    path = main_dir+f'lmax={lmax}_rmax={rmax}_ak1={ak1}.txt'
                    x, y = cl.read_1D_data_real(path)
                    X[f'{version}_{rmax}_{ak1}'] = x; Y[f'{version}_{rmax}_{ak1}'] = y
    
    fig = plt.figure(figsize=(40, 20))
    plt.rcParams.update({'font.size': 45, 'font.serif':"Times New Roman"})
    error_step = 5
    axs = []
    for i, ak1 in enumerate(AK1):
        ak1 = round(ak1/2/np.pi, 8) 
        ax = fig.add_subplot(2, 3, i+1)
        axs.append(ax)
        key = f'{RMAX[0]}_{ak1}'
        plt.gca().plot(X['3_'+key], Y['3_'+key], 'black', lw=5, label='SciPy')
        plt.gca().plot(X['3_cerf_'+key], Y['3_cerf_'+key], 'green', lw=15, alpha=0.3, label='CERF') 
        plt.gca().set_title(r'$ak_x$='+f'{AK1[i]}')
        plt.gca().set_ylabel('T')
        plt.gca().set_xlabel(r'$ak_0$')
        plt.gca().set_yscale('log')
        plt.gca().legend(frameon=False)
        plt.gca().tick_params('both', length=15, width=4, which='both')
        plt.gca().set_ylabel(r'$\mathbf{T}$', rotation=0)
        ax = fig.add_subplot(2, 3, i+4)
        axs.append(ax)
        error = np.abs(Y['3_'+key]-Y['3_cerf_'+key])
        plt.gca().scatter(X['3_'+key][::error_step], error[::error_step], color='black', s=30)
        # plt.gca().set_ylabel('abs(T_Faddeeva - T_cerf)')
        plt.gca().set_ylabel('error')
        plt.gca().set_xlabel(r'$ak_0$')
        plt.gca().set_yscale('log')
        plt.gca().tick_params('both', length=15, width=4, which='both')
    print(len(axs))
    print(axs[4])
    print(len(axs))

    fig.text(0.08, 0.97, r'a)', ha='center', fontsize=50)
    fig.text(0.4, 0.97, r'b)', ha='center', fontsize=50)
    fig.text(0.71, 0.97, r'c)', ha='center', fontsize=50)
    fig.text(0.08, 0.48, r'd)', ha='center', fontsize=50)
    fig.text(0.4, 0.48, r'e)', ha='center', fontsize=50)
    fig.text(0.71, 0.48, r'f)', ha='center', fontsize=50)



    print(axs[5])
    #setting ticks and labels
    # axs[2].set_xticks([3.73556, 3.73558, 3.73560, 3.73562])
    # axs[2].set_xticklabels([3.73556, 3.73558, 3.73560, 3.73562])
    axs[2].set_xticks([3.73556, 3.73560])
    axs[2].set_xticklabels([3.73556, 3.73560])


    axs[4].set_xticks([3.73558, 3.7355825, 3.735585])
    axs[4].set_xticklabels([3.73558, 3.7355825, 3.735585])
    axs[4].set_ylim([1e-5, 2])

    # axs[3].set_xticks([3.73556, 3.73558, 3.73560, 3.73562])
    # axs[3].set_xticklabels([3.73556, 3.73558, 3.73560, 3.73562])
    axs[3].set_xticks([3.73556, 3.73560])
    axs[3].set_xticklabels([3.73556, 3.73560])

    axs[5].set_xticks([3.73558, 3.7355825, 3.735585])
    axs[5].set_xticklabels([3.73558, 3.7355825, 3.735585])
    plt.subplots_adjust(left=0.08,
                        bottom=0.099,
                        right=0.95,
                        top=0.95,
                        wspace=0.25,
                        hspace=0.35)
    plt.savefig(f'fig6.pdf')
    plt.clf(); plt.close()



if 'fig5' in figures_to_plot:
    lmax = 4
    # RMAX = ['7', '9', '11', '13', '15', '17']
    RMAX = ['8', '22', '26', '27']
    main_dir = 'data/fig5/mode=3/version=with_lapack/zoomed/d=1.0/'
    fig = plt.figure(figsize=(15, 15))
    i = 0
    T = {}
    i += 1
    X = {}; Y = {}
    fig.add_subplot(1, 1, i)
    for rmax in RMAX:
        if rmax == '8':
            lw = 7
            alpha = 0.5
        elif rmax == '26':
            lw = 7
            alpha = 0.5
        elif rmax == '27':
            lw = 3
            alpha = 1
        else:
            lw = 1
            alpha = 1
        path = main_dir+f'lmax={lmax}_rmax={rmax}.txt'
        x, y = cl.read_1D_data_real(path)
        X[rmax] = x; Y[rmax] = y
        plt.gca().plot(X[rmax], Y[rmax], label=f'RMAX={rmax}', lw=lw, alpha=alpha)
    # plt.gca().set_yscale('log') 
    plt.gca().legend()
    plt.gca().set_ylabel('T')
    plt.gca().set_xlabel(r'$ak_0$')
    plt.savefig(f'fig5_3.pdf')
    plt.clf(); plt.close()

    # error = np.abs(Y[f'3_{rmax}']-Y[f'3_cerf_{rmax}'])
    # plt.gca().scatter(X[f'3_{rmax}'][::step], y[::step], label=f'max_error {np.max(y)}', s=10)

                    # if version == '3_cerf':
                    #     plt.gca().plot(X[f'{version}_{rmax}'], Y[f'{version}_{rmax}'], label=f'{version}', lw=5, alpha=0.7)
                    # else:
                    #     plt.gca().plot(X[f'{version}_{rmax}'], Y[f'{version}_{rmax}'], label=f'{version}')
                        
    # plt.gca().set_title(f'LMAX={LMAX[0]} RMAX={RMAX}')
    # plt.gca().set_ylabel('abs(T_Faddeeva - T_cerf)')
    # plt.gca().set_xlabel(r'$ak_0$')
    # for rmax in RMAX:
        # plt.gca().plot(X[f'3_{rmax}']*2*np.pi, np.abs(Y[f'3_{rmax}']-Y[f'3_cerf_{rmax}'])*factor, label=f'{rmax}')
        # plt.gca().plot(X[f'3_{rmax}'], Y[f'3_{rmax}']*factor, label=f'{rmax}')
    # plt.gca().set_yscale('log')
    # plt.gca().legend()
    # fig.text(0.25, 0.75, 'RMAX', ha='center')


    # fig.add_subplot(3, 1, 2)
    # plt.gca().plot(X['3_1.0']*2*np.pi, np.abs(Y['3_1.0']-Y['3_wo_faddeeva_1.0'])*factor, label=f'd=1.0')
    # plt.gca().set_yscale('log')
    # plt.gca().legend()

    # fig.add_subplot(3, 1, 3)
    # plt.gca().plot(X['3_1.8']*2*np.pi, np.abs(Y['3_1.8']-Y['3_wo_faddeeva_1.8'])*factor, label=f'd=1.8')
    # plt.gca().set_yscale('log')

    # plt.gca().plot(X['3'], np.abs(Y['3']-Y['3_wo_faddeeva'])*factor, label=f'v={version}_d={d}')
    # plt.gca().set_yscale('log')
    # fig.add_subplot(3, 1, 2)
    
    # fig.add_subplot(3, 1, 3)
    # plt.gca().plot(X['3'], np.abs(Y['with_lapack']-Y['3_wo_faddeeva'])*factor, label=f'v={version}_d={d}')
    # plt.gca().set_yscale('log')

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
    #                     top=0.93,
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
#         lines = f.readlines()
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
#     fig.text(0.73, 0.97, r'$\Im(w(z))$', ha='center')

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

# ------------- Fig. 2 errors lapack ---------------------------------------------------------------------------------------------
# fig = plt.figure(figsize=(16, 9))
# data_dir = 'figures/Maksimov_Bulgakov/data/lapack/'
# data = np.loadtxt(data_dir+'T.txt')
# x, T = prep_xy_data(data)
# x, T = x[:400], T[:400]
# error_wo_lapack = prep_xy_data(np.loadtxt(data_dir+'error_wo_lapack.txt'))[1][:400]
# error_with_lapack = prep_xy_data(np.loadtxt(data_dir+'error_with_lapack.txt'))[1][:400]
# ax1 = plt.subplot()
# l2, = ax1.plot(x, error_wo_lapack, lw=3)
# l1, = ax1.plot(x, np.array(error_with_lapack)*1e4)
# ax1.set_ylabel('error')
# ax2 = ax1.twinx()
# l3, = ax2.plot(x, T, 'r', lw=4, alpha=0.3, zorder=10)
# create_arc_and_arrow(ax2, 3.733282, 1.002)
# ax2.set_ylabel('T')
# plt.legend([l1, l2], ['LAPACK', 'W/O LAPACK'], frameon=False)
# plt.annotate(r"x$10^4$",
#             xy=(3.731546, 0.6), xycoords='data',    # 0.6 1 option; 0.34 2 option
#             xytext=(-19, 30), textcoords='offset points',
#             arrowprops=dict(arrowstyle="->", color='black')
#             )
# ax1.set_xlabel(r'$ak_0$')
# ax1.set_xticks([3.731, 3.732, 3.733, 3.734, 3.735], minor=False)
# ax1.set_yticks([-0.0235, 0.0235], minor=False)  # -0.0235 1 option; -0.005 2 option
# ax2.set_yticks([0, 1.], minor=False)
#
# plt.subplots_adjust(left=0.11,
#                     bottom=0.09,
#                     right=0.95,
#                     top=0.95,
#                     wspace=0.06,
#                     hspace=0.07)
#
# plt.subplot_tool()

# plt.show()

# plt.savefig('/home/ashalev/Projects/amos-try/multem2mod/multem3article/figures/ready_to_publish/fig2-1option.pdf')
# plt.clf(); plt.close()

# ----------------------------------------------------------------------------------------------------------



# ------------- Fig. 2 ---------------------------------------------------------------------------------------------
# fig = plt.figure(figsize=(10, 10))
# fig = plt.figure(figsize=(20, 7*3))
# gs = GridSpec(3, 2, figure=fig)
#
# # fig 3 a
# ax1 = fig.add_subplot(gs[0, :])
# core_dir = 'figures/Maksimov_Bulgakov/data/lmax4_all'
# #red
# data = np.loadtxt(core_dir+'/akxy=0.01_0.0.txt')
# x, y = prep_xy_data(data)
# ax1.plot(x, y, color='red', label=r'$0.01$'+'\t'+ r'$0$')
# #silver
# data = np.loadtxt(core_dir+'/akxy=0.0_0.049999999999999996.txt')
# x, y = prep_xy_data(data)
# ax1.plot(x, y, color='silver', lw=2, label=r'$0$'+'\t'+ '  0.05')
# #blue dashed
# data = np.loadtxt(core_dir+'/akxy=0.049999999999999996_0.0.txt')
# x, y = prep_xy_data(data)
# ax1.plot(x, y, 'b--', dashes=(15, 15), label=r'$0.05$'+'\t'+ r'$0$')
# #magenta dot-dashed
# data = np.loadtxt(core_dir+'/akxy=0.09999999999999999_0.0.txt')
# x, y = prep_xy_data(data)
# ax1.plot(x, y, color='magenta', linestyle= (0, (3, 5, 1, 5)), label=r'$0.1$'+'\t'+' 0')
# #------- customizing
# ax1.set_yticks([0.0, 1.0], minor=False)
# ax1.set_xlabel(r'$ak_0$')
# ax1.set_ylabel('$T$')
# ax1.set_xlim(3.7295, 3.74)
# ax1.set_ylim(-0.1, 1.1)
# ax1.legend(frameon=False, loc=2, bbox_to_anchor=(0.0,0.86) ,title='\t'+r'$ak_x$'+'\t'+r'$ak_y$')
#
# # fig.3 b
# core_dir = 'figures/Maksimov_Bulgakov/data/akx0.01'
# ax2 = fig.add_subplot(gs[1, 0])
# ax3 = fig.add_subplot(gs[1, 1], sharey=ax2)
# # --- red spectrum ----
# lmax = 4
# data = np.loadtxt(core_dir+f'/lmax={lmax}.txt')
# x, y = prep_xy_data(data)
# ax2.plot(x, y, color='red', label=str(lmax))
# ax3.plot(x, y, color='red', label=str(lmax))
# for lmax in [7, 10, 13]:
#     data = np.loadtxt(core_dir+f'/lmax={lmax}.txt')
#     x, y = prep_xy_data(data)
#     ax2.plot(x, y, label=str(lmax))
#     ax3.plot(x, y, label=str(lmax))
# ax2.set_ylim(-0.1, 1.1)
# ax2.set_xlim(3.7317, 3.7325)
# ax3.set_xlim(3.7348, 3.7356)
# ax2.legend(frameon=False, loc=2, title='\t'+r'$l_{max}$')
# # hiding the spines between ax2 and ax3
# ax2.spines['right'].set_visible(False)
# ax3.spines['left'].set_visible(False)
# ax3.yaxis.tick_right()
# ax3.tick_params(labelright=False)  # don't put tick labels at the top
# # drawing diagonal lines
# d = .015  # how big to make the diagonal lines in axes coordinates
# kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
# ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
# ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)
# kwargs.update(transform=ax3.transAxes)
# ax3.plot((-d, +d), (1 - d, 1 + d), **kwargs)
# ax3.plot((-d, +d), (-d, +d), **kwargs)
#
# #------- customizing
# ax2.set_ylabel('$T$')
# ax2.set_yticks([0.0, 1.0], minor=False)
# ax2.set_xticks([3.7318, 3.7321, 3.7324], minor=False)
# ax3.set_xticks([3.735, 3.7353, 3.7356], minor=False)
#
#
#
#
#
# # fig.3 c
# ax4 = fig.add_subplot(gs[2, :])
# data_dir = 'figures/Maksimov_Bulgakov/data/lapack/'
# data = np.loadtxt(data_dir+'T.txt')
# x, T = prep_xy_data(data)
# x, T = x[:400], T[:400]
# error_wo_lapack = prep_xy_data(np.loadtxt(data_dir+'error_wo_lapack.txt'))[1][:400]
# error_with_lapack = prep_xy_data(np.loadtxt(data_dir+'error_with_lapack.txt'))[1][:400]
# ax44 = ax4.twinx()
# # l1, = ax44.plot(x, np.array(error_with_lapack)*1e4, color='C5')
# l1, = ax44.plot(x, np.array(error_with_lapack)*1e4, color='C2')
#
# l2, = ax44.plot(x, error_wo_lapack, color='C0', lw=2)
# ax4.set_ylabel('error')
# l3, = ax4.plot(x, T, 'red', lw=4, alpha=0.3)
# # create_arc_and_arrow(ax11, 3.733282, 1.002)
# ax4.set_ylabel('T')
# plt.legend([l1, l2], ['LAPACK', 'W/O LAPACK'], frameon=False)
# plt.annotate(r"x$10^4$",
#             xy=(3.731546, 0.6), xycoords='data',    # 0.6 1 option; 0.34 2 option
#             xytext=(-19, 30), textcoords='offset points',
#             arrowprops=dict(arrowstyle="->", color='black')
#             )
#
#
#
#
# # customizing
# ax4.set_yticks([0.0, 1.0], minor=False)
# #ax4.tick_params(axis='y', colors='red')
# ax4.set_ylim([-0.25, 1.1])
# ax4.set_xlim([3.731, 3.735])
# ax4.set_xlabel(r'$ak_0$')
# ax4.set_xticks([3.731, 3.732, 3.733, 3.734, 3.735], minor=False)
# ax44.set_ylim([-0.009, 0.042])
# ax44.set_yticks([0, 0.04], minor=False)  # -0.0235 1 option; -0.005 2 option
# ax44.set_yticklabels(['0', '0.04'])
# ax44.set_ylabel('error', labelpad=-30)
#
#
# arrow_len = 0.00012
# create_arc_and_arrow(ax4, 3.7321, 0.9, -arrow_len, 0.0, 'pink', 1.0)
# create_arc_and_arrow(ax4, 3.7346, 0.2, arrow_len, 0.0, 'C2', 1.0)
# create_arc_and_arrow(ax4, 3.7333, 0.24, arrow_len, 0.0, 'C0', 1.0)
#
#
#
# # annotations on the whole figure
# fig.text(0.49, 0.34, r'$ak_0$', ha='center')
# fig.text(0.17, 0.28, r'$l_{max}=10$', ha='center')
# fig.text(0.17, 0.93, r'$l_{max}=4$', ha='center')
# fig.text(0.02, 0.97, 'a)', ha='center')
# fig.text(0.02, 0.65, 'b)', ha='center')
# fig.text(0.02, 0.33, 'c)', ha='center')
# fig.text(0.8, 0.14, r'$x10^4$', ha='center', color='C2', fontsize=26)
#
#
# # ---------- zooming
# zoom_effect01(ax1, ax3, 3.7349, 3.7355)
# zoom_effect01(ax2, ax4, 3.7321, 3.7322)
#
# # plt.subplot_tool()
# plt.subplots_adjust(left=0.07,
#                     bottom=0.07,
#                     right=0.915,
#                     top=0.963,
#                     wspace=0.033,
#                     hspace=0.297)
#
#
# plt.show()

# plt.savefig('/home/ashalev/Projects/amos-try/multem2mod/multem3article/figures/ready_to_publish/fig2.pdf')
# plt.clf(); plt.close()

# ----------------------------------------------------------------------------------------------------------
# ------ rmax issues
# import matplotlib.cm as cm

# def show_2D_map(x, y, z, zmin, zmax, logscale=False):
#     if logscale:
#         im = plt.imshow(z.T, extent = (np.min(x), np.max(x), np.min(y), np.max(y)), cmap=cm.hot, norm=LogNorm(vmin=zmin, vmax=1), aspect='auto', interpolation = 'none', origin='lower')
#     else:
#         im = plt.imshow(z.T, extent = (np.min(x), np.max(x), np.min(y), np.max(y)), cmap=cm.hot, vmin=zmin, vmax=zmax, aspect='auto', interpolation = 'none', origin='lower')

#     cb = plt.colorbar(im)
#     # plt.show()


# def multi_loadtxt(dir, filelist):
#     output = ()
#     for fname in filelist:
#         out = np.loadtxt(dir+"/"+fname)
#         output += (out,)
#     return output
# # r = f(kx) maps with several diffraction orders
# plt.rcParams.update({'font.size': 20, 'font.serif':"Times New Roman"})


# LMAX = [13]
# RMAX = [7, 29]
# main_dir = 'rmax_issue/'
# R = {}
# for lmax in LMAX:
#     fig = plt.figure(figsize=(30, 10))
#     i = 1
#     for rmax in RMAX:
#         fig.add_subplot(1, 3, i)
#         i += 1
#         dir = main_dir+f'lmax={lmax}_rmax={rmax}'
#         kx, omega, R[f'{rmax}'], T = multi_loadtxt(dir, ('kx.txt', 'omega.txt', 'R.txt', 'T.txt'))
#         show_2D_map(kx, omega, R[f'{rmax}'], 1e-5, 1, logscale=True)
#         plt.gca().set_title(f'lmax={lmax}_rmax={rmax}')
#         plt.gca().set_ylabel(r'$\omega a/{2\pi c}$')
#         plt.gca().set_xlabel(r'${k_x d/\pi}$')
#     fig.add_subplot(2, 3, i)
#     i += 1
#     show_2D_map(kx, omega, np.abs(R['30']-R['16']), 0, 1, logscale=False)
#     plt.gca().set_title('np.abs(R_rmax30 - R_rmax16)')
#     plt.gca().set_ylabel(r'$\omega a/{2\pi c}$')
#     plt.gca().set_xlabel(r'${k_x d/\pi}$')
#
#     fig.add_subplot(2, 3, i)
#     i += 1
#     show_2D_map(kx, omega, np.abs(R['30']-R['7']), 0, 1, logscale=False)
#     plt.gca().set_title('np.abs(R_rmax30 - R_rmax7)')
#     plt.gca().set_ylabel(r'$\omega a/{2\pi c}$')
#     plt.gca().set_xlabel(r'${k_x d/\pi}$')

    # plt.savefig(f'lmax={lmax}.pdf')
    # plt.clf(); plt.close()









