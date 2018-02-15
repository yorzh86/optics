from __future__ import division
import numpy as np
import pylab as pl


def plot_Eps4(ax, wl, args):
    ax.plot(wl, args[0], label= args[5], color='r')
    ax.plot(wl, args[1], label= args[6], color='b')
    ax.plot(wl, args[2], label= args[7], color='r', linestyle=':')
    ax.plot(wl, args[3], label= args[8], color='b', linestyle=':')
    ax.plot(wl, wl*0, linestyle = '-.', color='grey', linewidth = 0.3)

    ax.set_ylabel('Epsilon, 'r'$\epsilon$')
    ax.set_xlabel('log10 'r'$\lambda$')
    #ax.set_xlabel('Wavelength, 'r'$\lambda$ [nm]')

    xmajor_ticks = np.arange(2.5, 5.0, 0.5)
    xminor_ticks = np.arange(2.5, 4.55, 0.05)
    #xmajor_ticks = np.arange(500, 24000, 4000)
    #minor_ticks = np.arange(500, 21000, 1000)
    ymajor_ticks = np.arange(-20, 30, 5)
    yminor_ticks = np.arange(-20, 25.5, 0.5)

    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor = True)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor = True)
    ax.legend(loc=args[9], fancybox=True)

def doFigure_Eps4(wl, args):
    fig = pl.figure()
    axEps = fig.add_subplot(111)
    plot_Eps4(axEps, wl, args)
    fig.tight_layout()
    fig.savefig(args[4], dpi=500)
    return

def plot_Eps2(ax, wl, args):
    ax.plot(wl, args[0], label= args[2], color='r')
    ax.plot(wl, args[1], label= args[3], color='b')
    ax.plot(wl, wl*0, linestyle = '-.', color='grey', linewidth = 0.3)

    ax.set_ylabel('Epsilon, 'r'$\epsilon$')
    ax.set_xlabel('log10 'r'$\lambda$')
    #ax.set_xlabel('Wavelength, 'r'$\lambda$ [nm]')

    xmajor_ticks = np.arange(2.5, 5.0, 0.5)
    xminor_ticks = np.arange(2.5, 4.55, 0.05)
    #xmajor_ticks = np.arange(500, 24000, 4000)
    #minor_ticks = np.arange(500, 21000, 1000)
    ymajor_ticks = np.arange(-5, 25, 5)
    yminor_ticks = np.arange(-5, 21, 1)

    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor = True)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor = True)
    ax.legend(loc=args[5], fancybox=True)

def doFigure_Eps2(wl, args):
    fig = pl.figure()
    axEps = fig.add_subplot(111)
    plot_Eps2(axEps, wl, args)
    fig.tight_layout()
    fig.savefig(args[4], dpi=500)
    return

def plot_Rp_Tp(ax, xaxis, TMM, EMTi, EMTst, prop):
    ax.plot(xaxis[:], TMM[:],  label= prop+ ' TMM', color = 'r', linewidth = 0.4)
    ax.plot(xaxis[:], EMTi[:], label= prop+ ' EMTi', color = 'b', linewidth = 0.4)
    ax.plot(xaxis[:], EMTst[:],label= prop+ ' EMTst', color = 'g', linestyle=':')
    ax.yaxis.tick_left()
    ax.set_ylabel(prop +", " + prop[0])

    #ax.set_xlabel('Incidence angle theta, 'r'$\theta$')
    ax.set_xlabel('log10 'r'$\lambda$')
    ax.legend(loc=2, fancybox=True)

    ymajor_ticks = np.arange(0, 1.2,  0.2)
    yminor_ticks = np.arange(0, 1.02, 0.02)

    # 1 angle - many wl:
    xmajor_ticks = np.arange(2.5, 5.0, 0.5)
    xminor_ticks = np.arange(2.5, 4.55, 0.05)
    #xmajor_ticks = np.arange(500, 24000, 4000)
    #xminor_ticks = np.arange(500, 21000, 1000)

    # 1 wl - many angles:
    #xmajor_ticks = np.arange(0, 100, 10)
    #xminor_ticks = np.arange(0, 91, 1)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor = True)

    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor = True)

    #pl.text(0.5, 0.85,"$\lambda$ = 400 nm", horizontalalignment = 'center', verticalalignment = 'center',
    #        transform = ax.transAxes)
    #pl.text(0.7, 0.7,"$\ theta$ = 0", horizontalalignment = 'center', verticalalignment = 'center',
    #        transform = ax.transAxes)

def doFigure_RTA(xaxis, x, y, z, filename, prop):
    fig = pl.figure()
    axR = fig.add_subplot(111)
    plot_Rp_Tp(axR,xaxis, x, y, z, prop)
    fig.tight_layout()
    fig.savefig(filename, dpi=500)
    return

def basic_info(name, N_layers, N_periods):
    print
    print "Material:", name
    print "Overal number of layers:", N_layers
    print "Period number:", N_periods
    return

def writeToFile(fn, text1, text2):
    f = open(fn, 'w')
    for i in range(len(text1)):
        f.write(str(text1[i])+"\t" +str(text2[i])+"\n")
    f.close
    return