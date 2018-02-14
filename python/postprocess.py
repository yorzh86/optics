import numpy as np
import pylab as pl


def plot_Eps4(ax, wl, a, b, c, d):
    ax.plot(wl, a, label= "a", color='r')
    ax.plot(wl, b, label= "b", color='b')
    ax.plot(wl, c, label= 'c', color='r', linestyle=':')
    ax.plot(wl, d, label= 'd', color='b', linestyle=':')
    ax.plot(wl, wl*0, linestyle = '-.', color='grey', linewidth = 0.3)

    ax.set_ylabel('Epsilon, 'r'$\epsilon$')
    ax.set_xlabel('Wavelength, 'r'$\lambda$ [nm]')

    xmajor_ticks = np.arange(500, 4000, 500)
    xminor_ticks = np.arange(500, 3550, 50)
    #xmajor_ticks = np.arange(500, 24000, 4000)
    #minor_ticks = np.arange(500, 21000, 1000)
    ymajor_ticks = np.arange(-5, 25, 5)
    yminor_ticks = np.arange(-5, 25, 0.5)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor = True)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor = True)
    ax.legend(loc=1, fancybox=True)

def doFigure_Eps4(wl, a, b, c, d, filename):
    fig = pl.figure()
    axEps = fig.add_subplot(111)
    plot_Eps4(axEps, wl, a, b, c, d)
    fig.tight_layout()
    fig.savefig(filename, dpi=500)
    return


def plot_Rp_Tp(ax, xaxis, TMM, EMTi, EMTst):
    ax.plot(xaxis[:], TMM[:], label= 'Transmittance TMM', color = 'r', linewidth = 0.4)
    ax.plot(xaxis[:], EMTi[:], label= 'Transmittance EMTi', color = 'b', linewidth = 0.4)
    ax.plot(xaxis[:], EMTst[:], label= 'Transmittance EMTst', color = 'g', linestyle=':')#linewidth = 0.4)
    ax.yaxis.tick_left()
    ax.set_ylabel('Transmittance, $T$')#, color = 'r')

    #ax.set_xlabel('Incidence angle theta, 'r'$\theta$')
    ax.set_xlabel('Wavelength, 'r'$\lambda$')
    ax.legend(loc=2, fancybox=True)

    ymajor_ticks = np.arange(0, 1.01, 0.2)
    yminor_ticks = np.arange(0, 1.02, 0.02)

    # 1 angle - many wl:
    xmajor_ticks = np.arange(500, 4000, 500)
    xminor_ticks = np.arange(500, 3550, 50)
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
    pl.text(0.2, 0.7,"$\ theta$ = 0", horizontalalignment = 'center', verticalalignment = 'center',
            transform = ax.transAxes)

def doFigure_RTA(xaxis, T, Ti, Tst, filename):
    fig = pl.figure()
    axR = fig.add_subplot(111)
    plot_Rp_Tp(axR,xaxis, T, Ti, Tst)
    fig.tight_layout()
    fig.savefig(filename, dpi=500)
    return
