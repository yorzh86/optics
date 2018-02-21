from __future__ import division
import numpy as np
import pylab as pl
import glob, os
import sys


def plot_Eps4(ax, wl, args):
    ax.semilogx(wl, args[0], label= args[5], color='r')
    ax.semilogx(wl, args[1], label= args[6], color='b')
    ax.semilogx(wl, args[2], label= args[7], color='r', linestyle=':')
    ax.semilogx(wl, args[3], label= args[8], color='b', linestyle=':')
    ax.semilogx(wl, wl*0, linestyle = '-.', color='grey', linewidth = 0.3)

    ax.set_ylabel('Epsilon, 'r'$\epsilon$')
    ax.set_xlabel('Wavelength, 'r'$\lambda$ [nm]')

    ymajor_ticks = np.arange(-200, 85, 25)
    yminor_ticks = np.arange(-200, 80, 5)

    #ymajor_ticks = np.arange(-300, 75, 25)
    #yminor_ticks = np.arange(-300, 55, 5)

    #ax.set_xticks(xmajor_ticks)
    #ax.set_xticks(xminor_ticks, minor = True)
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
    ax.semilogx(wl, args[0], label= args[2], color='r')
    ax.semilogx(wl, args[1], label= args[3], color='b')
    ax.semilogx(wl, wl*0, linestyle = '-.', color='grey', linewidth = 0.3)

    ax.set_ylabel('Epsilon, 'r'$\epsilon$')
    ax.set_xlabel('Wavelength, 'r'$\lambda$ [nm]')
    #ax.set_xlabel('Wavelength, 'r'$\lambda$ [nm]')


    #xmajor_ticks = np.arange(500, 24000, 4000)
    #minor_ticks = np.arange(500, 21000, 1000)
    ymajor_ticks = np.arange(-60, 70, 10)
    yminor_ticks = np.arange(-60, 62, 2)

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
    ax.semilogx(xaxis[:], TMM[:],  label= prop+ ' TMM', color = 'r', linewidth = 0.4)
    ax.semilogx(xaxis[:], EMTi[:], label= prop+ ' EMTi', color = 'b', linewidth = 0.4)
    ax.semilogx(xaxis[:], EMTst[:],label= prop+ ' EMTst', color = 'g', linestyle=':')
    ax.yaxis.tick_left()
    ax.set_ylabel(prop +", " + prop[0])

    #ax.set_xlabel('Incidence angle theta, 'r'$\theta$')
    ax.set_xlabel('Wavelength, 'r'$\lambda$ [nm]')
    ax.legend(loc=1, fancybox=True)

    ymajor_ticks = np.arange(0, 1.2,  0.2)
    yminor_ticks = np.arange(0, 1.02, 0.02)


    # 1 wl - many angles:
    #xmajor_ticks = np.arange(0, 100, 10)
    #xminor_ticks = np.arange(0, 91, 1)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor = True)

    #ax.set_xticks(xmajor_ticks)!!!
    #ax.set_xticks(xminor_ticks, minor = True)!!!

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


#======================Messy Penetration Depth======================
def writeToFile2(fn, text1, text2):
    f = open(fn, 'w')
    for i in range(len(text1)):
        f.write(str(text1[i])+"\t" +str(text2[i])+"\n")
    f.close
    return

def writeToFile3(fn, text1, text2, text3, text4):
    f = open(fn, 'w')
    for i in range(len(text1)):
        f.write(str(text1[i])+"\t" + str(text2[i])+"\t"+str(text3[i])+"\t" + str(text4[i])+"\n")
    f.close
    return

def writeToFile4(fn, text1, text2, text3, text4, text5):
    f = open(fn, 'w')
    for i in range(len(text1)):
        f.write(str(text1[i])+"\t" + str(text2[i])+"\t" +str(text3[i])+"\t"+str(text4[i])+"\t"+ str(text5[i]) + "\n")
    f.close
    return

def readFile(fn):
    f = pl.loadtxt(fn, skiprows=0)
    wl = f[:,0]
    Pd = f[:,1]
    return wl, Pd

def plot_PD(ax, wl, PD, fn):
    ax.plot(wl, PD[0], label= fn[0][0:-12] , color='r')
    ax.plot(wl, PD[1], label= fn[1][0:-12], color='b')
    ax.plot(wl, PD[2], label= fn[2][0:-12], color='g')
    #ax.plot(wl, wl*0, linestyle = '-.', color='grey', linewidth = 0.3)

    ax.set_ylabel('Penetration depth [nm] ')
    ax.set_xlabel('log10 'r'$\lambda$')
    #ax.set_xlabel('Wavelength, 'r'$\lambda$ [nm]')

    xmajor_ticks = np.arange(500, 24000, 4000)
    xminor_ticks = np.arange(500, 21000, 1000)
    ymajor_ticks = np.arange(0, 550000, 50000)
    yminor_ticks = np.arange(0, 510000, 10000)

    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor = True)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor = True)
    ax.legend(loc='best', fancybox=True)

def doFigure_PD(wl, PD, fn):
    fig = pl.figure()
    axEps = fig.add_subplot(111)
    plot_PD(axEps, wl, PD, fn)
    fig.tight_layout()
    fig.savefig('PD.png', dpi=500)
    #fig.savefig(fn[0:-4]+'.png', dpi=500)
    return


def make_pd():
    directory = '../plots/updateFeb16/Bi2Te3/diel100/'

    wl = np.zeros((1, 500), dtype=float)
    fn = np.chararray(3, 32)
    pd = np.zeros((3, 500), dtype=float)

    j = 0
    os.chdir(directory)
    for file in glob.glob("*EMTi.txt"):
        a, b = readFile(file)
        wl[0] = a
        pd[j] = b
        fn[j]= file
        j=j+1

    doFigure_PD(wl[0], pd, fn)
    return

#make_pd()