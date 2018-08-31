from __future__ import division
import numpy as np
import pylab as pl
import glob, os
import sys
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def plot_Eps4(ax, wl, args):
    ax.semilogx(wl, args[0], label= args[5], color='r')
    ax.semilogx(wl, args[1], label= args[6], color='b')
    ax.semilogx(wl, args[2], label= args[7], color='r', linestyle=':')
    ax.semilogx(wl, args[3], label= args[8], color='b', linestyle=':')
    ax.semilogx(wl, wl*0, linestyle = '-.', color='grey', linewidth = 0.3)

    ax.set_ylabel('Epsilon, 'r'$\epsilon$')
    ax.set_xlabel('Wavelength, 'r'$\lambda$ [nm]')

    #ymajor_ticks = np.arange(-50, 60, 10)
    #yminor_ticks = np.arange(-50, 51, 1.0)

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
    ax.semilogx(wl, args[0], label= args[2], color='black')
    ax.semilogx(wl, args[1], label= args[3], color='black', linestyle='--')
    ax.semilogx(wl, wl*0, linestyle = '-.', color='grey', linewidth = 0.3)

    ax.set_ylabel('Epsilon, 'r'$\epsilon$')
    ax.set_xlabel('Wavelength, 'r'$\lambda$ [nm]')
    #ax.set_xlabel('Wavelength, 'r'$\lambda$ [nm]')
    ax.minorticks_on()

    #xmajor_ticks = np.arange(500, 24000, 4000)
    #minor_ticks = np.arange(500, 21000, 1000)
    #ymajor_ticks = np.arange(-5, 25, 5)
    #yminor_ticks = np.arange(-5, 26, 1)

    #ax.set_yticks(ymajor_ticks)
    #ax.set_yticks(yminor_ticks, minor = True)
    ax.legend(loc=args[5], fancybox=True)

def doFigure_Eps2(wl, args):
    fig = pl.figure()
    axEps = fig.add_subplot(111)
    plot_Eps2(axEps, wl, args)
    fig.tight_layout()
    fig.savefig(args[4], dpi=500)
    return

def plot_Rp_Tp(ax, xaxis, TMM, EMTi, EMTst, prop, l):
    ax.semilogx(xaxis[:], TMM[:],  label= prop+ ' TMM', color = 'r', linewidth = 0.4)
    ax.semilogx(xaxis[:], EMTi[:], label= prop+ ' EMTi', color = 'b', linewidth = 0.4)
    ax.semilogx(xaxis[:], EMTst[:],label= prop+ ' EMTst', color = 'g', linestyle=':')
    ax.yaxis.tick_left()
    ax.set_ylabel(prop +", " + prop[0])

    #ax.set_xlabel('Incidence angle theta, 'r'$\theta$')
    ax.set_xlabel('Wavelength, 'r'$\lambda$ [nm]')
    ax.legend(loc=l, fancybox=True)

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

def doFigure_RTA(xaxis, x, y, z, filename, prop, loc='best'):
    fig = pl.figure()
    axR = fig.add_subplot(111)
    plot_Rp_Tp(axR,xaxis, x, y, z, prop, loc)
    fig.tight_layout()
    fig.savefig(filename, dpi=500)
    return

def basic_info(name, N_layers, N_periods):
    print
    print "Material:", name
    print "Overal number of layers:", N_layers
    print "Period number:", N_periods
    return

def doContourPlot(wl, theta_i, R, filename):
    fig = pl.figure()
    ax = fig.add_subplot(111)
    x,y = np.meshgrid(wl, theta_i)
    #pl.yscale("log")
    ax.set_yscale("log")
    p = ax.contourf(y,x,R)
    pl.xlabel('Angle, 'r'$\Theta$')
    pl.ylabel('Wavelength, 'r'$\lambda$ [nm]')
    pl.title(filename)
    #plt.colobar(p,shrink=0.8, extend='both')
    plt.colorbar(p)
    #ax.set_ybound(0, 20000)
    
    xmajor_ticks = np.arange(0, 105, 15)
    xminor_ticks = np.arange(0, 95, 5)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor = True)
    
    #ymajor_ticks = np.arange(0, 25000, 5000)
    #yminor_ticks = np.arange(0, 95, 5)
    #ax.set_yticks(ymajor_ticks)
    #ax.set_xticks(yminor_ticks, minor = True)
    ax.axes.text(-12,18000, r'$2\cdot 10^{4}$', fontsize=10, transform = ax.transData)
    ax.axes.text(-8,470, '500', fontsize=10, transform = ax.transData)
    fig.savefig('../plots/September/2/'+filename, dpi=500)
    pl.show()
    return
    
    

#================== Write to file ======================#

def writeToFile(fn, title, data):
    f = open(fn, 'w')

    f.write(title+"\n")

#    for k in data:
#        print k

    for j in range(len(data[0])):
        for i in range(len(data)):
            f.write(str(data[i][j]) + '\t')
            if (i == len(data)-1):
                f.write('\n')
        f.close
    return

#================== Read from file =====================#

def readBulk(fn):
    D = pl.loadtxt(fn,skiprows=1)
    wl = D[:,0]
    er = D[:,1]
    ei = D[:,2]
    return wl, er, ei

def readCond(fn):
    D = pl.loadtxt(fn,skiprows=1)
    wl = D[:,0]
    er = D[:,1]
    o_r = D[:,2]
    ei = D[:,3]
    o_i= D[:,4]
    return wl, er, o_r, ei, o_i


def plot_Fig2Bulk(ax, wl, er, ei, loc1, num):
    ax.semilogx(wl, er, label= "Epsilon real", color='r')
    ax.semilogx(wl, ei, label= "Epsilon real", color='r', linestyle=':')
    ax.semilogx(wl, wl*0, linestyle = '-.', color='grey', linewidth = 0.3)

    ax.set_ylabel('Epsilon, 'r'$\epsilon$')
    ax.set_xlabel('Wavelength, 'r'$\lambda$ [nm]')

    ymajor_ticks = np.arange(-10, 30, 5)
    yminor_ticks = np.arange(-10, 26, 1)

    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor = True)
    ax.legend(loc=loc1, fancybox=True, fontsize='xx-small')

    ax.text(180,25, num)

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize('x-small')

def plot_Fig2Cond(ax, wl, er, o_r, ei, o_i, loc1, num):
    ax.semilogx(wl, er, label="Extraordinary, real", color='r')
    ax.semilogx(wl, o_r, label= "Ordinary, real", color='b')
    ax.semilogx(wl, ei, label= "Extraordinary, imaginary", color='r', linestyle=':')
    ax.semilogx(wl, o_i, label= "Ordinary, imaginary", color='b', linestyle=':')
    ax.semilogx(wl, wl*0, linestyle = '-.', color='grey', linewidth = 0.3)

    ax.set_ylabel('Epsilon, 'r'$\epsilon$')
    ax.set_xlabel('Wavelength, 'r'$\lambda$ [nm]')

    ymajor_ticks = np.arange(-280, 120, 40)
    yminor_ticks = np.arange(-280, 90, 10)

    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor = True)
    ax.legend(loc=loc1, fancybox=True, fontsize='xx-small')

    ax.text(180, 75, num)

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize('x-small')

def doFigure2(BiSe_B, BiTe_B, BiSe_C, BiTe_C):

    wl,er,ei = readBulk(BiSe_B)
    wl1,er1,ei1 = readBulk(BiTe_B)

    wl2,er2,o_r2,ei2,o_i2 = readCond(BiSe_C)
    wl3,er3,o_r3,ei3,o_i3 = readCond(BiTe_C)

    fig = pl.figure()

    bulk_BiSe = fig.add_subplot(221)
    plot_Fig2Bulk(bulk_BiSe, wl, er, ei, 2, '(a)')

    bulk_BiTe = fig.add_subplot(223)
    plot_Fig2Bulk(bulk_BiTe, wl, er, ei, 2, '(c)')

    cond_BiSe = fig.add_subplot(222)
    plot_Fig2Cond(cond_BiSe,wl2,er2,o_r2,ei2,o_i2, 3, '(b)')

    cond_BiTe = fig.add_subplot(224)
    plot_Fig2Cond(cond_BiTe,wl3,er3,o_r3,ei3,o_i3, 3, '(d)')


    fig.tight_layout()
    fig.savefig('../plots/April/Fig2/test.png',dpi=1000)
    return

directory = '../plots/April/Fig2/'
BiSe_B = directory + 'BiSe_Bulk.xls'
BiTe_B = directory + 'BiTe_Bulk.xls'

BiSe_C = directory + 'BiSe_Cond.xls'
BiTe_C = directory + 'BiTe_Cond.xls'

#doFigure2(BiSe_B, BiTe_B, BiSe_C, BiTe_C)


#args=[er, ei,"Epsilon real","Epsilon imaginary", title, 3]



