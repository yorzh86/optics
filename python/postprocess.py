from __future__ import division
import glob, os
import sys
import math
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import matplotlib


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

def doContourPlot(wl, theta_i, R, foldername, filename, prop, style='gray'):
    fig = pl.figure()
    ax = fig.add_subplot(111)
    x,y = np.meshgrid(wl[:]/1000,theta_i)
    ax.set_yscale("log")

    p = ax.contourf(y,x,R, cmap=style)
    norm = colors.Normalize(vmin=0, vmax=p.cvalues.max())

    l1 = np.linspace(0.0, p.cvalues.max(), num=6, endpoint=True)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=p.cmap)
    sm.set_array([])

    # We either have always same step (e.g. 0.1) or same number of steps.
    if prop ==1: # Reflectance is always up to 100%
        num_ticks = 6
        ticks_labels = ['0.00', '0.20', '0.40', '0.60', '0.80', '1.00']
    else:
        a = round(p.cvalues.max(), 2)
        if (int(repr(a)[-1])%2 == 1):
            a = a - 1.0/pow(10,len(repr(a))-2)

        if (int(repr(a)[-1]) < 5):
            num_ticks = int(repr(a)[-1])+6 #ATTENTION!!!!!!!!!!!!!!!!
        else:
            num_ticks = int(repr(a)[-1])+1
        labels = np.linspace(0, a, num= num_ticks)
        ticks_labels = [[] for _ in range(num_ticks)]
        for i in range(len(labels)):
            labels[i] = round(labels[i], 2)
            ticks_labels[i] = "%.2f" % labels[i]
        #ticks_labels = [repr(i) for i in labels]


# Creating colobar with 6 ticks and applying user-created labels
    cbar = fig.colorbar(sm)
    tick_locator = ticker.LinearLocator(num_ticks)
    cbar.locator = tick_locator
    cbar.ax.set_yticklabels([''])
    cbar.update_ticks()
    cbar.ax.set_yticklabels(ticks_labels)


    plt.xlabel('Angle, 'r'$\Theta$ (deg)')
    #pl.ylabel('Wavelength, 'r'$\lambda$ ($\mu$m)') #sits on 5*10^0
    plt.text(-11.4, 3.2, 'Wavelength, 'r'$\lambda$ ($\mu$m)', horizontalalignment='left',
            verticalalignment='center', transform=ax.transData, rotation = 'vertical')
    #plt.title(filename[:-11])

    xmajor_ticks = np.arange(0, 105, 15)
    xminor_ticks = np.arange(0, 95, 5)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor = True)
    ax.set_yticklabels([])

    ax.axes.text(-5.8, 18.6, '20', fontsize=10, transform = ax.transData)
    ax.axes.text(-5.8, 9.4, '10', fontsize=10, transform = ax.transData)
    ax.axes.text(-4.9, 4.72, '5', fontsize=10, transform = ax.transData)
    ax.axes.text(-4.9, 0.94, '1', fontsize=10, transform = ax.transData)
    ax.axes.text(-7,  0.47, '0.5', fontsize=10, transform = ax.transData)
    #ax.axes.text(-12,4.7, r'$5\cdot10^{0}$', fontsize=10, transform = ax.transData)


    fig.savefig(foldername+filename+'.png', dpi=600)
    plt.show()
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



