from __future__ import division
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.patches as patches


#=============INPUT==================
n = 3 # N layers
m = ['$SiO_{2}$', 'Si', '$SiO_{2}$']           #Names and order of layers
d = [200, 400, 200]                       #thickness in [nm]

c = ['#e9ebee','#b5e5da','#e9ebee', '#b5e5da'] #colors (optional)
r = [6.0, 2.0, 6.0]                            #refractive index (optional)
h = ['','/','','/']                            # patterns (optional)
#===================================



scale = n*420.0                    #to fit the figure
st = np.zeros((1,n), dtype=float)  #auxiliary array
for i in range(n-1):
    st[0][i+1] = st[0][i]+ d[i]

fig = pl.figure()
ax = fig.add_subplot(111)
# Draw squares
for p in [
        patches.Rectangle(

           (0.25, 0.02 +st[0][i]/scale),
           0.5,
           d[i]/scale,
           edgecolor = 'black',
           facecolor = c[i],
           #alpha = r[i]/(r[0]+1.5),
           hatch = h[i],
           linewidth= 2.0,
           )
        for i in range(n)]:
            ax.add_patch(p)

params = {'mathtext.default': 'regular' }
pl.rcParams.update(params)
# Draw text (not scalable)
for i in range(n):
    ax.text(0.5, 0.02+(st[0][i]+d[i]/2.0)/scale,
    m[i]+' '+ "("+str(d[i])+' nm)',
    fontsize=15,
    color='black',
    horizontalalignment='center',
    verticalalignment='center',
    transform=ax.transAxes)


# Draw lines
a = 0.02+(st[0][-1]-d[-1])/scale-0.05
arrow1 = ax.arrow(0.5, 1.0, 0.0, -a, head_length= 0.05, head_width = 0.03,
         linewidth=0.5)
ax.set_axis_off()
pl.show()


# NOTES:
#1. to use latin/greek letters: r'$\alpha_{2}$'

# Sources:
#http://matthiaseisen.com/pp/patterns/p0203/
#https://matplotlib.org/api/_as_gen/matplotlib.patches.Rectangle.html
#https://stackoverflow.com/questions/27698377/how-do-i-make-sans-serif-superscript-or-subscript-text-in-matplotlib
#https://matplotlib.org/examples/pylab_examples/arrow_demo.html
#https://matplotlib.org/users/index_text.html
