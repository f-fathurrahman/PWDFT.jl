from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import sys

from operator import sub
def get_aspect(ax):
    # Total figure size
    figW, figH = ax.get_figure().get_size_inches()
    # Axis size on figure
    _, _, w, h = ax.get_position().bounds
    # Ratio of display units
    disp_ratio = (figH * h) / (figW * w)
    # Ratio of data units
    # Negative over negative because of the order of subtraction
    data_ratio = sub(*ax.get_ylim()) / sub(*ax.get_xlim())

    return disp_ratio / data_ratio

plt.style.use("classic")

# Load band structure data
dat = np.loadtxt(sys.argv[1])
Nkpt = dat.shape[0]
Nstates = dat.shape[1] - 1
print("Nkpt = ", Nkpt)
print("Nstates = ", Nstates)

# Load xticks
f = open(sys.argv[1], "r")
Nkpt_spec = int(f.readline().replace("#",""))
print("Nkpt_spec = ", Nkpt_spec)
xticks_pos = [] # XXX change to numpy array ?
xticks_label = []
for i in range(Nkpt_spec):
    l = f.readline().replace("#","")
    xticks_pos.append( float(l.split()[0]) )
    xticks_label.append( l.split()[1] )
print(xticks_pos)
print(xticks_label)

plt.clf()

fig = plt.gcf()
scal = 1.2
fig.set_size_inches(6*scal,8*scal)
#fig.set_size_inches(8*scal,6*scal)  # use landscape

x = dat[:,0]
for ist in range(1,Nstates+1):
    plt.plot( x, dat[:,ist], marker='o' )

plt.xticks(xticks_pos,xticks_label)

ymin, ymax = plt.gca().get_ylim()
print("ymin ymax = ", ymin, ymax)

for ik in range(Nkpt_spec):
    xt = xticks_pos[ik]
    plt.plot( [xt,xt], [ymin,ymax], color="gray" )

plt.gca().set_ylim([ymin,ymax])

filplot = sys.argv[1].replace(".dat",".pdf")
plt.savefig(filplot)

#import os
#os.system("pdfcrop " + filplot + " " + filplot)

