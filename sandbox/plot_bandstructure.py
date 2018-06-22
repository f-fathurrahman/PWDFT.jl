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

x = dat[:,0]
for ist in range(1,Nstates+1):
    plt.plot( x, dat[:,ist], marker='o' )

plt.xticks(xticks_pos,xticks_label)

ymin, ymax = plt.gca().get_ylim()
print("ymin ymax = ", ymin, ymax)

for ik in range(Nkpt_spec):
    x = xticks_pos[ik]
    plt.plot( [x,x], [ymin,ymax], color="gray" )

#w = plt.xlim()[1] - plt.xlim()[0]
#h = plt.ylim()[1] - plt.ylim()[0]
#aspect_orig = get_aspect(plt.axes())
#h_w = h/w
#print("h_w = ", h_w)
#plt.axes().set_aspect(h_w/aspect_orig*aspect_orig)

filplot = sys.argv[1].replace(".dat",".pdf")
plt.savefig(filplot)

#import os
#os.system("pdfcrop " + filplot + " " + filplot)

