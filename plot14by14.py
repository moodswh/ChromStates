#!/usr/bin/env python

import sys
import math
import os
import gzip
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.stats.mstats import mquantiles

infilename=sys.argv[1]
NumOfStates=int(sys.argv[2])
method=sys.argv[3]
day=sys.argv[4]
StatesMatrix=np.zeros((NumOfStates, NumOfStates))

infile=open(infilename,'r')

for line in infile:
        words=line.rstrip().split()
        i=int(words[0])-1
        j=int(words[1])-1
        counts=float(words[2]) 
        StatesMatrix[i,j]=counts
        StatesMatrix[j,i]=counts

StatesMatrix =np.transpose(StatesMatrix)

vmaxLim=mquantiles(StatesMatrix,[1])[0]
print StatesMatrix.max()
print np.shape(StatesMatrix)
print vmaxLim
fig, ax = plt.subplots()
m = ax.matshow(StatesMatrix, origin="bottom", #norm=colors.LogNorm(),  #norm=colors.SymLogNorm(1),
               cmap="afmhot_r", vmax=vmaxLim)

ax.axhline(-0.5, color="#000000", linewidth=1, linestyle="--")
ax.axvline(-0.5, color="#000000", linewidth=1, linestyle="--")

cb = fig.colorbar(m)
cb.set_label("Obs/Exp Contact Counts")

ax.set_ylim((-0.5, len(StatesMatrix) - 0.5))
ax.set_xlim((-0.5, len(StatesMatrix) - 0.5))

#x=[]
#for i in ax.get_xticks():
#    x.append(i+1)
#ax.set_xticklabels(x)

xy=[]
for i in ax.get_xticks():
    xy.append(int(i+1))
ax.set_xticklabels(xy)
ax.set_yticklabels(xy)

#ax.set_title("%s (%s)" % (organism_name, sample))

#fig.savefig(outfilename+".png",dpi=300)
fig.savefig(method+"_"+day+"ObsExp", figsize=(10,10), dpi=150)

