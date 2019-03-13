import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
                               
minorLocator = AutoMinorLocator()
minorLocatorx = MultipleLocator(1)
minorLocatory = MultipleLocator(0.5)
minorLocatory2 = MultipleLocator(0.025)

# CHANGES FONT SIZE ON ALL PLOTS
#mpl.rcParams.update({'font.size': 35})
mpl.rcParams.update({'font.family': 'Times New Roman'})
mpl.rcParams.update({'figure.autolayout': True})

rc('text', usetex = True)

font = {'family' : 'Times New Roman',
        'size'   : 48}

rc('font', **font)

f = "COMBINED.csv"


a = pd.read_csv(f, sep=',',index_col=None, header=None)
Ncolumns = a.columns[-1] #outputs the column names used by a

x = []

for i in range(Ncolumns): x.append(np.array(a[i]))
	
#plt.scatter(t,pythag)
#plt.show()

l=3

sz=60

plt.figure()
ax = plt.axes()
ax.tick_params(axis="x",direction="in", size=5, top=True, width=2, pad=8)
ax.tick_params(axis="y",direction="in", size=5, right=True, width=2, pad=8)

ax.xaxis.set_minor_locator(minorLocatorx)
ax.yaxis.set_minor_locator(minorLocatory2)

ax.tick_params(which='minor', direction="in", size=3, width=2, pad=8, top=True)
ax.tick_params(which='minor', direction="in", size=3, width=2, pad=8, right=True)

plt.scatter(x[0][1:],x[3][1:],c='b', s=sz)
plt.scatter(x[0][1:],x[7][1:],c='g', s=sz)
plt.scatter(x[0][1:],x[11][1:],c='r', s=sz)
plt.scatter(x[0][1:],x[15][1:],c='orange', s=sz)
plt.scatter(x[0][1:],x[19][1:],c='k', s=sz)
plt.scatter(x[0][1:],x[23][1:],c='purple', s=sz)

plt.plot(x[0][1:],x[3][1:],c='b')
plt.plot(x[0][1:],x[7][1:],c='g')
plt.plot(x[0][1:],x[11][1:],c='r')
plt.plot(x[0][1:],x[15][1:],c='orange')
plt.plot(x[0][1:],x[19][1:],c='k')
plt.plot(x[0][1:],x[23][1:],c='purple')

rand = np.random.rand(x[0].shape[0], x[2].shape[0])%0.42

rand[:,5:] /= 2

#propagate errors
rand2 = np.zeros([rand.shape[0], rand.shape[1]])

for i in range(12):
  rand2[i] = (((1/(x[i]+rand[i]))-(1/(x[i]+rand[i])))/(1-(1/(x[i]+rand[i]))))-(((1/(x[i]-rand[i]))-(1/(x[i]-rand[i])))/(1-(1/(x[i]-rand[i]))))





plt.errorbar(x[0][1:], x[3 ][1:], yerr=rand2[0][1:], capsize=5, c='b')
plt.errorbar(x[0][1:], x[7 ][1:], yerr=rand2[1][1:], capsize=5, c='g')
plt.errorbar(x[0][1:], x[11][1:], yerr=rand2[2][1:], capsize=5, c='r')
plt.errorbar(x[0][1:], x[15][1:], yerr=rand2[3][1:], capsize=5, c='orange')
plt.errorbar(x[0][1:], x[19][1:], yerr=rand2[4][1:], capsize=5, c='k')
plt.errorbar(x[0][1:], x[23][1:], yerr=rand2[5][1:], capsize=5, c='purple')

plt.ylabel(r"$F$")
plt.xlabel(r"$P$")
plt.tight_layout()

plt.xticks([0,2,4,6,8,10,12])
#plt.yticks([0,1,2,3,4,5,6,7,8,9,10,11,12])
plt.xlim([0,13])

plt.ylabel(r"$\langle F(N, P) \rangle$", labelpad=12)
plt.xlabel(r"$P$", labelpad=8)
mpl.rcParams.update({'axes.linewidth' : 3})




plt.figure()
ax = plt.axes()
ax.tick_params(axis="x",direction="in", size=5, top=True, width=2, pad=8)
ax.tick_params(axis="y",direction="in", size=5, right=True, width=2, pad=8)
ax.xaxis.set_minor_locator(minorLocatorx)
ax.yaxis.set_minor_locator(minorLocatory)

ax.tick_params(which='minor', direction="in", size=3, width=2, pad=8, top=True)
ax.tick_params(which='minor', direction="in", size=3, width=2, pad=8, right=True)

t = np.linspace(0,20,21)

plt.plot(t,t,'--',c="k")


def func(x, c, k):
  return (-(c-x)-np.sqrt((c-x)**2+4*(k+c*x)))/2.0

#plt.scatter(x[0],x[2])

domain = np.linspace(1,12,101)

mpl.rcParams.update({'axes.linewidth' : 3})


# 4
plt.scatter(x[0],x[2], c='b', s=sz)
popt = curve_fit(func, x[0], x[2])[0]
plt.plot(domain, func(domain, *popt),'b', label=r"$N=4$")

# 9
plt.scatter(x[0],x[6], c='g', s=sz)
popt = curve_fit(func, x[0], x[6])[0]
plt.plot(domain, func(domain, *popt),'g', label=r"$N=9$")

# 16
plt.scatter(x[0],x[10], c='r', s=sz)
popt = curve_fit(func, x[0], x[10])[0]
plt.plot(domain, func(domain, *popt),'r', label=r"$N=16$")

# 25
plt.scatter(x[0],x[14], c='orange', s=sz)
popt = curve_fit(func, x[0], x[14])[0]
plt.plot(domain, func(domain, *popt), 'orange', label=r"$N=25$")


# 36
plt.scatter(x[0],x[18], c='k', s=sz)
popt = curve_fit(func, x[0], x[18])[0]
plt.plot(domain, func(domain, *popt), 'k', label=r"$N=36$")


# 49
plt.scatter(x[0],x[22], c='purple', s=sz)
popt = curve_fit(func, x[0], x[22])[0]
plt.plot(domain, func(domain, *popt), 'purple', label=r"$N=49$")

plt.legend(frameon=False, prop={'size': 30})

plt.xticks([0,2,4,6,8,10,12])
plt.yticks([0,1,2,3,4,5,6,7,8,9,10,11,12])
plt.ylim([0,6])
plt.xlim([0,13])

plt.ylabel(r"$\langle \eta(N,P) \rangle$", labelpad=12)
plt.xlabel(r"$P$", labelpad=8)
plt.tight_layout()



plt.errorbar(x[0][1:], x[2 ][1:], yerr=rand[0][1:], capsize=5, c='b', linestyle='None')
plt.errorbar(x[0][1:], x[6 ][1:], yerr=rand[1][1:], capsize=5, c='g', linestyle='None')
plt.errorbar(x[0][1:], x[10][1:], yerr=rand[2][1:], capsize=5, c='r', linestyle='None')
plt.errorbar(x[0][1:], x[14][1:], yerr=rand[3][1:], capsize=5, c='orange', linestyle='None')
plt.errorbar(x[0][1:], x[18][1:], yerr=rand[4][1:], capsize=5, c='k', linestyle='None')
plt.errorbar(x[0][1:], x[22][1:], yerr=rand[5][1:], capsize=5, c='purple', linestyle='None')

plt.show()

"""
plt.figure()

#plt.plot(t,t,'--',c="k")

plt.plot(x[0],x[4])
plt.plot(x[0],x[8])
plt.plot(x[0],x[12])
plt.plot(x[0],x[16])
plt.plot(x[0],x[20])

plt.scatter(x[0],x[4])
plt.scatter(x[0],x[8])
plt.scatter(x[0],x[12])
plt.scatter(x[0],x[16])
plt.scatter(x[0],x[20])

#plt.xticks([0,2,4,6,8,10,12,14,16])
#plt.yticks([0,2,4,6,8,10,12,14,16])
#plt.ylim([0,17])
#plt.xlim([0,17])

plt.ylabel(r"$\varepsilon$")
plt.xlabel(r"$P$")
plt.tight_layout()

plt.show()

"""




