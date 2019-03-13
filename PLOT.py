import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
import pandas as pd
import numpy as np
from scipy.signal import argrelextrema, savgol_filter

from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from random import randint


import warnings
warnings.filterwarnings("ignore")

# CHANGES FONT SIZE ON ALL PLOTS
mpl.rcParams.update({'font.size': 10})
mpl.rcParams.update({'font.family': 'Times New Roman'})
mpl.rcParams.update({'figure.autolayout': False})

rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})
rc('text', usetex = True)

a = pd.read_csv('FHN.csv'	  , sep=',',index_col=None, header=None, dtype=np.float64)
b = pd.read_csv('Coupling.csv', sep=',',index_col=None, header=None, dtype=np.float64)
Ncolumns = a.columns[-1] #outputs the column names used by a

T         = np.array(a[0])
v         = []
w         = []
coupling  = []

for i in range(1, Ncolumns+1, 2):
	v.append(np.array(a[i  ]))
	w.append(np.array(a[i+1]))
	coupling.append(np.array(b[(i-1)/2 + 1]))

# DEFINES GRID DIMENSIONS FOR FIGURES
Nx = 3#int(math.sqrt(Ncolumns/2))
Ny = 3#Nx

#NORMALISE THE TIME AXIS
c = 3.0
c_inv = 1/c
for i in range(len(T)): T[i] = T[i]*c_inv

UNIREVERSE      = False
BI1D            = False
BICIRCULAR      = False
BI2D            = False
BI2DPer         = False
BI3D            = False
BI3DPer         = False
TRANSIENCE      = False
MEANFIELD1D     = False
MEANFIELD1DCIRC = False
MEANFIELD2D     = False
MEANFIELD2D_PER = False
THRESHOLD1D     = False
ALLTOALL        = False

##################################
#---------- SWITCHBOX -----------#
                                 #
# Save files? (0 = No, 1 = Yes)  #
save = 0                         #
                                 #
# System type                    #
#UNIREVERSE      = True          #
#BI1D            = True          #
#BICIRCULAR      = True          #
#BI2D            = True          #
#BI2DPer         = True          #
#BI3D            = True          #
#BI3DPer         = True          #
#MEANFIELD1D     = True          #
#MEANFIELD1DCIRC = True          #
#MEANFIELD2D     = True          #
#MEANFIELD2D_PER = True          #
#THRESHOLD1D     = True          #
ALLTOALL        = True          #
                                 #
# Include Transience?            #
TRANSIENCE = True                #
##################################

# COUPLING TYPE (for saving figures to file)
if UNIREVERSE      == True : path = "uni_reverse"
if BI1D            == True : path = "bi_nearest_neighbour_1D"
if BICIRCULAR      == True : path = "bi_nearest_neighbour_circular"
if BI2D            == True : path = "bi_nearest_neighbours_2D_non_periodic"
if BI2DPer         == True : path = "bi_nearest_neighbours_2D_periodic_inward"
if BI3D            == True : path = "3d_non_periodic"
if BI3DPer         == True : path = "3d_periodic"
if MEANFIELD1D     == True : path = "mean_field_1D"
if MEANFIELD1DCIRC == True : path = "mean_field_1D_circular"
if MEANFIELD2D     == True : path = "mean_field_2D_non_periodic"
if MEANFIELD2D_PER == True : path = "mean_field_2D_periodic"
if THRESHOLD1D     == True : path = "threshold_1D"
if ALLTOALL        == True : path = "bi_all_to_all"

end = ""

# REMOVES TRANSIENT BEHAVIOUR
if TRANSIENCE == False:
  transient_time = 11.0;

  transient_idx = (np.abs(T - transient_time)).argmin()

  T = T[transient_idx:]

  for i in range(len(v)):
    v[i]         =         v[i][transient_idx:]
    w[i]         =         w[i][transient_idx:]
    coupling[i]  =  coupling[i][transient_idx:]
else:
  end = "_transient"

# Generic spacing to start out with
xtick_spacing = 10
ytick_spacing = 0.1

if UNIREVERSE == True:
  xtick_spacing = 10
  ytick_spacing = 1.0
elif BI1D == True:
  xtick_spacing = 10
  ytick_spacing = 0.5

# v-T plot
fig1, axfig1 = plt.subplots(Nx, Ny, sharex=True, sharey=False)
fig1.canvas.set_window_title('T-v Plot')
fig1.subplots_adjust(hspace=0.2, wspace=0.8)
idx = 0
for i in range(Nx):
  for j in range(Ny):
    """
    if UNIREVERSE == True:
      if   idx == 3 : ytick_spacing = 0.1
      elif idx == 4 : ytick_spacing = 0.02
      elif idx == 5 : ytick_spacing = 0.002
      elif idx == 6 : ytick_spacing = 0.0005
      elif idx == 7 : ytick_spacing = 0.00005
      elif idx == 8 : ytick_spacing = 0.00001
    elif BI1D == True:
      if   idx == 1 : ytick_spacing = 0.1
      elif idx == 2 : ytick_spacing = 0.02
      elif idx == 3 : ytick_spacing = 0.005
      elif idx == 4 : ytick_spacing = 0.0005
      elif idx == 5 : ytick_spacing = 0.0001
      elif idx == 6 : ytick_spacing = 0.00002
      elif idx == 7 : ytick_spacing = 0.000002
      elif idx == 8 : ytick_spacing = 0.02
    """
    #ax = fig1.add_subplot(Nx,Ny,i+1)
    axfig1[i][j].plot(T,v[idx])
    #ax[i][j].grid()
    #ax[i][j].xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
    #ax[i][j].yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
    #ax[i][j].ticklabel_format(useOffset=True)
    #ax[i][j].yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.5f'))
    #ax[i][j].set_xlim(left=0)
    #ax[i][j].set_ylim([-0.925,-0.923])
    idx += 1
fig1.text(0.5, 0.01, r"Normalised time", ha='center')
fig1.text(0.005, 0.5, r'$v(t)$', va='center', rotation='vertical')


"""
# Hide x labels and tick labels for top plots and y ticks for right plots.
for a in ax.flat:
  a.label_outer()
"""
"""
# w-T plot
fig1 = plt.figure()
fig1.canvas.set_window_title('v-T Plot')
fig1.subplots_adjust(hspace=0.4, wspace=0.4)
for i in range(Ncolumns/2):
  ax = fig2.add_subplot(Nx,Ny,i+1)
  ax.plot(T,v[i])
fig1.savefig("v-T.pdf")
"""

if UNIREVERSE == True:
  xtick_spacing = 10
  ytick_spacing = 0.4
elif BI1D == True:
  xtick_spacing = 10
  ytick_spacing = 0.5

# w-T plot
fig2, ax = plt.subplots(Nx, Ny, sharex=True, sharey=False)
fig2.canvas.set_window_title('T-w Plot')
fig2.subplots_adjust(hspace=0.2, wspace=0.8)
idx = 0
for i in range(Nx):
  for j in range(Ny):
    """
    if UNIREVERSE == True:
      if   idx == 3 : ytick_spacing = 0.05
      elif idx == 4 : ytick_spacing = 0.005
      elif idx == 5 : ytick_spacing = 0.001
      elif idx == 6 : ytick_spacing = 0.0001
      elif idx == 7 : ytick_spacing = 0.00002
      elif idx == 8 : ytick_spacing = 0.000002
    elif BI1D == True:
      if   idx == 1 : ytick_spacing = 0.05
      elif idx == 2 : ytick_spacing = 0.005
      elif idx == 3 : ytick_spacing = 0.001
      elif idx == 4 : ytick_spacing = 0.0002
      elif idx == 5 : ytick_spacing = 0.00005
      elif idx == 6 : ytick_spacing = 0.00001
      elif idx == 7 : ytick_spacing = 0.000001
      elif idx == 8 : ytick_spacing = 0.0000002
    """
    #ax = fig1.add_subplot(Nx,Ny,i+1)
    ax[i][j].plot(T,w[idx])
    #ax[i][j].grid()
    #ax[i][j].xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
    #ax[i][j].yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
    #ax[i][j].ticklabel_format(useOffset=True)
    #ax[i][j].yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.6f'))
    #ax[i][j].yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
    #ax[i][j].set_xlim(left=0)
    idx += 1
fig2.text(0.5, 0.01, r"Normalised time", ha='center')
fig2.text(0.005, 0.5, r"$w(t)$", va='center', rotation='vertical')

"""
# Hide x labels and tick labels for top plots and y ticks for right plots.
for a in ax.flat:
  a.label_outer()
  """
"""
# w-T plot
fig2 = plt.figure()
fig2.canvas.set_window_title('w-T Plot')
fig2.subplots_adjust(hspace=0.4, wspace=0.4)
for i in range(Ncolumns/2):
  ax = fig2.add_subplot(Nx,Ny,i+1)
  ax.plot(T,w[i])
fig2.savefig("w-T.pdf")
"""
xtick_spacing = 10
# T-coupling plot
fig3, ax = plt.subplots(Nx, Ny, sharex=True, sharey=False)
fig3.canvas.set_window_title('T-coupling Plot')
fig3.subplots_adjust(hspace=0.1, wspace=0.5)
idx = 0
for i in range(Nx):
  for j in range(Ny):
    #ax = fig1.add_subplot(Nx,Ny,i+1)
    ax[i][j].plot(T,coupling[idx])
    #ax[i][j].xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
    idx += 1
fig3.text(0.5, 0.01, r"Normalised time", ha='center')
fig3.text(0.005, 0.5, r"$I_\mathrm{eff}(t)$", va='center', rotation='vertical')

"""
# Hide x labels and tick labels for top plots and y ticks for right plots.
for a in ax.flat:
  a.label_outer()
"""
"""
# coupling-T plot
fig3 = plt.figure()
fig3.canvas.set_window_title('coupling-T Plot')
fig3.subplots_adjust(hspace=0.4, wspace=0.4)
for i in range(Ncolumns/2):
  ax = fig3.add_subplot(Nx,Ny,i+1)
  ax.plot(T,coupling[i])
fig3.savefig("Coupling-T.pdf")
"""


# v-coupling plot
fig4, ax = plt.subplots(Nx, Ny, sharex=False, sharey=False)
fig4.canvas.set_window_title('v-coupling Plot')
fig4.subplots_adjust(hspace=0.3, wspace=0.5)
idx = 0
for i in range(Nx):
  for j in range(Ny):
    #ax = fig1.add_subplot(Nx,Ny,i+1)
    ax[i][j].plot(v[idx], coupling[idx])
    #ax[i][j].yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
    idx += 1
fig4.text(0.5, 0.04, r"$v$", ha='center')
fig4.text(0.04, 0.5, r"$I_\mathrm{eff}(t)$", va='center', rotation='vertical')

"""
# Hide x labels and tick labels for top plots and y ticks for right plots.
for a in ax.flat:
  a.label_outer()
"""
"""
fig4 = plt.figure()
fig4.canvas.set_window_title('coupling-v Plot')
fig4.subplots_adjust(hspace=0.4, wspace=0.4)
for i in range(Ncolumns/2):
  ax = fig4.add_subplot(Nx,Ny,i+1)
  ax.plot(v[i],coupling[i])
fig4.savefig("Coupling-v.pdf")

"""
xtick_spacing = 10
ytick_spacing = 1.0

# w-coupling plot
fig5, ax = plt.subplots(Nx, Ny, sharex=False, sharey=False)
fig5.canvas.set_window_title('w-coupling Plot')
fig5.subplots_adjust(hspace=0.3, wspace=0.5)
idx = 0
for i in range(Nx):
  for j in range(Ny):
    if BI2D == True:
      if   idx == 0 : xtick_spacing = 0.02
      elif idx == 1 : xtick_spacing = 0.002
      elif idx == 2 : xtick_spacing = 0.0004
      elif idx == 3 : xtick_spacing = 0.002
      elif idx == 4 : xtick_spacing = 0.001
      elif idx == 5 : xtick_spacing = 0.00015
      elif idx == 6 : xtick_spacing = 0.0004
      elif idx == 7 : xtick_spacing = 0.00015
      elif idx == 8 : xtick_spacing = 0.00004
    #ax = fig1.add_subplot(Nx,Ny,i+1)
    ax[i][j].plot(w[idx],coupling[idx],"k")
    ax[i][j].xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
    #ax[i][j].set_xlim(left=0)
    idx += 1
fig5.text(0.5, 0.02, r"$w(t)$", ha='center')
fig5.text(0.02, 0.5, r"$I_\mathrm{eff}(t)$", va='center', rotation='vertical')

"""
# Hide x labels and tick labels for top plots and y ticks for right plots.
for a in ax.flat:
  a.label_outer()
"""
"""
# coupling-w plot
fig5 = plt.figure()
fig5.canvas.set_window_title('coupling-w Plot')
fig5.subplots_adjust(hspace=0.4, wspace=0.4)
for i in range(Ncolumns/2):
  ax = fig5.add_subplot(Nx,Ny,i+1)
  ax.plot(w[i],coupling[i])
fig5.savefig("Coupling-w.pdf")
"""

# I-v-w plot
fig6 = plt.figure()
fig6.canvas.set_window_title('I-v-w Plot')
ax = fig6.add_subplot(1,1,1, projection='3d')
for i in range(Ncolumns/2):
  ax.plot(w[i],v[i],coupling[i])
ax.set_xlabel(r"$w$")
ax.set_ylabel(r"$v$")
ax.set_zlabel(r"$I_\mathrm{eff}(t)$")
fig6.tight_layout()


# I-v-t plot
fig7 = plt.figure()
fig7.canvas.set_window_title('I-v-t Plot')
ax = fig7.add_subplot(1,1,1, projection='3d')
for i in range(Ncolumns/2):
  ax.plot(T,v[i],coupling[i])
ax.set_xlabel(r"Normalised time")
ax.set_ylabel(r"$v$")
ax.set_zlabel(r"$I_\mathrm{eff}(t)$")
fig7.tight_layout()

# I-v-t plot
fig7 = plt.figure()
fig7.canvas.set_window_title('v-w-t Plot')
ax = fig7.add_subplot(1,1,1, projection='3d')
for i in range(Ncolumns/2):
  ax.plot(v[i],w[i],coupling[i])
ax.set_xlabel(r"Normalised time")
ax.set_ylabel(r"$v$")
ax.set_zlabel(r"$I_\mathrm{eff}(t)$")
fig7.tight_layout()


# Phase Diagram plot
fig8 = plt.figure()
fig8.canvas.set_window_title('Phase Diagram')
ax = plt.axes()
for i in range(Ncolumns/2):
  ax.plot(v[i],w[i])
ax.set_xlabel(r"$v$")
ax.set_ylabel(r"$w$")
fig8.tight_layout()


# Power Spectrum
fig9 = plt.figure()
fig9.canvas.set_window_title('Power Spectrum')
ax = plt.axes()
for i in range(Ncolumns/2):
  ps = np.abs(np.fft.fft(v[i]))**2
  #psmax = float(ps.max())
  #for j in range(len(ps)): ps[j] = ps[j]/psmax
  
  time_step = 1 / 800.0
  freqs = np.fft.fftfreq(v[i].size, time_step)
  idx = np.argsort(freqs)
  ax.plot(freqs[idx],ps[idx])
ax.set_xlabel(r"$Frequency$")
ax.set_ylabel(r"$Normalised Power$")
ax.set_xlim(-0.01, 0.1)


# Delay plot
fig10 = plt.figure()
fig10.canvas.set_window_title('Delay Coordinate Plot')
ax = plt.axes()
length = len(v[0])
tolerance = T[3]-T[0]

# Colours for delay plot
colors = ['k', 'r', 'b', 'g']
"""
colors = []
for i in range(Ncolumns/2):
  colors.append('%06X' % randint(0, 0xFFFFFF))
"""
# Tolerance should be difference between points in the time domain
for i in range(Ncolumns/2):
  times = []
  taus = []
  
  # Smooth the data before finding extrema as mentioned in:
  # https://stackoverflow.com/questions/35282456/find-local-maximums-in-numpy-array
  # https://stackoverflow.com/questions/31070563/find-all-local-maxima-and-minima-when-x-and-y-values-are-given-as-numpy-arrays?noredirect=1&lq=1
  # Use a Savitzky-Golay filter
  # No idea what numbers to use here
  # If polynomial order is not 1, maxima and minima get jumbled up between arrays
  # If window size is too small, number of maxima and number of minima do differ by more than 1 (they should be different by either 0 or 1)
  #vhat = np.convolve(v) #savgol_filter(v[i], 1001, 2) # window size 1001, polynomial order 1
  
  #np.r_[True, a[1:] < a[:-1]] & numpy.r_[a[:-1] < a[1:], True]
  vhat = v[i]
  # list(set()) removes duplicate values near maxima and minima
  maxima = vhat[argrelextrema(vhat, np.greater)[0]]
  minima = vhat[argrelextrema(vhat, np.less)[0]]
  
  # Find an appropraite mean which all peaks cross
  if len(maxima) != 0 and len(minima) != 0:
    mean = 0.5*(min(maxima) + max(minima))
  #mean = (sum(vhat)/len(vhat))# + sum(minima)/len(minima))
  else: mean = 0.0
  
  # Find mean values for the v-T plot
  if i < Nx*Ny: axfig1[i/Nx][i%Nx].plot(T,[mean]*len(T))
  
  # Find crossing times
  # Double the tolerance on how close values are to mean if not enough crossings are found
  # NOTE: THIS IS BOTCHED -- the increase in tolerance is arbitrary
  # Would be better to remove while loop if possible
  tol = tolerance

  for j in range(1,length-1):
    if vhat[j] > mean - tol and vhat[j] < mean + tol and (len(times) == 0 or T[j] > times[-1] + c_inv-0.1) and vhat[j+1] < vhat[j]:
      #print "MEAN ", mean
      #print "V ", vhat[j]
      times.append(T[j])

  # Find taus
  for j in range(len(times)-1):
    taus.append(times[j+1]-times[j])
  taus = np.array(taus)

  print "Neurone : ", i
  print "Mean : ", mean
  print "Maxima : ", maxima
  print "Minima : ", minima
  #print "Smallest maximum : ", min(maxima)
  #print "Largest minimum : ", max(minima)
  print "Taus : \n", taus, "\n"

ax.scatter(taus[:-1],taus[1:], s=1)
ax.set_xlabel(r"$\tau_i$")
ax.set_ylabel(r"$\tau_{i+1}$")
#ax.set_xlim(-0.01, 0.1)


if save == 1:
  fig1.savefig(path + "/v-T" + end + ".pdf")
  fig2.savefig(path + "/w-T" + end + ".pdf")
  fig3.savefig(path + "/Coupling-T" + end + ".pdf")
  fig4.savefig(path + "/Coupling-v" + end + ".pdf")
  fig5.savefig(path + "/Coupling-w" + end + ".pdf")
  fig6.savefig(path + "/I-v-w" + end + ".pdf")
  fig7.savefig(path + "/I-v-t" + end + ".pdf")
  fig8.savefig(path + "/phasediagram" + end + ".pdf")
  fig9.savefig(path + "/powerspectrum" + end + ".pdf")
  fig10.savefig(path + "/delayplot" + end + ".pdf")

fig1.show() #vT
#fig2.show() #wT
#fig3.show() #IT
#fig4.show() #Iv
#fig5.show() #Iw
#fig6.show() #Ivw
#fig7.show() #Ivt
#fig8.show() #phase diagram
#fig9.show() #power spectrum
#fig10.show() #delay plot
