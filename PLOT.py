import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
import pandas as pd
import numpy as np

from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

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
Nx = 3
Ny = 3

#NORMALISE THE TIME AXIS
c = 3.0
for i in range(len(T)): T[i] = T[i]/c

UNIREVERSE      = False
BI1D            = False
BICIRCULAR      = False
BI2D            = False
BI2DPer         = False
TRANSIENCE      = False
MEANFIELD1D     = False
MEANFIELD2D_PER = False
THRESHOLD1D     = False

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
#MEANFIELD1D     = True          #
#MEANFIELD2D_PER = True          #
THRESHOLD1D     = True          #
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
if MEANFIELD1D     == True : path = "mean_field_1D"
if MEANFIELD2D_PER == True : path = "mean_field_2D_periodic"
if THRESHOLD1D     == True : path = "threshold_1D"

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
ytick_spacing = 1.0

if UNIREVERSE == True:
  xtick_spacing = 10
  ytick_spacing = 1.0
elif BI1D == True:
  xtick_spacing = 10
  ytick_spacing = 0.5

# v-T plot
fig1, ax = plt.subplots(Nx, Ny, sharex=True, sharey=False)
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
    ax[i][j].plot(T,v[idx])
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

fig1.show() #vT
#fig2.show() #wT
#fig3.show() #IT
#fig4.show() #Iv
#fig5.show() #Iw
#fig6.show() #Ivw
#fig7.show() #Ivt
#fig8.show() #phase diagram
#fig9.show() #power spectrum
