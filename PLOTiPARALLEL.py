import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cmx
import matplotlib as mpl
import pandas as pd
import numpy as np
import csv
from scipy.signal import hilbert
import scipy as sp
from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rcParams
from scipy import signal
import pymp

"""
  import warnings
  warnings.filterwarnings("ignore")
"""

"""
def plot_function(ax1, ax2, ax3, t, v_av, N, neigh, ndriv, k, thresh, length_ratio, rho_t, rho_local, rho_global):
  
  # rho_t vs t
  ax1.plot(t, rho_t)

  # S(t) vs iH(t)
  
  #ax2.scatter(length_ratio, ndriv)
  
  #x = np.linspace(0, N-1, N)
  #for i in range(N): ax3.plot(t, v.iloc[:,i])

  X, Y, Z = (x, v_av, v)

  # Normalize to [0,1]
  norm = plt.Normalize(Z.min(), Z.max())
  colors = cm.viridis(norm(Z))
  rcount, ccount, _ = colors.shape

  ax3 = fig.gca(projection='3d')
  surf = ax3.plot_surface(X, Y, Z, rcount=rcount, ccount=ccount,
                        facecolors=colors, shade=False)
  surf.set_facecolor((0,0,0,0))
  plt.show()

  # Phase Diagram plot
  fig4 = plt.figure(4)
  fig4.canvas.set_window_title('Phase Diagram')
  ax4 = plt.axes()
  ax4 = plt.axes()
  for i in range(N):
    ax4.plot(v[i],w[i])
    ax4.set_xlabel(r"$v$")
    ax4.set_ylabel(r"$w$")
    fig4.tight_layout()

    # Power Spectrum
    fig3 = plt.figure()
    fig3.canvas.set_window_title('Power Spectrum')
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

  
  #fig4.savefig(string + "powerspectrum.pdf")
  
  return fig1#, fig2, fig3
"""
def frequencyCounter(idx, size, stem, driv_locs, v, N, neigh, ndriv, k, thresh):
  
  fig = plt.figure(100*idx, figsize=(8,6), dpi=100)
  ax = Axes3D(fig)
  for i in range(v.shape[1]):
    time_series = v.iloc[:,i]
    T=150/3.0
    N=600
    x = np.linspace(0.0, N*T, N)
    fft = sp.fftpack.fft(time_series)
    psd = np.abs(fft)**2
    psd /= float(max(psd))
    # Sampling window = c_inv
    fftfreq = sp.fftpack.fftfreq(len(psd))
    l = fftfreq > 0
    result = len(l)
    xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
    fields = [N, neigh, ndriv, k, thresh, result]
    with open("FrequencyCounter.csv", 'a') as f:
      writer = csv.writer(f)
      writer.writerow(fields)
    print psd
    
    #xf[xf>0.004]= np.nan
    if float(i) in list(driv_locs): ax.plot(xf, [i+1]*len(xf), psd[:len(xf)], color='red')
    else: ax.plot(xf, [i+1]*len(xf), psd[:len(xf)], color='k')
    #print xf
  #ax.set_xlim3d([0, 0.004])
  ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0.004))
  ax.set_xlabel(r"$f$", fontsize=size, labelpad=20)
  ax.set_zlabel(r"Normalised power", fontsize=size, labelpad=10)
  ax.set_ylabel(r"Neurone", fontsize=size, labelpad=20)

  rcParams.update({'figure.autolayout': True})
  ax.tick_params(axis='both', labelsize=size)
  #ax.legend(loc=0, ncol=2, shadow=False, fontsize=size*0.7)
  fig.tight_layout()
  #plt.gcf().subplots_adjust(bottom=0.5, right = 0.5)


  fig.savefig("Graphs/" + stem + "_fft.pdf", bbox_inches = "tight")
  #ax.clear()
  #plt.show()
  return

def delay_plot(t, v, ax14):
  #fig13 = plt.figure(13)
  #ax13 = plt.axes()

  ax14.grid(linestyle='-')
  ax14.set_axisbelow(True)

  istheredata = False
  cols = [ "red", "blue", "green", "teal", "purple", "orange", "black", "sienna", "violet", "royalblue", "silver", "cyan", "lightcoral", "lightgreen", "darkgreen", "deeppink", "deepskyblue", "gold", "lime", "slategrey", "darkturquoise", "tomato", "navy", "maroon", "tan"]
  for i in range(v.shape[1]):
    time_series = v[:,i]
    m = np.mean(time_series)

    idx = []
    crossing_points = []
    valid_idx = []
    
    for point in time_series:
      if point >= m: idx.append(1)
      else:          idx.append(0)

    for ii in range(1, len(idx)):
      if idx[ii]==1 and idx[ii-1]==0:
        valid_idx.append(ii)
        crossing_points.append(t[ii])

    if len(crossing_points) != 0:
      #print len(crossing_points)
      #ax13.plot(t, time_series)
      #ax13.scatter(crossing_points, time_series[valid_idx])
      taus = np.zeros(len(crossing_points)-1)

      for j in range(len(taus)): taus[j] = crossing_points[j+1] - crossing_points[j]
      
      if len(taus) != 0:
        ax14.scatter(taus[:-1], taus[1:], s=5, c=cols[i])
        istheredata = True
  

  #fig13.savefig("Graphs/Delay Plots/TimeSeries_" + filename + ".pdf")
  t = np.linspace(-100,1000,1101)
  xlim, ylim = ax14.get_xlim(), ax14.get_ylim()
  if xlim>ylim: ylim=xlim
  else:         xlim=ylim
  ax14.plot(t,t, 'k', linewidth=1)
  mpl.rcParams.update({'font.weight': 'bold'})
  ax14.set_xlim(xlim)
  ax14.set_ylim(ylim)
  ax14.set_xlabel(r"$\tau_n$",fontsize=20)
  ax14.set_ylabel(r"$\tau_{n+1}$",fontsize=20)
  ax14.tick_params(labelsize=20)
  
  
  locsx = ax14.get_xticks()
  locsy = ax14.get_yticks()

  """
  order = abs(locsx[-1]-locsx[-2])

  if order<=1e-3:
    ax14.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.4f'))
    ax14.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.4f'))
  elif order<=1e-2:
    ax14.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    ax14.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
  elif order<=1e-1:
    ax14.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    ax14.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
  elif order<=1:
    ax14.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    ax14.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
  elif order<=10:
    ax14.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    ax14.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
  elif order<=100:
    ax14.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    ax14.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
  elif order<=1000:
    ax14.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    ax14.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
  """
  ax14.xaxis.set_major_locator(ticker.MaxNLocator(6))
  ax14.yaxis.set_major_locator(ticker.MaxNLocator(6))

  ax14.set_aspect('equal', adjustable='box')
  
  return istheredata

def scatter6d():
  # Visualizing 6-D mix data using scatter charts
  # leveraging the concepts of hue, size, depth and shape
  fig = plt.figure(figsize=(8, 6))
  t = fig.suptitle('Wine Residual Sugar - Alcohol Content - Acidity - Total Sulfur Dioxide - Type - Quality', fontsize=14)
  ax = fig.add_subplot(111, projection='3d')

  xs = list(wines['residual sugar'])
  ys = list(wines['alcohol'])
  zs = list(wines['fixed acidity'])
  data_points = [(x, y, z) for x, y, z in zip(xs, ys, zs)]

  ss = list(wines['total sulfur dioxide'])
  colors = ['red' if wt == 'red' else 'yellow' for wt in list(wines['wine_type'])]
  markers = [',' if q == 'high' else 'x' if q == 'medium' else 'o' for q in list(wines['quality_label'])]

  for data, color, size, mark in zip(data_points, colors, ss, markers):
      x, y, z = data
      ax.scatter(x, y, z, alpha=0.4, c=color, edgecolors='none', s=size, marker=mark)

  ax.set_xlabel('Residual Sugar')
  ax.set_ylabel('Alcohol')
  ax.set_zlabel('Fixed Acidity')
  
  return

def scatter3d(x, y, z, cs, nfig, colorsMap='jet'):
  """
    cs:It's the values for each 3D point.
    If colors_data is a 3D matrix, with a value for each 3D point, X is a Nx3 matrix of list of
    x,y,z coordinates, so cs=[colors_data[x, y, z] for x,y,z in zip(X[:, 0], X[:, 1], X[:, 2])]
    
  """
  # Work out how to rotate camera before saving
  cm = plt.get_cmap(colorsMap)
  cNorm = mpl.colors.Normalize(vmin=0, vmax=1)
  scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
  fig = plt.figure(nfig, figsize=(10, 8), dpi=150)
  ax = Axes3D(fig)
  ax.scatter(x, y, z, c=scalarMap.to_rgba(cs), s=50)
  #ax.xaxis._axinfo['label']['space_factor'] = 10.0
  #ax.yaxis._axinfo['label']['space_factor'] = 10.0
  #ax.zaxis._axinfo['label']['space_factor'] = 10.0
  scalarMap.set_array(cs)

  if nfig == 5:
    ax.set_xlabel(r"$\frac{m}{N}$", fontsize=15, labelpad=20)
    ax.set_ylabel(r"$\frac{D}{N}$", fontsize=15, labelpad=20)
    ax.set_zlabel(r"$k$", fontsize=15, labelpad=20)
    ax.tick_params(axis='both', pad=10, labelsize=15)
  elif nfig == 6:
    ax.set_xlabel(r"$\mathrm{k}$", fontsize=15, labelpad=20)
    ax.set_ylabel(r"$\frac{D}{N}$", fontsize=15, labelpad=20)
    ax.set_zlabel("Firing threshold", fontsize=15, labelpad=20)
    ax.tick_params(axis='both', pad=10, labelsize=15)
  
  cbar = fig.colorbar(scalarMap, shrink=0.7)
  cbar.set_label(r"$\rho_\mathrm{global}$",size=15)
  cbar.ax.tick_params(labelsize=15)
  
  return fig

def PT_2Dcut(pt_k, pt_ndrivs, pt_thresh, pt_rhos, colorsMap = 'jet'):
  """
    cs:It's the values for each 3D point.
    If colors_data is a 3D matrix, with a value for each 3D point, X is a Nx3 matrix of list of
    x,y,z coordinates, so cs=[colors_data[x, y, z] for x,y,z in zip(X[:, 0], X[:, 1], X[:, 2])]
    
  """
  # Work out how to rotate camera before saving
  cm = plt.get_cmap(colorsMap)
  cNorm = mpl.colors.Normalize(vmin=min(pt_rhos), vmax=max(pt_rhos))
  scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)


  mpl.rcParams.update({'figure.autolayout': True})

  fig1 = plt.figure(8, figsize=(10, 8), dpi=150)
  ax1 = plt.axes()
  fig2 = plt.figure(9, figsize=(10, 8), dpi=150)
  ax2 = plt.axes()
  fig3 = plt.figure(10, figsize=(10, 8), dpi=150)
  ax3 = plt.axes()
  
  ax1.scatter(pt_k, pt_ndrivs, c=scalarMap.to_rgba(pt_rhos), s=50)
  ax2.scatter(pt_ndrivs, pt_thresh, c=scalarMap.to_rgba(pt_rhos), s=50)
  ax3.scatter(pt_k, pt_thresh, c=scalarMap.to_rgba(pt_rhos), s=50)
  #ax.xaxis._axinfo['label']['space_factor'] = 10.0
  #ax.yaxis._axinfo['label']['space_factor'] = 10.0
  #ax.zaxis._axinfo['label']['space_factor'] = 10.0
  scalarMap.set_array(pt_rhos)

  ax1.set_xlabel(r"$k$", fontsize=15, labelpad=20)
  ax1.set_ylabel(r"$\frac{D}{N}$", fontsize=15, labelpad=20)
  ax1.tick_params(axis='both', pad=10, labelsize=15)
  
  ax2.set_xlabel(r"$\frac{D}{N}$", fontsize=15, labelpad=20)
  ax2.set_ylabel("Firing threshold", fontsize=15, labelpad=20)
  ax2.tick_params(axis='both', pad=10, labelsize=15)

  ax3.set_xlabel(r"$\mathrm{k}$", fontsize=15, labelpad=20)
  ax3.set_ylabel(r"Firing Threshold", fontsize=15, labelpad=20)
  ax3.tick_params(axis='both', pad=10, labelsize=15)
  
  cbar1 = fig1.colorbar(scalarMap, shrink=0.7)
  cbar1.set_label(r"$\rho_\mathrm{global}$",size=15)
  cbar1.ax.tick_params(labelsize=15)

  cbar2 = fig2.colorbar(scalarMap, shrink=0.7)
  cbar2.set_label(r"$\rho_\mathrm{global}$",size=15)
  cbar2.ax.tick_params(labelsize=15)

  cbar3 = fig3.colorbar(scalarMap, shrink=0.7)
  cbar3.set_label(r"$\rho_\mathrm{global}$",size=15)
  cbar3.ax.tick_params(labelsize=15)
  
  return fig1, fig2, fig3

def lengthscale(v, driv, driv_locs, N_inv, dim):
  """ Finds ratio of signal death to box size """
  # Test neurones away from the driver
  # Loop over neurones
  # Find most dead neurone
  # Measure distance to nearest driver
  # If many neurones dead, find closest to a driver
  # repeat this for each driver
  double_amplitude = np.zeros(N)
  for neurone in range(N): double_amplitude[neurone] = max(v[:, neurone])-min(v[:, neurone])
  dead_neurones = np.where(double_amplitude == double_amplitude.min())[0]
  
  # Array with all combinations of dead neurones and drivers
  neurone_driver_pair = np.empty([dead_neurones.shape[0]*driv_locs.shape[0], 3])
  
  idx = 0
  for neurone in dead_neurones:
    for driverloc in driv_locs:
      neurone_driver_pair[idx][0] = neurone
      neurone_driver_pair[idx][1] = driverloc
      neurone_driver_pair[idx][2] = 0
    
    neurone_no = neurone_driver_pair[idx][0]
    driver_no  = neurone_driver_pair[idx][1]
    
    # Find a dead neurone that is closest to the driver
    if dim==1:
      dist = np.abs(neurone_no-driver_no)
    
    elif dim==2:
      no_rows    = abs(driver_no-neurone_no) // sqrtN
      no_columns = abs(driver_no-no_rows*sqrtN-neurone_no)
      
      dist = np.sqrt(no_rows*no_rows + no_columns*no_columns)
    
    elif dim==3:
      group_no_of_driver   = driver_no  // cbrtN
      group_pos_of_driver  = driver_no  % cbrtN
      
      group_no_of_neurone  = neruone_no // cbrtN
      group_pos_of_neurone = neurone_no % cbrtN
      
      no_x = abs(group_pos_of_driver - group_pos_of_neurone)
      no_y = abs(group_no_of_driver%cbrtN - group_no_of_neurone%cbrtN)
      no_z = abs((driver_no/cbrtN)%cbrtN-(neurone_no/cbrtN)%cbrtN)
      
      dist = np.sqrt(no_x*no_x + no_y*no_y + no_z*no_z)
    
    neurone_driver_pair[idx][2] = dist
    idx += 1
  
  # Find a dead neurone that is closest to the driver
  l_dead = np.argmin(neurone_driver_pair[:,2])
  
  length_ratio = l_dead*N_inv

  return length_ratio

def find_expip(v, N, neigh, ndriv, k, thresh, neurone):
  """ Utility function for OP calculation """

  data = v[:,neurone]
  mean = np.mean(data)

	#Subtract the mean from data to ensure that Argand diagram centred on 0
  data = np.subtract(data,mean)
  
  #calculates the analytic signal and then exp(i*phi)
  analytic_signal = hilbert(data)
  i_p = np.unwrap(np.angle(analytic_signal))

  return list(map(lambda x: np.exp(1j*x),i_p))

def OP(vFile, v, N, neigh, ndriv, k, thresh):
  """ Calculates global OP """
  expips = []
  
  if v.shape[1] != 0:
    for neurone in range(v.shape[1]): expips.append(find_expip(v, N, neigh, ndriv, k, thresh, neurone))
  else: return -1, -1
  # Find time varying OP and its time average
  orderParam = np.array([np.abs(np.sum(x)/float(N)) for x in zip(*expips)])
  OPavg = np.mean(orderParam)

  """
  # Time varying order parameter
  np.savetxt("/media/jay/0FD90FF80FD90FF8/PROJECTDATA/2DPERIODIC/OP/OP_" + vFile, orderParam,delimiter = ",\n")
  
  # Write global OP for the system to file
  fields = [N, neigh, ndriv, k, thresh, OPavg]
  with open("/media/jay/0FD90FF80FD90FF8/PROJECTDATA/2DPERIODIC/OrderParameter.csv", 'a') as f:
    writer = csv.writer(f)
    writer.writerow(fields)
  """
  return OPavg, orderParam
  
def localOP(v, N, neigh, ndriv, k, thresh):

  """ Calculates local OP """
  data = []
  for neurone in v:
    temp = []
    temp.append(find_expip(v.iloc[:,neurone-1])[0])
    temp.append(find_expip(v.iloc[:,neurone  ])[0])
    temp.append(find_expip(v.iloc[:,neurone+1])[0])
    
    #average local order parameter for each neuron, calculated only from the nearest neighbours
    data.append(np.mean([np.abs(np.sum(x)*0.333333) for x in zip(*temp)]))
  
  fields = [N,neigh,locy,k,thresh]
  with open("LocalOrderParameter.csv", 'a') as f:
    writer = csv.writer(f)
    writer.writerow(fields)
    writer.writerow(data)
  
  return data

""" DOES NOT WORK WITH PARALLELISATION
def F1(ax1, fig, t, rho_t):
  
  # rho_t vs t
  ax1.plot(t, rho_t)
  # Pausing no longer works with parallelisation
  plt.pause(1)

  return fig1

def F2(ax2, fig, length_ratio, ndriv):
  
  ax2.scatter(length_ratio, ndriv)
  plt.pause(0.02)

  return fig

def F3(ax3, fig, t, v):
  
  # Amplitude Plot
  for i in range(N): ax3.plot(t, v.iloc[:,i])
  plt.pause(0.02)

  return fig
"""
def spectrogram():

  fs, data = wavfile.read('../Audio Files/monoguitar.wav')

  # Get first column of wav data
  x = data[:,0]

  # Get data for spectrogram
  f, t, Sxx = signal.spectrogram(x, fs, nperseg=8192)

  # Obtain the melody (highest frequencies that are not overtones)
  # (axis 0 is first column)
  high_f = f[np.argmax(Sxx, axis=0)]

def contour_plot(pt_neighN, pt_ndrivs, ndrivs_list, pt_thresh, rho_matrix, length_ratio, loop_number):
  pi_1 = pt_neighN
  pi_2 = ndrivs_list
  #pi_3 = length_ratio
  #pi_4 = 0
  pi_5 = []

  X = pi_1
  Y = pi_2
  Z = rho_matrix
  
  fig11 = plt.figure(11)
  ax11 = Axes3D(fig11)
  ax11.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)
  cset = ax11.contour(X, Y, Z, zdir='z', offset=-100, cmap=cmx.coolwarm)
  cset = ax11.contour(X, Y, Z, zdir='x', offset=-40, cmap=cmx.coolwarm)
  cset = ax11.contour(X, Y, Z, zdir='y', offset=40, cmap=cmx.coolwarm)
  ax11.set_xlabel(r"$\frac{m}{N}$")
  ax11.set_ylabel(r"$\frac{D}{N}$")
  ax11.set_ylabel(r"$\rho$")

  #fig12 = plt.figure(12)
  #ax12 = plt.axes()
  #ax12.plot(pi_2, pi_5)
  #ax12.set_xlabel(r"$\frac{D}{N}$")
  #ax12.set_ylabel(r"$\frac{\mathrm{Firing Threshold}}{k}$")

  #fig13 = plt.figure(13)
  #ax13 = plt.axes()
  #ax13.plot(pi_3, pi_4)

  #fig14 = plt.figure(14)
  #ax14 = plt.axes()
  #ax14.plot(pi_4, pi_5)

  #fig15 = plt.figure(15)
  #ax15 = plt.axes()
  #ax15.plot(pi_5, pi_1)
  #ax15.set_xlabel(r"$\frac{\mathrm{Firing Threshold}}{k}$")
  #ax15.set_ylabel(r"$\frac{m}{N}$")

  return

def OPvsNoise(op):
  #op = orderparameter csv dataframe
  orderParam = op.iloc[:,-1]
  thresh = op.iloc[:,-2]
  noise = list(map(lambda x: np.sqrt(x*(1-x)),thresh))
######################################################################

if __name__ == "__main__":
  # CHANGES FONT SIZE ON ALL PLOTS
  mpl.rcParams.update({'font.size': 10})
  mpl.rcParams.update({'font.family': 'Times New Roman'})
  mpl.rcParams.update({'figure.autolayout': False})

  rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})
  rc('text', usetex = True)

  fileloc = "/media/jay/0FD90FF80FD90FF83/PROJECTDATA/" + "2DPERIODIC_2/"
  #fileloc = "/media/jay/0FD90FF80FD90FF83/PROJECTDATA/" + "2DPERIODIC/N=25/"

  # Initial values
  # LOOPS ARE INCLUSIVE OF HIGHEST VALUE
  nrootmin = 5
  nrootmax = 6
  neighmin = 1
  neighmax = 2
  drivmin = 4
  k_idx_min = 0
  k_idx_max = 3
  thresh_idx_min = 0
  thresh_idx_max = 3

  timestep = 1.5e-3
  transient_time = 12

  WITHOUT_DEAD = False

  #---------------#
  #---SWITCHBOX---#
  dimension = 2
  #WITHOUT_DEAD = True
  #---------------#

  start_row = int(transient_time/timestep)+1

  size = 20
  
  """
  # Amplitude plot
  fig1 = plt.figure(1, figsize=(8, 6), dpi=100)
  ax1 = plt.axes()

  ax1.set_xlabel(r"$t$", fontsize=size)
  ax1.set_ylabel(r"$\rho(t)$", fontsize=size)
  ax1.tick_params(axis='both', labelsize=size)

  
  # L_death/L vs ndrivers vs Neigh
  fig2 = plt.figure(2, figsize=(12, 10), dpi=100)
  ax2 = plt.axes()

  #ax2.set_zlabel(r"$m$", fontsize=size, labelpad=30)
  ax2.set_xlabel("Length Ratio", fontsize=size, labelpad=30)
  ax2.set_ylabel(r"$D$", fontsize=size, labelpad=30)
  ax2.tick_params(axis='both', pad=10, labelsize=size)

  fig3 = plt.figure(3, figsize=(8, 6), dpi=100)
  ax3 = plt.axes()

  ax3.set_xlabel("Neurone number", fontsize=size)
  ax3.set_ylabel(r"$v_{av}$", fontsize=size)
  ax3.tick_params(axis='both', labelsize=size)
  """
  figs = []
  axes = []
  
  for g in range(thresh_idx_min, thresh_idx_max):
    fig1 = plt.figure(g)
    ax1 = plt.axes()
    figs.append(fig1)
    axes.append(ax1)  
  
  idx = 0
  loop_number = 0

  # Loop over N
  for nroot in range(nrootmin, nrootmax):
    
    N = 1;
    for dim in range(dimension): N *= nroot
    neqn = 2*N
    N_inv = 1/float(N)
    
    # Define a sensible maximum number of drivers
    if N%2==0:
      drivmax = int(N*0.5)
    #neighmax = N*0.5;
    else:
      drivmax = 7
    #neighmax = (N-1)*0.5;
    
    array_length = 20000#(drivmax-1)*(k_idx_max-k_idx_min)*(thresh_idx_max-thresh_idx_min)

    print array_length

    # Shared arrays for phase transition plot
    pt_neighN = pymp.shared.array((array_length,), dtype='float32')
    pt_ndrivs = pymp.shared.array((array_length,), dtype='float32')
    pt_k      = pymp.shared.array((array_length,), dtype='float32')
    pt_thresh = pymp.shared.array((array_length,), dtype='float32')
    pt_rhos   = pymp.shared.array((array_length,), dtype='float32')
    pt_length_ratio   = pymp.shared.array((array_length,), dtype='float32')
    k_rho_matrix = pymp.shared.array((k_idx_max-k_idx_min, array_length), dtype='float32')
    k_list = pymp.shared.array((k_idx_max-k_idx_min,), dtype='float32')
    ndrivs_rho_matrix = pymp.shared.array((drivmax-drivmin, array_length), dtype='float32')
    ndrivs_list = pymp.shared.array((drivmax-drivmin,), dtype='float32')

    # Get Computation Times
    #compTime = pd.read_csv(str(fileloc+"Times.csv"), sep=',',index_col=None, header=None, dtype=np.float64)
    
    # Dimensions of systems
    sqrtN = np.sqrt(N)
    cbrtN = np.cbrt(N)
    
    for neigh in range(neighmin, neighmax):
      for ndriv in range(drivmin, drivmax):
        for k_idx in range(k_idx_min, k_idx_max):
          k = 0.1*k_idx
          with pymp.Parallel(pymp.config.num_threads[0]-2) as p:
            for thresh_idx in p.xrange(thresh_idx_min, thresh_idx_max):
              thresh = thresh_idx*0.1

              idx_private = idx + thresh_idx
              
              # Filenames
              stem = str(N) + "_" + str(neigh) + "_" + str(ndriv) + "_" + str(k) + "00_" + str(thresh)
              vFile     = "v_"     + stem + "00.csv"
              drivFile  = "driv_"  + stem + "00.csv"
              couplFile = "coupl_" + stem + "00.csv"
              
              print vFile, idx_private
              
              v     = pd.read_csv(str(fileloc+"Data/"+vFile    ), sep=',', index_col=None, header=None, dtype=np.float64, skiprows=start_row)
              driv  = pd.read_csv(str(fileloc+"Data/"+drivFile ), sep=',', index_col=None, header=None, dtype=np.float64)
              #coupl = pd.read_csv(str(fileloc+"Data/"+couplFile), sep=',', index_col=None, header=None, dtype=np.float64, skiprows=start_row)
              
              # Get time points, membrane potential and locations of drivers from dataframes
              t         = v.iloc[:,0]
              #driv_locs = driv.iloc[:,0]
              
              if any(v.iloc[:,-1].isna()):
                v_all = np.asarray(v.iloc[:,1:-1])
              else:
                v_all = np.asarray(v.iloc[:,1:])
              
              v = v_all
              
              Nf=N
              
              # Remove dead neurones from analysis
              if WITHOUT_DEAD == True:
                for i in range(N):
                  if np.max(v_all[:,i])-np.min(v_all[:,i]) < 3e-4:
                    v[:,i] = np.nan
                    Nf-=1

                # Remove all nan columns (dead neurones)
                v = v[:,~np.all(np.isnan(v), axis=0)]

              # Get the smallest length scale on which neurones in the system die
              #length_ratio      = lengthscale(v_all, ndriv, driv_locs, N_inv, dimension)
              #rho_global, rho_t = OP(vFile, v, Nf, neigh, ndriv, k, thresh)
              #rho_local         = localOP(v, N, neigh, ndriv, k, thresh)
              rho_local = 0
              
              # Initialise average amplitude array
              v_av = np.zeros(N)

              # Calculate average amplitude of each neurone
              #for i in range(v.shape[1]): v_av[i] = np.sum(v[:,i])/float(len(t))
              
              # Append data to arrays for phase transition plot
              #pt_neighN[idx_private] = neigh/float(v.shape[1])
              #pt_ndrivs[idx_private] = ndriv/float(v.shape[1])
              #pt_k[idx_private]      = k
              #pt_thresh[idx_private] = thresh
              #pt_rhos[idx_private]   = rho_global
              #pt_length_ratio[idx_private] = length_ratio
              
              # ALTER THE NUMBER OF NEIGHBOURS FOR EACH SYSTEM TYPE AND CHANGE CSV NAME FROM C TO MATCH
              ax14 = axes[thresh_idx]
              fig14 = figs[thresh_idx]

              istheredata = delay_plot(t, v, ax14)
              
              if istheredata == True: fig14.savefig("Graphs/Delay Plots/DelayPlot_"  + stem + ".pdf")
              ax14.clear()

              """ DOES NOT WORK WITH PARALLELISATION
              # Send data off for plotting
              #fig1 = F1(ax1, fig1, t, rho_t)
              #plt.close()
              #fig2 = F2(ax2, fig2, length_ratio, ndriv)
              #plt.close()
              #fig3 = F3(ax3, fig3, t, v)
              #plt.close()
              #frequencyCounter(idx, size, stem, driv_locs, v, N, neigh, ndriv, k, thresh)
              """
          
          idx += thresh_idx_max-thresh_idx_min
          
          #k_rho_matrix[k_idx, :] = pt_rhos
          #k_list[k_idx] = k
        
        #ndrivs_rho_matrix[ndriv-1, :] = pt_rhos
        #ndrivs_list[ndriv-1] = ndriv

    # Remove the values for dead neurones (rho_t = -1 from OP())
    """
    if WITHOUT_DEAD == True:
      # Note: Not elegant, but it works!
      #new_pt_neighN = []
      new_pt_ndrivs = []
      new_pt_k      = []
      new_pt_thresh = []
      new_pt_rhos   = []
      #new_pt_length_ratio = []

      for i in range(array_length):
        if pt_rhos[i] != -1:
          #new_pt_neighN.append(pt_neighN[i])
          new_pt_ndrivs.append(pt_ndrivs[i])
          new_pt_k.append(pt_k[i])
          new_pt_thresh.append(pt_thresh[i])
          new_pt_rhos.append(pt_rhos[i])
          #new_pt_length_ratio.append(pt_length_ratio[i])
      
      #pt_neighN = new_pt_neighN
      pt_ndrivs = new_pt_ndrivs
      pt_k      = new_pt_k
      pt_thresh = new_pt_thresh
      pt_rhos   = new_pt_rhos
      #pt_length_ratio = new_pt_length_ratio
    
    # Phase transition plot
    #fig5 = scatter3d(pt_neighN, pt_ndrivs, pt_k, pt_rhos, 5, colorsMap = 'jet')
    fig6 = scatter3d(pt_k, pt_ndrivs, pt_thresh, pt_rhos, 6, colorsMap = 'jet')
    #fig8, fig9, fig10 = PT_2Dcut(pt_k, pt_ndrivs, pt_thresh, pt_rhos, colorsMap = 'jet')
    """
    #fig11 = contour_plot(pt_neighN, pt_ndrivs, ndrivs_list, pt_thresh, ndrivs_rho_matrix, pt_length_ratio, loop_number)
    loop_number += 1

  #fig1.savefig("Graphs/" + stem + "_rho_t.pdf")
  #fig2.savefig("Graphs/" + stem + "_el_death.pdf")
  #fig3.savefig("Graphs/" + stem + "_amplitudeplot.pdf")
  #fig5.savefig("Graphs/" + stem + "_m_D_k.pdf")
  #fig6.savefig("Graphs/2DNONPERIODIC_N=" + str(N) + "_k_D_thresh.pdf")
  #plt.show()

  #plt.close()

