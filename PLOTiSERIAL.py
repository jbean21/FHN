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
import matplotlib as mpl

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

def delay_plot(v):
  for i in range(v.shape[1]):
    time_series = v.iloc[:,i]
    m = time_series/len(time_series)

    valid_idx = nearest(time_series, m)
    crossing_points = time_series[valid_idx]

    if len(crossing_points != 0):
      #print len(crossing_points)
      plt.plot(t, v)
      plt.scatter(t[valid_idx], crossing_points)
      plt.show()
      taus = np.zeros(len(crossing_points)-1)

      for j in range(len(taus)):
        taus[j] = crossing_points[j+1] - crossing_points[j]
      
      plt.figure(10)
      plt.scatter(taus[:-1], taus[1:])
  plt.show()

  return

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
  cNorm = mpl.colors.Normalize(vmin=min(cs), vmax=max(cs))
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
  for neurone in range(N): double_amplitude[neurone] = max(v.iloc[:, neurone])-min(v.iloc[:, neurone])
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
  data = v[neurone]
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
    for neurone in v: expips.append(find_expip(v, N, neigh, ndriv, k, thresh, neurone))
  
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

def F1(ax1, fig, t, rho_t, ndriv):
  
  # rho_t vs t
  if ndriv == 4:
    #c='k'
    zorder=1
    l = r"D=4"
  elif ndriv == 6:
    #c='C0'
    zorder=0
    l = r"D=6"
  elif ndriv == 11:
    #c='-.k'
    zorder=0
    l = r"D=11"
  elif ndriv == 5:
    #c='--k'
    zorder=1
    l = r"D=5"
  
  ax1.plot(t[8334:22885], rho_t[8334:22885], c='k', zorder=zorder, label=l)
  ax1.xaxis.set_major_locator(plt.MaxNLocator(7))
  #plt.pause(1)
  #plt.legend(frameon=False, prop={'size': 30})

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

def spectrogram():

  fs, data = wavfile.read('../Audio Files/monoguitar.wav')

  # Get first column of wav data
  x = data[:,0]

  # Get data for spectrogram
  f, t, Sxx = signal.spectrogram(x, fs, nperseg=8192)

  # Obtain the melody (highest frequencies that are not overtones)
  # (axis 0 is first column)
  high_f = f[np.argmax(Sxx, axis=0)]


######################################################################

if __name__ == "__main__":
  # CHANGES FONT SIZE ON ALL PLOTS
  """
  mpl.rcParams.update({'font.size': 20})
  mpl.rcParams.update({'font.family': 'Times New Roman'})
  mpl.rcParams.update({'figure.autolayout': False})

  rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})
  rc('text', usetex = True)
  """
  
  mpl.rcParams.update({'font.family': 'Times New Roman'})
  mpl.rcParams.update({'figure.autolayout': True})
  mpl.rcParams.update({'axes.linewidth' : 3})

  rc('text', usetex = True)

  font = {'family' : 'Times New Roman',
          'size'   : 40}

  rc('font', **font)
  rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})
  rc('text', usetex = True)

  fileloc = "/media/jay/0FD90FF80FD90FF83/PROJECTDATA/"+ "2DPERIODIC_2/"

  # Initial values
  # LOOPS ARE INCLUSIVE OF HIGHEST VALUE
  nrootmin = 5
  nrootmax = 6
  neighmin = 1
  neighmax = 2
  drivmin = 1
  k_idx_min = 1
  k_idx_max = 2
  thresh_idx_min = 3
  thresh_idx_max = 4

  timestep = 1.5e-3
  transient_time = 14

  #------------#
  dimension = 2
  #------------#

  start_row = int(transient_time/timestep)+1

  # Arrays for phase transition plot
  pt_neighN = []
  pt_ndrivs = []
  pt_k      = []
  pt_thresh = []
  pt_rhos   = []

  idx = 0
  size = 50

  WITHOUT_DEAD = True

  """
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

  # Loop over N
  for nroot in range(nrootmin, nrootmax):
    
    N = 1;
    for dim in range(dimension): N *= nroot
    neqn = 2*N
    N_inv = 1/float(N)
    
    # Define a sensible maximum number of drivers
    if N%2==0:
      drivmax = int(N*0.5);
    #neighmax = N*0.5;
    else:
      drivmax = int((N-1)*0.5);
    #neighmax = (N-1)*0.5;
    
    # Get Computation Times
    #compTime = pd.read_csv(str(fileloc+"Times.csv"), sep=',',index_col=None, header=None, dtype=np.float64)
    
    # Dimensions of systems
    sqrtN = np.sqrt(N)
    cbrtN = np.cbrt(N)
    
    for neigh in range(neighmin, neighmax):
      for ndriv in range(drivmin, drivmax):  # Amplitude plot
        fig1 = plt.figure(idx, figsize=(11, 7), dpi=100)
        ax1 = plt.axes()

        #ax1.set_ylim([0.0, 1.0])
        ax1.set_xlabel(r"$t$  /c", fontsize=size, labelpad=8)
        ax1.set_ylabel(r"$\rho(t)$", fontsize=size, labelpad=12)
        ax1.tick_params(axis='both', labelsize=size)
        ax1.yaxis.set_minor_locator(plt.MultipleLocator(0.125))
        ax1.xaxis.set_minor_locator(plt.MultipleLocator(50))
        ax1.tick_params(axis="x", which='major', direction="in", size=6, top=True, width=2, pad=8, length=3)
        ax1.tick_params(axis="y", which='major', direction="in", size=6, right=True, width=2, pad=8, length=3)
        ax1.tick_params(axis="x", which='minor', direction="in", size=3, top=True, width=2, pad=8, length=0.1)
        ax1.tick_params(axis="y", which='minor', direction="in", size=3, right=True, width=2, pad=8, length=0.1)
        #plt.tight_layout()
        if ndriv == 4 or ndriv == 5 or ndriv == 6:
          for k_idx in range(k_idx_min, k_idx_max):
            k = 0.1*k_idx
            for thresh_idx in range(thresh_idx_min, thresh_idx_max):
              thresh = thresh_idx*0.1
              
              # Filenames
              stem = str(N) + "_" + str(neigh) + "_" + str(ndriv) + "_" + str(k) + "00_" + str(thresh) + "00.csv"
              vFile     = "v_"     + stem
              drivFile  = "driv_"  + stem
              couplFile = "coupl_" + stem
              print vFile
            
              v     = pd.read_csv(str(fileloc+"Data/"+vFile    ), sep=',', index_col=None, header=None, dtype=np.float64)
              driv  = pd.read_csv(str(fileloc+"Data/"+drivFile ), sep=',', index_col=None, header=None, dtype=np.float64)
              #coupl = pd.read_csv(str(fileloc+"Data/"+couplFile), sep=',', index_col=None, header=None, dtype=np.float64, skiprows=start_row)
              
              # Get time points, membrane potential and locations of drivers from dataframes
              t         = v.iloc[:,0]
              driv_locs = driv.iloc[:,0]
              
              if any(v.iloc[:,-1].isna()):
                v = v.iloc[:,1:-1]
              else:
                v = v.iloc[:,1:]

              
              Nf=N

              if WITHOUT_DEAD == True:
                for i in range(N):
                  if np.max(v.iloc[:,i])-np.min(v.iloc[:,i]) < 3e-4:
                    v.iloc[:,i] = np.nan
                    Nf-=1

                # Remove all nan columns (dead neurones)
                v = v.iloc[:,~np.all(np.isnan(np.asarray(v)), axis=0)]
              
              # Get the smallest length scale on which neurones in the system die
              #length_ratio      = lengthscale(v, ndriv, driv_locs, N_inv, dimension)
              rho_global, rho_t = OP(vFile, v, Nf, neigh, ndriv, k, thresh)
              #rho_local         = localOP(v, N, neigh, ndriv, k, thresh)
              rho_local = 0
              
              # Initialise average amplitude array
              #v_av = np.zeros(N)

              # Calculate average amplitude of each neurone
              #for i in range(N): v_av[i] = np.sum(v.iloc[:,i])/float(len(t))
              
              # Append data to arrays for phase transition plot
              #pt_neighN.append(neigh/float(N))
              #pt_ndrivs.append(ndriv/float(N))
              #pt_k.append(k)
              #pt_thresh.append(thresh)
              #pt_rhos.append(rho_global)
              
              # ALTER THE NUMBER OF NEIGHBOURS FOR EACH SYSTEM TYPE AND CHANGE CSV NAME FROM C TO MATCH
              
              #delay_plot(v)

              # Send data off for plotting
              fig1 = F1(ax1, fig1, t, rho_t, ndriv)
              #plt.close()
              #fig2 = F2(ax2, fig2, length_ratio, ndriv)
              #plt.close()
              #fig3 = F3(ax3, fig3, t, v)
              #plt.close()

              idx += 1
              #frequencyCounter(idx, size, stem, driv_locs, v, N, neigh, ndriv, k, thresh)

    # Phase transition plot
    #fig5 = scatter3d(pt_neighN, pt_ndrivs, pt_k, pt_rhos, 5, colorsMap = 'jet')
    #fig6 = scatter3d(pt_k, pt_ndrivs, pt_thresh, pt_rhos, 6, colorsMap = 'jet')
    #fig8, fig9, fig10 = PT_2Dcut(pt_k, pt_ndrivs, pt_thresh, pt_rhos, colorsMap = 'jet')

  #fig1.savefig("Graphs/" + stem + "_rho_t.pdf")
  #fig2.savefig("Graphs/" + stem + "_el_death.pdf")
  #fig3.savefig("Graphs/" + stem + "_amplitudeplot.pdf")
  #fig5.savefig("Graphs/" + stem + "_m_D_k.pdf")
  #fig6.savefig("Graphs/2DNONPERIODIC_N=" + str(N) + "_k_D_thresh.pdf")
  plt.show()

  #plt.close()

