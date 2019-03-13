from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cmx
import matplotlib as mpl
import pandas as pd
import numpy as np
import csv
from scipy.signal import hilbert
import threading

from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

"""
  import warnings
  warnings.filterwarnings("ignore")
"""

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
  fig = plt.figure(nfig)
  ax = Axes3D(fig)
  ax.scatter(x, y, z, c=scalarMap.to_rgba(cs), s=50)
  scalarMap.set_array(cs)
  if nfig == 5:
    ax.set_xlabel(r"$\frac{\mathrm{neigh}}{N}$")
    ax.set_ylabel(r"$\frac{\mathrm{ndrivs}}{N}$")
    ax.set_zlabel(r"$k$")
  elif nfig == 6:
    ax.set_xlabel(r"$\frac{\mathrm{neigh}}{N}$")
    ax.set_ylabel(r"$\frac{\mathrm{ndrivs}}{N}$")
    ax.set_zlabel(r"$\mathrm{thresh}$")
  
  fig.colorbar(scalarMap)
  
  return

def plot_function(t, v_av, N, neigh, ndriv, k, thresh, length_ratio, rho_t, rho_local, rho_global):
  
  # rho_t vs t
  fig1 = plt.figure(1)
  ax1 = plt.axes()
  ax1.plot(t, rho_t)
  
  # S(t) vs iH(t)
  
  # L_death/L vs ndrivers vs Neigh
  fig2 = plt.figure(2)
  ax2 = fig2.add_subplot(1,1,1, projection='3d')
  ax2.scatter(length_ratio, ndriv, neigh)
  ax2.set_xlabel(r"$w$")
  ax2.set_ylabel(r"$v$")
  ax2.set_zlabel(r"$I_\mathrm{eff}(t)$")
  
  # Amplitude plot
  fig3 = plt.figure(3)
  ax3 = plt.axes()
  x = np.linspace(0, N-1, N)
  ax3.plot(x, v_av)
  
  """
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
  """
  
  """
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
    """
  
  string = str(N) + "_" + str(neigh) + "_" + str(ndriv) + "_" + str(k) + "00_" + str(thresh) + "00_"
  
  #fig1.savefig(string + "rho_t.pdf")
  #fig2.savefig(string + "el_death.pdf")
  #fig3.savefig(string + "amplitudeplot.pdf")
  #fig4.savefig(string + "powerspectrum.pdf")
  
  return


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
  
  for neurone in v: expips.append(find_expip(v, N, neigh, ndriv, k, thresh, neurone))
  
  # Find time varying OP and its time average
  orderParam = np.array([np.abs(np.sum(x)/N) for x in zip(*expips)])
  OPavg = np.mean(orderParam)
  
  # Time varying order parameter
  np.savetxt("OP_" + vFile, orderParam,delimiter = ",\n")
  
  # Write global OP for the system to file
  fields = [N, neigh, ndriv, k, thresh, OPavg]
  with open("OrderParameter.csv", 'a') as f:
    writer = csv.writer(f)
    writer.writerow(fields)
  
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

class myThread (threading.Thread):
  def __init__(self, threadID, N_inv, dimension, N, neigh, ndriv, k, thresh):
    threading.Thread.__init__(self)
    self.threadID  = threadID
    self.N_inv     = N_inv
    self.dimension = dimension
    self.N         = N
    self.neigh     = neigh
    self.ndriv     = ndriv
    self.k         = k
    self.thresh    = thresh

    # Initialise average amplitude array
    self.v_av = np.zeros(self.N)

    self.rho_t        = 0
    self.rho_local    = 0
    self.rho_global   = 0
    self.length_ratio = 0
  def run(self):
    # Arrays for phase transition plot

    print "Thread " + str(self.threadID)

    # Filenames
    stem = str(self.N) + "_" + str(self.neigh) + "_" + str(self.ndriv) + "_" + str(self.k) + "00_" + str(self.thresh) + "00.csv"
    self.vFile     = "v_"     + stem
    self.drivFile  = "driv_"  + stem
    self.couplFile = "coupl_" + stem
    
    # Read data and remove transience by skipping to an appropriate row
    v     = pd.read_csv(str(fileloc+"Data/"+self.vFile    ), sep=',', index_col=None, header=None, dtype=np.float64, skiprows=start_row)
    driv  = pd.read_csv(str(fileloc+"Data/"+self.drivFile ), sep=',', index_col=None, header=None, dtype=np.float64)
    #coupl = pd.read_csv(str(fileloc+"Data/"+self.couplFile), sep=',', index_col=None, header=None, dtype=np.float64, skiprows=start_row)
  
    # Get time points, membrane potential and locations of drivers from dataframes
    self.t         = v.iloc[:,0]
    self.v         = v.iloc[:,1:-1]
    self.driv_locs = driv.iloc[:,0]

    # Run the program
    # Get the smallest length scale on which neurones in the system die
    self.length_ratio = lengthscale(self.v, self.ndriv, self.driv_locs, self.N_inv, self.dimension)
    self.rho_global, self.rho_t = OP(self.vFile, self.v, self.N, self.neigh, self.ndriv, self.k, self.thresh)
    #rho_local         = localOP(v, N, neigh, ndriv, k, thresh)
    self.rho_local = 0

    # Calculate average amplitude of each neurone
    for i in range(self.N): self.v_av[i] = np.sum(self.v.iloc[:,i])/len(self.t)
    
  def return_values(self):
    return self.t, self.v_av, self.N, self.neigh, self.ndriv, self.k, self.thresh, self.length_ratio, self.rho_t, self.rho_local, self.rho_global

######################################################################
""" Main Function """

if __name__ == "__main__":
  # Make plots look nice
  mpl.rcParams.update({'font.size': 20})
  mpl.rcParams.update({'font.family': 'Times New Roman'})
  mpl.rcParams.update({'figure.autolayout': False})

  rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})
  rc('text', usetex = True)

  fileloc = "/media/jay/0FD90FF80FD90FF8/PROJECTDATA/" + "2DPERIODIC/"

  # Initial values
  nrootmin = 3
  nrootmax = 4
  neighmin = 1
  neighmax = 2
  drivmin = 1
  k_idx_min = 0
  k_idx_max = 5
  thresh_idx_min = 0
  thresh_idx_max = 10

  timestep = 1.5e-3
  transient_time = 11

  #------------#
  dimension = 2
  #------------#

  start_row = int(transient_time / timestep) + 1

  # Loop over N
  for nroot in range(nrootmin, nrootmax):
    
    N = 1
    for dim in range(dimension): N *= nroot
    neqn = 2*N
    N_inv = 1/N
    
    # Define a sensible maximum number of drivers
    if N%2==0:
      drivmax  = int(N*0.5)
      #neighmax = N*0.5;
    else:
      drivmax  = int((N-1)*0.5)
      #neighmax = (N-1)*0.5;
    
    # Get Computation Times
    #compTime = pd.read_csv(str(fileloc+"Times.csv"), sep=',',index_col=None, header=None, dtype=np.float64)
    
    # Dimensions of systems
    sqrtN = np.sqrt(N)
    cbrtN = np.cbrt(N)
    
    
    for neigh in range(neighmin, neighmax):
      for ndriv in range(drivmin, drivmax):
        for k_idx in range(k_idx_min, k_idx_max):
          k = 0.1*k_idx

          # Initialise parallelisation
          threads = []
          #threadLock = threading.Lock()

          # Get a thread
          for thresh_idx in range(thresh_idx_min, thresh_idx_max):
            thresh = thresh_idx*0.1
            thread_obj = myThread(thresh_idx, N_inv, dimension, N, neigh, ndriv, k, thresh)
            thread_obj.start()

            # Keep track
            threads.append(thread_obj)
          
          for t in threads:
            t.join()
            t, v_av, N, neigh, ndriv, k, thresh, length_ratio, rho_t, rho_local, rho_global = t.return_values()
            
            # Send data off for plotting
            plot_function(t, v_av, N, neigh, ndriv, k, thresh, length_ratio, rho_t, rho_local, rho_global)
    
            pt_neighN = []
            pt_ndrivs = []
            pt_k      = []
            pt_thresh = []
            pt_rhos   = []
            
            # Append data to arrays for phase transition plot
            pt_neighN.append(neigh/N)
            pt_ndrivs.append(ndriv/N)
            pt_k.append(k)
            pt_thresh.append(thresh)
            pt_rhos.append(rho_global)
          


    # Phase transition plot
    #scatter3d(pt_neighN, pt_ndrivs, pt_k, pt_rhos, 5, colorsMap = 'jet')
    #scatter3d(pt_neighN, pt_ndrivs, pt_thresh, pt_rhos, 6, colorsMap = 'jet')
    
    plt.show()

