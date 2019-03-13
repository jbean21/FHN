import pandas as pd
import numpy as np
import csv

fileloc = "/media/jay/0FD90FF80FD90FF83/PROJECTDATA/" + "2DPERIODIC_2/"

rho = False
Nf = False

#rho = True
Nf = True

parameters = {'N'     : 0,
              'neigh' : 1,
              'D'     : 2,
              'k'     : 3,
              'thresh': 4
              }

incr       = {'N'     : -1,
              'neigh' : 1,
              'D'     : -1,
              'k'     : 10,
              'thresh': 10
              }

# N=0, neigh=1, D=2, k=3, thresh=4
# increment=how many N/how many neigh etc there are
var = 'D'
column = parameters[var]
increment = incr[var]

avgs_over_thresh = []
Ns = []
N_idx = []

File     = "OrderParameter_ordered100_no_nans.csv"
v = pd.read_csv(str(fileloc+File), sep=',', index_col=None, header=None, dtype=np.float64, skiprows=1)
v = np.asarray(v)

i=0
N_idx.append(i)

for i in range(1,v.shape[0]):
  if v[i,column]!=v[i-1,column]:
    N_idx.append(i)

N_idx.append(len(v)-1)

N_idx = list(set(N_idx))

for i in range(len(N_idx)-1):
  avgs_over_thresh.append(sum(v[N_idx[i]:N_idx[i+1],5])/(N_idx[i+1]-N_idx[i]))
  Ns.append(v[N_idx[i],column])
  print i, Ns[i], N_idx[i+1], avgs_over_thresh[i]


array = np.zeros([len(Ns),3])

array[:,0] = Ns
array[:,1] = avgs_over_thresh


if var=='N': 
  np.savetxt("/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/AVERAGES_"+var+"_rhot_no_nans_v2.csv", array, delimiter=',', newline='\n', fmt="%0.3f")

avgs_over_thresh = []
Ns = []
N_idx = []

File     = "OrderParameter_ordered100_no_nans.csv"
v = pd.read_csv(str(fileloc+File), sep=',', index_col=None, header=None, dtype=np.float64, skiprows=1)
v = np.asarray(v)
i=0
N_idx.append(i)

for i in range(1,v.shape[0]):
  if v[i,column]!=v[i-1,column]:
    N_idx.append(i)

N_idx.append(len(v)-1)

N_idx = list(set(N_idx))

for i in range(len(N_idx)-1):
  avgs_over_thresh.append(sum(v[N_idx[i]:N_idx[i+1],6])/(N_idx[i+1]-N_idx[i]))
  Ns.append(v[N_idx[i],column])
  print i, Ns[i], N_idx[i+1], avgs_over_thresh[i]

array[:,2] = avgs_over_thresh


if var == 'D':
  array2 = np.zeros([49,3])
  #np.savetxt("/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/AVERAGES_"+var+"_rhot_no_nans_Nf.csv", array, delimiter=',', newline='\n', fmt="%0.3f")
  for i in range(49):
    rho_mean = 0
    Nf_mean = 0
    counter = 0
    array2[i,0] = i+1
   
    for j in range(array.shape[0]):
      if array[j,0] == i:
        rho_mean += np.sum(array[i,1])
        Nf_mean  += np.sum(array[i,2])
        counter += 1
        
    if counter!=0:
      rho_mean /= counter
      Nf_mean /= counter
      
      array2[i,1] = rho_mean
      array2[i,2] = Nf_mean
  
  np.savetxt("/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/AVERAGES_"+var+"_rhot_no_nans.csv", array2, delimiter=',', newline='\n', fmt="%0.3f")
else:
  array2 = np.zeros([10,3])

  array2[:,0] = array[:increment,0]
  for i in range(increment):
    array2[i,1] = sum(array[i::increment,1])/(array.shape[0]/float(increment))
    array2[i,2] = sum(array[i::increment,2])/(array.shape[0]/float(increment))


  np.savetxt("/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/AVERAGES_"+var+"_rhot_no_nans.csv", array2, delimiter=',', newline='\n', fmt="%0.3f")


"""
idx = 0

for ndriv in range(1,49):
  v_without_drivers = np.zeros([v.shape[0],v.shape[1]])
  for i in range(v.shape[0]):
    if v[i, 2] == ndriv:
      v_without_drivers[idx] = v[i]
      idx += 1
  
  v_without_drivers = v_without_drivers[v_without_drivers[:,0]!=0]
  
  ind = np.lexsort((v_without_drivers[:,3],v_without_drivers[:,2],v_without_drivers[:,1],v_without_drivers[:,0]))

  #v_without_drivers = v_without
  np.savetxt("/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC/ndriv/OrderParameter_ordered_all_"+str(ndriv)+"_drivers.csv", v_without_drivers[ind], delimiter=',', newline='\n', fmt="%0.3f")
"""
