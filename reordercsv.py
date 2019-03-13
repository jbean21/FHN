import pandas as pd
import numpy as np
import csv

fileloc = "/media/jay/0FD90FF80FD90FF83/PROJECTDATA/" + "2DPERIODIC/"
File     = "OrderParameter.csv"
v = pd.read_csv(str(fileloc+File), sep=',', index_col=None, header=None, dtype=np.float64, skiprows=1)
v = np.asarray(v)

ind = np.lexsort((v[:,4],v[:,3],v[:,2],v[:,1],v[:,0]))

np.savetxt("/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/OrderParameter_ordered.csv", v[ind], delimiter=',', newline='\n', fmt="%0.3f")


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
  np.savetxt("/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/ndriv/OrderParameter_ordered_without_dead_"+str(ndriv)+"_drivers.csv", v_without_drivers[ind], delimiter=',', newline='\n', fmt="%0.3f")
