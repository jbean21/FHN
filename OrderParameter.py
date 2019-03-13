"""
Calculate the order parameter for a system as a function of time,
record this in a csv, and also record the average order parameter of
each system as a function of the system's input parameters.
"""
from scipy.signal import hilbert
import numpy as np 
import csv


def find_expip(v, N, neigh, ndriv, k, thresh, neurone):
  """ Utility function for OP calculation """

  data = v.iloc[:,neurone]
  mean = np.mean(data)

  #Subtract the mean from data to ensure that Argand diagram centred on 0
  data = np.subtract(data,mean)
  
  #print(data)
  
  #calculates the analytic signal and then exp(i*phi)
  analytic_signal = hilbert(data)
  i_p = np.unwrap(np.angle(analytic_signal))
  

  return list(map(lambda x: np.exp(1j*x),i_p))

def OP(vFile, v, N, neigh, ndriv, k, thresh):
  #filename,df,params  
  """ Calculates global OP """
  expips = []
  #Number of firing neurons
  Nf = N
  #Reduce count of firing neurons
  
  #create
  alive_list = []
  
  WITHOUT_DEAD = True
  
  if WITHOUT_DEAD == True:
    for i in range(1, v.shape[1]-1): 
      if np.max(v.iloc[:,i])-np.min(v.iloc[:,i]) < 5e-3:
          Nf = Nf - 1
      else:
        alive_list.append(i)
  else:
    alive_list=[i for i in range(1,v.shape[1]-1)]
  
  for neurone in alive_list: expips.append(find_expip(v, Nf, neigh, ndriv, k, thresh, neurone))
  
  
  # Find time varying OP and its time average
  
  
  orderParam = np.array([np.abs(np.sum(x)/Nf) for x in zip(*expips)]) #divide by number of firing neurons, not total number of neurons
  OPavg = np.mean(orderParam)
  if len(alive_list) == 1:
    OPavg = 0
    
  life_ratio = Nf/N #write this to a csv, perhaps with the order param in same csv
  
  # Time varying order parameter, again index 17 comes from filepath, this must be changed for different computers
  np.savetxt("/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/OP/OP_" + vFile[59:], orderParam,delimiter = ",\n")
  
  # Write global OP for the system to file
  fields = [N, neigh, ndriv, k, thresh, OPavg, life_ratio]
  with open("/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/OrderParameter.csv", 'a') as f:
    writer = csv.writer(f, lineterminator = "\n")
    writer.writerow(fields)



