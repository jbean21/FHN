"""
Load a csv file using a given filename.
Check if the filename is a time series of x.
Create a dataframe from the csv.
Return the dataframe and the parameters of the system
"""
import pandas as pd
import csv

def FileLoader(filename):
  #create a dataframe from the csv
  df = pd.read_csv(filename)
  #Remove the transient in the system by removing the first 12 "seconds" worth of data
  df = df.drop(range(6000))
  #Read in the system parameters by splitting the filename up by '_'
  filename = filename[57:-4]
  splitt = filename.split("_")
  #filename is in order of N, neighbours, drivers, coupling strength, threshold value
  #use float() method to convert strings into floats
  #print(splitt)
  N = float(splitt[1])
  neighbours = float(splitt[2])
  drivers = float(splitt[3])
  k = float(splitt[4])
  threshold = float(splitt[5]) 
  return [df, [N,neighbours,drivers,k,threshold]]
        
def OPinit():
  with open("/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/OrderParameter.csv", 'w') as f:
    writer = csv.writer(f,lineterminator = '\n')
    writer.writerow(['N','neigh','drivers','k','thresh','p', 'Nf/N'])

def OPread():
  df = pd.read_csv("/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/OrderParameter.csv")
  df.dropna()
  return df

