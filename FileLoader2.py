"""
Load a csv file using a given filename.
Check if the filename is a time series of x.
Create a dataframe from the csv.
Return the dataframe and the parameters of the system
"""
import pandas as pd
import csv

def FileLoader2(filename):
  #create a dataframe from the csv
  df = pd.read_csv(filename,header=None)
  #Remove the transient in the system by removing the first 12 "seconds" worth of data
  #Read in the system parameters by splitting the filename up by '_'
  filename = filename[63:-4]
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


